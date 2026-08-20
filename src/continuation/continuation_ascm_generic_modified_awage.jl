export AttractorSeedContinueMatch
import ProgressMeter
using Random: MersenneTwister

"""
    AttractorSeedContinueMatch(bmap, matcher = MatchBySSSetDistance(); seeding)

A global continuation method for [`global_continuation`](@ref).
`bmap` is any subtype of [`BasinMap`](@ref) which implements
[`extract_attractors`](@ref), i.e., it finds the actual attractors.
`matcher` is a configuration of how to match attractor IDs, see [`IDMatcher`](@ref)
for more options.

## Description

This is a general and composable global continuation method outlined in our article
[Datseris2026](@cite) and based on a 4-step process:

1. Seed initial conditions from previously found attractors
2. Propagate those forwards to "continue" previous attractors
3. Estimate basin fractions and potentially find new attractors
4. Match attractors

### Step 0 - Finding initial attractors

At the first parameter slice of the global continuation process, attractors and their fractions
are found using the given `bmap` and [`basins_fractions`](@ref).
See the `bmap` documentation and [`BasinMap`](@ref)
for details on how this works. Then, from the second parameter onwards the continuation occurs.

### Step 1 - Seeding initial conditions

Initial conditions can be seeded from previously found attractors.
This is controlled by the `seeding` keyword, which must be a function that given
a `StateSpaceSet` (an attractor), it returns an iterator of initial conditions.
By default the first point of an attractor is provided as the only seed.

Seeding can be turned off by providing the dummy function `seeding = A -> []`,
i.e., it always returns an empty iterator and hence no seeds and we skip to step 2.

### Step 2 - Continuing the seeds

The dynamical system referenced by the `bmap` is now set to the new parameter value.
The seeds are run through the `bmap` to converge to attractors at the new parameter value.
Seeding initial conditions close to previous attractors increases the probability
that if an attractor continues to exist in the new parameter, it is found.
Additionally, for some `mappers` this seeding process improves the accuracy as well as
performance of finding attractors, see e.g. discussion in [Datseris2023](@cite).

This seeding works for any `bmap`, regardless of if they can map individual initial conditions
with the `bmap(u0)` syntax! If this syntax isn't supported, steps 2 and 3 are done together.

### Step 3 - Estimate basins fractions

After the special seeded initial conditions are mapped to attractors,
attractor basin fractions are computed by sampling additional initial conditions
using the provided `ics` in [`global_continuation`](@ref).
I.e., exactly as in [`basins_fractions`](@ref).
Naturally, during this step new attractors may be found, besides those found
using the "seeding from previous attractors".

### Step 4 - Matching

Normally the ID an attractor gets assigned is somewhat a random integer.
Therefore, to ensure a logical output of the global continuation process,
attractors need to be "matched".
This means: attractor and fractions must have their _IDs changed_,
so that attractors that are "similar" to those at a
previous parameter get assigned the same ID.

What is "similar enough" is controlled by the `matcher` input.
The default `matcher` [`MatchBySSSetDistance`](@ref) matches
sets which have small distance in state space.
The matching algorithm itself can be quite involved,
so read the documentation of the `matcher` for how matching works.

A note on matching: the [`MatchBySSSetDistance`](@ref) can also be used
after the continuation is completed, as it only requires as input
the state space sets (attractors), without caring at which parameter each attractor
exists at. If you don't like the final matching output,
you may use a different instance of [`MatchBySSSetDistance`](@ref)
and call [`match_sequentially!`](@ref) again on the output,
without having to recompute the whole global continuation!

### Step 5 - Finish

After matching the parameter(s) is incremented.
Steps 1-4 repeat until all parameter values are exhausted.

### Further note

This global continuation method is a generalization of the "RAFM" continuation
described in [Datseris2023](@cite). This continuation method is still exported
as [`RecurrencesFindAndMatch`](@ref).
"""
struct AttractorSeedContinueMatch{A, M, S} <: GlobalContinuationAlgorithm
    bmap::A
    matcher::M
    seeding::S
end

const ASCM = AttractorSeedContinueMatch

ASCM(bmap, matcher = MatchBySSSetDistance(); seeding = _default_seeding) =
    ASCM(bmap, matcher, seeding)

# TODO: This is currently not used, and not sure if it has to be.
function _scaled_seeding(attractor::AbstractStateSpaceSet, rng = MersenneTwister(1))
    max_possible_seeds = 9 # ≈ log(10,000)
    seeds = round(Int, log(length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, attractor) for _ in 1:seeds)
end

# This is the one used
function _default_seeding(attractor::AbstractStateSpaceSet)
    return (attractor[1],) # must be iterable
end

function global_continuation(
        ascm::AttractorSeedContinueMatch, pcurve::Vector, icsampler::InitialConditionsSampler;
        show_progress = true,
    )
    # Setup generic variables:
    progress = ProgressMeter.Progress(
        length(pcurve);
        desc = "Global continuation:", PMKWARGS..., enabled = show_progress
    )
    bmap = ascm.bmap
    ds = referenced_dynamical_system(bmap)
    additional_ics = typeof(current_state(ds))[]
    # Setup matching variables:
    # Matching happens inside the loop, one parameter at a time; this is because
    # the `update_sampler!` function needs to have matched labels already.
    # Each matching step is essentially the inner loop of `match_sequentially!`.
    prev_p = first(pcurve)
    use_vanished = _use_vanished(ascm.matcher)
    latest_ghosts = empty(extract_attractors(bmap))
    prev_attractors = empty(extract_attractors(bmap))
    # Setup output containers:
    attractors_cont = typeof(prev_attractors)[]
    fractions_cont = Dict{Int, Float64}[]
    quantifiers_cont = []
    other_cont = Dict{String, Any}()
    # Loop over parameters
    for (i, p) in enumerate(pcurve)
        set_parameters!(ds, p)
        reset_mapper!(bmap)
        # Seed initial conditions from previous attractors.
        empty!(additional_ics)
        for att in values(prev_attractors)
            for u0 in ascm.seeding(att)
                push!(additional_ics, u0)
            end
        end
        counts = Dict{Int, Float64}()
        local matched_attractors, rmap
        while true
            # call basin fractions; it knows how to do all calculations given the bmap,
            # including generating the initial conditions from the sampler
            fs, labels = basins_fractions_labels(
                bmap, icsampler; params = p, additional_ics, show_progress, offset = 2
            )
            # `fs` was normalized by the sampled ics, which is what `labels` counts, plus
            # the seeded ones
            n_round = length(labels) + length(additional_ics)
            for (k, v) in fs
                counts[k] = get(counts, k, 0.0) + v * n_round
            end
            empty!(additional_ics)
            # one step of `match_sequentially!`: the attractors just found, against the
            # previous parameter (or against the ghosts, if the matcher wants them)
            matched_attractors = deepcopy(extract_attractors(bmap))
            prev = use_vanished ? latest_ghosts : prev_attractors
            rmap = only(match_sequentially!(
                [prev, matched_attractors], ascm.matcher;
                retract_keys = false, ds, pcurve = [prev_p, p],
               ))
            replace!(labels, rmap...)
            update_sampler!(icsampler, labels)
            # Now, check if the sampler requires us to re-sample at the current parameter
            resampling_required(icsampler) || break
        end
        # `counts` was pooled in the IDs the basin map issued, so it is relabelled here
        swap_dict_keys!(counts, rmap)
        # All the computations are done, and now we just store the result(s). The
        # attractors are matched already, so there is no matching pass at the end.
        # The fractions go through the sampler, because it is the only one that knows how
        # densely it covered each part of the region; see `weighted_fractions`.
        push!(fractions_cont, weighted_fractions(icsampler, counts))
        push!(attractors_cont, matched_attractors)
        # the quantifiers of the accumulator are also keyed in the IDs the basin map
        # issued, so they are relabelled here as well
        if bmap isa StabilityQuantifiersAccumulator
            quantifiers = finalize_accumulator(bmap)
            for dict in values(quantifiers)
                swap_dict_keys!(dict, rmap)
            end
            push!(quantifiers_cont, quantifiers)
        end
        if use_vanished
            for (k, A) in matched_attractors
                latest_ghosts[k] = A
            end
        end
        prev_attractors, prev_p = matched_attractors, p
        # any extras that need to be updated can be done so here:
        add_extra_continuation_info!(other_cont, icsampler)
        # update progress bar
        showvalues = i < length(pcurve) ? [("pcurve index", i + 1)] : []
        ProgressMeter.next!(progress; showvalues)
    end
    # everything has been matched already inside the loop, so all that is left is to
    # transform the tracked quantities to the agreed output
    fractions, quantifiers = generate_continuation_output(bmap, fractions_cont, quantifiers_cont)
    out = GlobalContinuationOutput(attractors_cont, fractions, quantifiers, other_cont, pcurve)
    return out
end

# This function has a generic form that just forwards the sampled fractions, and a more
# technical form that collects various quantifiers, taken care off by the accumulator
function generate_continuation_output(bmap, fractions_cont, quantifiers_cont)
    return fractions_cont, Dict{String, Vector}()
end
