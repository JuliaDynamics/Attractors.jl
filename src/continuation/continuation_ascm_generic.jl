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

This is a general/composable global continuation method based on a 4-step process:

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
        ascm::AttractorSeedContinueMatch, pcurve::Vector, icsampler::InitialConditionSampler;
        show_progress = true,
    )
    progress = ProgressMeter.Progress(
        length(pcurve);
        desc = "Global continuation:", PMKWARGS..., enabled = show_progress
    )
    bmap = ascm.bmap
    ds = referenced_dynamical_system(bmap)
    additional_ics = typeof(current_state(ds))[]
    prev_attractors = empty(extract_attractors(bmap))
    prev_p = nothing
    # The matching happens here, one parameter at a time, and its running state is kept
    # across the loop. This is the same chain `match_sequentially!` would rebuild at the
    # end, but the loop cannot wait for it: `update_sampler!` is handed labels, and those
    # have to name the same attractor at every parameter for the sampler to compare them.
    # Doing it once, here, is what keeps a single set of IDs in play.
    use_vanished = _use_vanished(ascm.matcher)
    latest_ghosts = empty(extract_attractors(bmap))
    next_id = 1
    # At each parameter `p`, a dictionary mapping attractor ID to fraction is created.
    attractors_cont = Dict[]
    fractions_cont = Dict{Int, Float64}[]
    # Continue loop over all remaining parameters
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
        rmap = Dict{Int, Int}()
        counts = Dict{Int, Float64}()
        local attractors
        while true
            # prepare the initial conditions
            pics = generate_ics(icsampler, p)
            # and finally call basin fractions; it knows how to do all calculations given the bmap
            fs, labels = basins_fractions_labels(bmap, pics; additional_ics, show_progress, offset = 2)
            # `fs` was normalized by the sampled ics, which is what `labels` counts, plus
            # the seeded ones
            n_round = length(labels) + length(additional_ics)
            for (k, v) in fs
                counts[k] = get(counts, k, 0.0) + v * n_round
            end
            empty!(additional_ics)
            attractors = deepcopy(extract_attractors(bmap))
            next_id = _extend_matching!(
                rmap, attractors, use_vanished ? latest_ghosts : prev_attractors,
                ascm.matcher, next_id; ds, p, pprev = prev_p,
            )
            replace!(labels, rmap...)
            update_sampler!(icsampler, labels)
            # Now, check if the sampler requires us to re-sample at the current parameter
            resampling_required(icsampler) || break
        end
        # `attractors` keeps the IDs the basin map issued; `rmap` is applied to a copy
        matched_attractors = copy(attractors)
        swap_dict_keys!(matched_attractors, rmap)
        # `counts` was pooled in those same raw IDs, so it is relabelled here, once
        swap_dict_keys!(counts, rmap)
        # All the computations are done, and now we just store the result(s). The
        # attractors are matched already, so there is no matching pass at the end.
        # The fractions go through the sampler, because it is the only one that knows how
        # densely it covered each part of the region; see `weighted_fractions`.
        push!(fractions_cont, weighted_fractions(icsampler, counts))
        push!(attractors_cont, matched_attractors)
        if use_vanished
            for (k, A) in matched_attractors
                latest_ghosts[k] = A
            end
        end
        prev_attractors, prev_p = matched_attractors, p
        # update progress bar
        showvalues = i < length(pcurve) ? [("pcurve index", i + 1)] : []
        ProgressMeter.next!(progress; showvalues)
    end
    return fractions_cont, attractors_cont
end


"""
    _extend_matching!(rmap, attractors, prev, matcher, next_id; kw...) → next_id

Give an ID in `prev`'s space to every attractor of `attractors` that does not have one in
`rmap` yet, leaving the entries `rmap` already holds alone. Returns `next_id`, an ID that
no attractor has been given yet.

Attractors found in a later round of the same parameter are matched only against the
previous attractors that nothing has claimed yet, so extending the map can never
invalidate the labels already handed out.
"""
function _extend_matching!(rmap, attractors, prev, matcher, next_id; kw...)
    fresh = Dict(k => A for (k, A) in attractors if !haskey(rmap, k))
    isempty(fresh) && return next_id
    available = Dict(k => A for (k, A) in prev if k ∉ values(rmap))
    for ids in (keys(attractors), keys(prev), values(rmap))
        isempty(ids) || (next_id = max(next_id, maximum(ids) + 1))
    end
    if isempty(available)
        isempty(prev) && return next_id
        for k in keys(fresh)
            rmap[k] = next_id
            next_id += 1
        end
    else
        # `matching_map` names all of `fresh`, the unmatched ones from `next_id` upwards;
        # how many it took is read back off `rmap` the next time around.
        merge!(rmap, matching_map(fresh, available, matcher; next_id, kw...))
    end
    return next_id
end
