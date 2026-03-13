export AttractorSeedContinueMatch
import ProgressMeter
using Random: MersenneTwister

"""
    AttractorSeedContinueMatch(mapper, matcher = MatchBySSSetDistance(); seeding)

A global continuation method for [`global_continuation`](@ref).
`mapper` is any subtype of [`AttractorMapper`](@ref) which implements
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
are found using the given `mapper` and [`basins_fractions`](@ref).
See the `mapper` documentation and [`AttractorMapper`](@ref)
for details on how this works. Then, from the second parameter onwards the continuation occurs.

### Step 1 - Seeding initial conditions

Initial conditions can be seeded from previously found attractors.
This is controlled by the `seeding` keyword, which must be a function that given
a `StateSpaceSet` (an attractor), it returns an iterator of initial conditions.
By default the first point of an attractor is provided as the only seed.

Seeding can be turned off by providing the dummy function `seeding = A -> []`,
i.e., it always returns an empty iterator and hence no seeds and we skip to step 2.

### Step 2 - Continuing the seeds

The dynamical system referenced by the `mapper` is now set to the new parameter value.
The seeds are run through the `mapper` to converge to attractors at the new parameter value.
Seeding initial conditions close to previous attractors increases the probability
that if an attractor continues to exist in the new parameter, it is found.
Additionally, for some `mappers` this seeding process improves the accuracy as well as
performance of finding attractors, see e.g. discussion in [Datseris2023](@cite).

This seeding works for any `mapper`, regardless of if they can map individual initial conditions
with the `mapper(u0)` syntax! If this syntax isn't supported, steps 2 and 3 are done together.

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
    mapper::A
    matcher::M
    seeding::S
end

const ASCM = AttractorSeedContinueMatch

ASCM(mapper, matcher = MatchBySSSetDistance(); seeding = _default_seeding) =
    ASCM(mapper, matcher, seeding)

# TODO: This is currently not used, and not sure if it has to be.
function _scaled_seeding(attractor::AbstractStateSpaceSet; rng = MersenneTwister(1))
    max_possible_seeds = 6
    seeds = round(Int, log(max_possible_seeds, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, vec(attractor)) for _ in 1:seeds)
end

# This is the one used
function _default_seeding(attractor::AbstractStateSpaceSet)
    return (attractor[1],) # must be iterable
end

function global_continuation(
        ascm::AttractorSeedContinueMatch, pcurve, ics;
        samples_per_parameter = 100, show_progress = true,
    )
    N = samples_per_parameter
    progress = ProgressMeter.Progress(
        length(pcurve);
        desc = "Continuing attractors and basins:", enabled = show_progress
    )
    mapper = ascm.mapper
    prev_attractors = empty(extract_attractors(mapper))
    additional_ics = typeof(current_state(referenced_dynamical_system(mapper)))[]
    # At each parameter `p`, a dictionary mapping attractor ID to fraction is created.
    attractors_cont = Dict[]
    fractions_cont = Dict[]
    ProgressMeter.next!(progress)
    # Continue loop over all remaining parameters
    for p in pcurve
        set_parameters!(referenced_dynamical_system(mapper), p)
        reset_mapper!(mapper)
        # Seed initial conditions from previous attractors.
        # Here we utilize the interal keyword `additional_ics` of `basins_fractions`.
        # The seeding process finds attractors with really small basins, and we need
        # to take this into account when creating the basin fractions, otherwise there
        # could be attractor IDs with 0 fractions in the final output (if the small basin
        # attractors are not found from the random sampling)
        # collect seeds
        empty!(additional_ics)
        for att in values(prev_attractors)
            for u0 in ascm.seeding(att)
                push!(additional_ics, u0)
            end
        end
        # now prepare the initial conditions if per-parameter is requested
        if ics isa PerParameterInitialConditions
            pics = ics.generator(p, N)
        else
            pics = ics
        end
        # and finally call basin fractions; it knows how to do all calculations given the mapper
        # TODO: Would be nice to enable nested progress meters here!
        ret = basins_fractions(mapper, pics; N, additional_ics, show_progress = false)
        fs = pics isa AbstractVector ? ret[1] : ret # if fractions also return labels.
        # deepcopy is important here as attractor container always referrenced
        prev_attractors = deepcopy(extract_attractors(mapper))
        # we don't match attractors here, this happens directly at the end.
        # here we just store the result
        push!(fractions_cont, fs)
        push!(attractors_cont, prev_attractors)
        ProgressMeter.next!(progress)
    end
    rmaps = match_sequentially!(
        attractors_cont, ascm.matcher; pcurve, ds = referenced_dynamical_system(mapper)
    )
    match_sequentially!(fractions_cont, rmaps)
    return fractions_cont, attractors_cont
end
