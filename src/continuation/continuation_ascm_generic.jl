export AttractorSeedContinueMatch
import ProgressMeter
using Random: MersenneTwister

"""
    AttractorSeedContinueMatch(mapper; matcher = MatchBySSSetDistance(), seeding)

A global continuation method for [`global_continuation`](@ref).
`mapper` is any subtype of [`AttractorMapper`](@ref) which implements
[`extract_attractors`](@ref), i.e., it finds the actual attractors.
`matcher` is a configuration of how to match attractor IDs,
and at the moment can only be an instance of [`MatchBySSSetDistance`](@ref).

## Description

This is a general global continuation method is based on a 4-step process:

1. Seed initial conditions from previously found attractors
2. Propagate those forwards to "continue" previous attractors
3. Estimate basin fractions and potentially find new attractors
4. Match attractors

### Step 0 - Finding attractors

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
i.e., it always returns an empty iterator and hence no seeds.

### Step 2 - Continuing the seeds

The seeds are run through the `mapper` to converge to attractors at the new parameter value.
Seeding initial conditions close to previous attractors increases the probability
that if an attractor continues to exist in the new parameter, it is found.
Additionally, for some `mappers` this seeding process improves the accuracy as well as
performance of finding attractors, see e.g., discussion in [Datseris2023](@cite).

This seeding works for any `mapper`, regardless of if they can map individual initial conditions
with the `mapper(u0)` syntax! If this syntax isn't supported, steps 2 and 3 are done together.

### Step 3 - Estimate basins fractions

After the special seeded initial conditions are mapped to attractors,
attractor basin fractions are computed by sampling additional initial conditions
using the provided `ics` in [`global_continuation`](@ref).
I.e., exactly as in [`basins_fractions`](@ref).
Naturally, during this step new attractors may be found, besides those found
using the "seeding from previous attractors".
Once the basins fractions are computed,
the parameter is incremented again and we perform steps 1-3 as before.

### Step 4 - Matching

Steps 1-3 continue for all parameter values. After all parameters are exhausted,
the found attractors (and their fractions) are "matched" to the previous ones
starting from the second parameter slice onwards.
This means: their _IDs are changed_, so that attractors that are "similar" to those at a
previous parameter get assigned the same ID.

At the moment matching is done by distance of state space sets, as described
in [`MatchBySSSetDistance`](@ref).
The matching itself happes via a trivial call to [`match_sequentially!`](@ref)
at the end of the `global_continuation` call.

If you don't like the final matching output,
you may use a different `matcher` and call [`match_sequentially!`](@ref) again,
without having to recompute the whole global continuation!

### Further notes

This global continuation method is a generalization of the "RAFM" continuation
described in [Datseris2023](@cite).
"""
struct AttractorSeedContinueMatch{A, M, S} <: GlobalContinuationAlgorithm
    mapper::A
    matcher::M
    seeding::S
end

const ASCM = AttractorSeedContinueMatch

ASCM(mapper; matcher = MatchBySSSetDistance(), seeding = _default_seeding_process) =
ASCM(mapper, matcher, seeding)

# TODO: This is currently not used, and not sure if it has to be.
function _default_seeding_process_10(attractor::AbstractStateSpaceSet; rng = MersenneTwister(1))
    max_possible_seeds = 6
    seeds = round(Int, log(max_possible_seeds, length(attractor)))
    seeds = clamp(seeds, 1, max_possible_seeds)
    return (rand(rng, vec(attractor)) for _ in 1:seeds)
end

# This is the one used
function _default_seeding_process(attractor::AbstractStateSpaceSet)
    return (attractor[1],) # must be iterable
end

function global_continuation(acam::AttractorSeedContinueMatch, prange, pidx, ics;
        samples_per_parameter = 100, show_progress = true,
    )
    N = samples_per_parameter
    progress = ProgressMeter.Progress(length(prange);
        desc = "Continuing attractors and basins:", enabled=show_progress
    )
    mapper = acam.mapper
    reset_mapper!(mapper)
    # first parameter is run in isolation, as it has no prior to seed from
    set_parameter!(referenced_dynamical_system(mapper), pidx, prange[1])
    if ics isa Function
        fs = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    else # we ignore labels in this continuation algorithm
        fs, = basins_fractions(mapper, ics; show_progress = false)
    end

    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_cont = [fs]
    # The attractors are also stored (and are the primary output)
    prev_attractors = deepcopy(extract_attractors(mapper))
    attractors_cont = [prev_attractors]
    ProgressMeter.next!(progress; showvalues = [("previous parameter", prange[1]),])
    # Continue loop over all remaining parameters
    for p in prange[2:end]
        set_parameter!(referenced_dynamical_system(mapper), pidx, p)
        reset_mapper!(mapper)
        # Seed initial conditions from previous attractors
        # Notice that one of the things that happens here is some attractors have
        # really small basins. We find them with the seeding process here, but the
        # subsequent random sampling in `basins_fractions` doesn't. This leads to
        # having keys in `mapper.bsn_nfo.attractors` that do not exist in the computed
        # fractions. The fix is easy: we add the initial conditions mapped from
        # seeding to the fractions using an internal argument.
        fs = if allows_mapper_u0(mapper)
            seed_attractors_to_fractions_individual(mapper, prev_attractors, ics, N, acam.seeding)
        else
            seed_attractors_to_fractions_grouped(mapper, prev_attractors, ics, N, acam.seeding)
        end
        # We do not match attractors here; the matching is independent step done at the end
        current_attractors = deepcopy(extract_attractors(mapper))
        push!(fractions_cont, fs)
        push!(attractors_cont, current_attractors)
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    # Match attractors (and basins)
    rmaps = match_sequentially!(attractors_cont, acam.matcher)
    match_sequentially!(fractions_cont, rmaps)
    return fractions_cont, attractors_cont
end

function seed_attractors_to_fractions_individual(mapper, prev_attractors, ics, N, seeding)
    # actual seeding
    for att in values(prev_attractors)
        for u0 in seeding(att)
            # We map the initial condition to an attractor, but we don't care
            # about which attractor we go to. This is just so that the internal
            # array of `AttractorsViaRecurrences` registers the attractors
            label = mapper(u0; show_progress = false)
            seeded_fs[label] = get(seeded_fs, label, 0) + 1
        end
    end
    # Now perform basin fractions estimation as normal, utilizing found attractors
    # (the function comes from attractor_mapping.jl)
    if ics isa Function
        fs = basins_fractions(mapper, ics;
            additional_fs = seeded_fs, show_progress = false, N
        )
    else
        fs, = basins_fractions(mapper, ics;
            additional_fs = seeded_fs, show_progress = false
        )
    end
    return seeded_fs
end

function seed_attractors_to_fractions_grouped(mapper, prev_attractors, ics, N, seeding)
    # what makes this version different is that we can't just use `mapper(u0)`,
    # so we need to store the seeded initial conditions and then combine them with
    # the the ones generated from `ics`.
    u0s = typeof(current_state(referenced_dynamical_system(mapper)))[]
    # collect seeds
    for att in values(prev_attractors)
        for u0 in seeding(att)
            push!(u0s, u0)
        end
    end
    # now combine these with the rest of the initial conditions
    if ics isa Function
        for _ in 1:N
            push!(u0s, copy(ics()))
        end
    else
        append!(u0s, vec(ics))
    end
    # with these extra u0s we now perform fractions estimation as normal
    fs, = basins_fractions(mapper, StateSpaceSet(u0s); show_progress = false)
    return fs
end

allows_mapper_u0(::AttractorMapper) = true
function allows_mapper_u0(mapper::AttractorsViaFeaturizing)
    if mapper.group_config isa GroupViaClustering
        return false
    elseif mapper.group_config isa GroupViaPairwiseComparison
        return false
    else
        return true
    end
end