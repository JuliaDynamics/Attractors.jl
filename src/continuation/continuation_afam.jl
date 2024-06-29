export AttractorsContinueAndMatch
import ProgressMeter
using Random: MersenneTwister

"""
    AttractorsContinueAndMatch(mapper, matcher = MatchBySSSetDistance())

A global_continuation method for [`global_continuation`](@ref).
`mapper` is any subtype of [`AttractorMapper`](@ref) which implements
[`extract_attractors`](@ref).
`matcher` is a configuration of how to match attractor IDs,
and at the moment can only be an instance of [`MatchBySSSetDistance`](@ref).

## Description

This global_continuation method is a generalization of the "RAFM" global_continuation
described in [Datseris2023](@cite). It continues attractors by
"seeding" initial conditions from previously found attractors.
The generalization here is that the method works for any valid `mapper`.

At the first parameter slice of the global_continuation process, attractors and their fractions
are found using the given `mapper`.
At each subsequent parameter slice,
new attractors are found by selecting initial conditions from the previously found
attractors and then running these initial conditions through the `mapper`.
This process is called "seeding".
Seeding initial conditions close to previous attractors increases the probability
that if an attractor continues to exist in the new parameter, it is found.
Seeding is controlled by the `seeding` input. It is a function that given
a state space set (an attractor) it returns an _iterator_ of initial conditions
sampled from the attractor. This iterator generates the seeds.
By default the first point of an attractor is provided as the only seed.

After the special seeded initial conditions are mapped to attractors,
attractor basin fractions are computed by sampling additional initial conditions
using the provided `ics` in [`global_continuation`](@ref)).
I.e., exactly as in [`basins_fractions`](@ref).
Naturally, during this step new attractors may be found, besides those found
using the "seeding from previous attractors".
Once the basins fractions are computed,
the parameter is incremented again and we perform the steps as before.

This process continues for all parameter values. After all parameters are exhausted,
the found attractors (and their fractions) are "matched" to the previous ones.
This means: their _IDs are changed_, so that attractors that are "similar" to those at a
previous parameter get assigned the same ID.
Matching is done by using the provided `matcher` in [`match_sequentially!`](@ref).
If you don't like the final matching output,
you may use a different `matcher` and call [`match_sequentially!`](@ref) again,
without having to recompute the whole global_continuation!

The matching algorithm is a bit involved, so it is best to read the documentation
of `matcher` for how it works in detail.
"""
struct AttractorsContinueAndMatch{A, M, S} <: GlobalContinuationAlgorithm
    mapper::A
    matcher::M
    seeding::S
end

const ACAM = AttractorsContinueAndMatch

ACAM(mapper, matcher) = ACAM(mapper, matcher, _default_seeding_process)

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

function global_continuation(acam::AttractorsContinueAndMatch, prange, pidx, ics;
        samples_per_parameter = 100, show_progress = true,
    )
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
        fs, = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
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
        seeded_fs = Dict{Int, Int}()
        for att in values(prev_attractors)
            for u0 in acam.seeding(att)
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
                additional_fs = seeded_fs, show_progress = false, N = samples_per_parameter
            )
        else
            fs, = basins_fractions(mapper, ics;
                additional_fs = seeded_fs, show_progress = false, N = samples_per_parameter
            )
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
