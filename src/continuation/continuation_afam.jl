"""
    AttractorsContinueAndMatch(mapper, matcher [, seeding])

A continuation method for [`continuation`](@ref).
`mapper` is any subtype of [`AttractorMapper`](@ref) which implements
[`extract_attractors`](@ref), while `matcher` is any subtype of [`SSSetMatcher`](@ref).

## Description

This global continuation method has a simple basis: it continues attractors by
"seeding" initial conditions from previously found attractors.

At the first parameter slice of the continuation process, attractors and their fractions
are found using the given `mapper`.
At each subsequent parameter slice,
new attractors are found by seeding initial conditions from the previously found
attractors and then running these initial conditions through the `mapper`.
Seeding initial conditions close to previous attractors increases the probability
that if an attractor continues to exist in the new parameter, it is found.
Seeding is controlled by the `seeding` input. It is a function that given
a state space set (an attractor) it returns an _iterator_ of initial conditions
sampled from the attractor. This generates the seeds.

After the special seeded initial conditions are mapped to attractors,
attractor basin fractions are computed by sampling additional initial conditions
using the provided `ics` in [`continuation`](@ref)).
I.e., exactly as in [`basins_fractions`](@ref).
Naturally, during this step new attractors may be found, besides those found
using the "seeding from previous attractors".
Once the basins fractions are computed,
the parameter is incremented again and we perform the steps as before.

This process continues for all parameter values. After all parameters are exhausted,
the found attractors (and their fractions) are "matched" to the previous ones.
This means: their _IDs are changed_, so that attractors that are "similar" to those at a
previous parameter get assigned the same ID.
Matching is done by the provided `matcher`.
In code, matching is a rather trivial call to [`match_continuation!`](@ref)
at the end of the `continuation` function.
If you don't like the final matching output,
you may use arbitrarily different `matcher` in [`match_continuation!`](@ref)
without having to recompute the whole continuation.

The matching algorithms are rather sophisticated,
so how matching works is described in the docstrings of each `matcher`.
"""
struct AttractorsContinueAndMatch{A, M, S} <: AttractorsBasinsContinuation
    mapper::A
    matcher::M
    seeding::S
end

function continuation(acam::AttractorsContinueAndMatch, prange, pidx, ics;
        samples_per_parameter = 100, show_progress = true,
    )
    progress = ProgressMeter.Progress(length(prange);
        desc = "Continuing attractors and basins:", enabled=show_progress
    )
    mapper = acam.mapper
    reset_mapper!(mapper)
    # first parameter is run in isolation, as it has no prior to seed from
    set_parameter!(referenced_dynamical_system(mapper), pidx, prange[1])
    fs = basins_fractions(mapper, ics; show_progress = false, N = samples_per_parameter)
    # At each parmaeter `p`, a dictionary mapping attractor ID to fraction is created.
    fractions_curves = [fs]
    # The attractors are also stored (and are the primary output)
    prev_attractors = deepcopy(extract_attractors(mapper))
    attractors_info = [prev_attractors]
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
        fs = basins_fractions(mapper, ics;
            additional_fs = seeded_fs, show_progress = false, N = samples_per_parameter
        )
        # We do not match attractors here; the matching is independent step done at the end
        current_attractors = deepcopy(extract_attractors(mapper))
        push!(fractions_curves, fs)
        push!(attractors_info, current_attractors)
        overwrite_dict!(prev_attractors, current_attractors)
        ProgressMeter.next!(progress; showvalues = [("previous parameter", p),])
    end
    # Match attractors (and basins)
    rmaps = match_continuation!(attractors, matcher)
    match_continuation!(rmaps, fractions_curves)
    return fractions_curves, attractors_info
end
