using BlackBoxOptim: bboptimize, best_candidate
using Random: GLOBAL_RNG
export minimal_fatal_shock, MFSBruteForce, MFSBlackBoxOptim
export excitability_threshold

"""
    minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, algorithm; kw...)

Return the _minimal fatal shock_ `mfs` (also known as _excitability threshold_)
for the initial point `u0` according to the
specified `algorithm` given a `mapper` that satisfies the `id = mapper(u0)` interface
(see [`AttractorMapper`](@ref) if you are not sure which mappers do that).
The `mapper` contains a reference to a [`DynamicalSystem`](@ref).
The options for `algorithm` are: [`MFSBruteForce`](@ref) or [`MFSBlackBoxOptim`](@ref).
For high dimensional systems [`MFSBlackBoxOptim`](@ref) is likely more accurate.

The `search_area` dictates the state space range for the search of the `mfs`.
It can be a 2-tuple of (min, max) values, in which case the same values are used
for each dimension of the system in `mapper`. Otherwise, it can be a vector of 2-tuples,
each for each dimension of the system. The search area is defined w.r.t. to `u0`
(i.e., it is the search area for perturbations of `u0`).

An alias to `minimal_fata_shock` is `excitability_threshold`.

## Keyword arguments

- `metric = Euclidean()`: a metric function that gives the norm of a perturbation vector.
  This keyword is ignored for the [`MFSBruteForce`](@ref) algorithm.
- `target_id = nothing`: when not `nothing`, it should be an integer or a vector of
  integers corresponding to target attractor label(s).
  Then, the MFS is estimated based only on perturbations that lead to the target
  attractor(s).

## Description

The minimal fatal shock is defined as the smallest-norm perturbation of the initial
point `u0` that will lead it a different basin of attraction. It is inspired by the paper
"Minimal fatal shocks in multistable complex networks" [Halekotte2020](@cite),
however the implementation here is generic: it works for _any_ dynamical system.

The _excitability threshold_ is a concept nearly identical, however, instead of looking
for a perturbation that simply brings us out of the basin, we look for the smallest
perturbation that brings us into specified basin(s). This is enabled via the keyword
`target_id`.
"""
function minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, algorithm;
        metric = Euclidean(), target_id = nothing
    )
    id_u0 = mapper(u0)
    # generate a function that returns `true` for ids that that are in the target basin
    idchecker = id_check_function(id_u0, target_id)
    dim = dimension(mapper.ds)
    if typeof(search_area) <: Tuple{Any,Any}
        search_area  = [search_area for _ in 1:dim]
    elseif length(search_area) == 1
        search_area = [search_area[1] for _ in 1:dim]
    elseif length(search_area) != dim
        error("Input search area does not match the dimension of the system")
    end
    return _mfs(algorithm, mapper, u0, search_area, idchecker, metric)
end
const excitability_threshold = minimal_fatal_shock

id_check_function(id::Int, ::Nothing) = i -> i != id
function id_check_function(id::Int, ti::Int)
    id == ti && error("target id and attractor id of u0 are the same.")
    return i -> i == ti
end
function id_check_function(id::Int, ti::AbstractVector{Int})
    id ∈ ti && error("target id and attractor id of u0 are the same.")
    return i -> i ∈ ti
end


"""
    MFSBruteForce(; kwargs...)

The brute force randomized search algorithm used in [`minimal_fatal_shock`](@ref).

It consists of
two steps: random initialization and sphere radius reduction. On the first step,
the algorithm generates random perturbations within the search area and records
the perturbation that leads to a different basin but with the smallest magnitude.
With this
obtained perturbation it proceeds to the second step. On the second step, the algorithm
generates random perturbations on the surface of the hypersphere with radius equal to the
norm of the perturbation found in the first step.
It reduces the radius of the hypersphere and continues searching for the better result
with a smaller radius. Each time a better result is found, the radius is reduced further.

The algorithm records the perturbation with smallest radius that leads to a different basin.

## Keyword arguments

- `initial_iterations = 10000`: number of random perturbations to try in the first step of the
  algorithm.
- `sphere_iterations = 10000`: number of steps while initializing random points on hypersphere and
  decreasing its radius.
- `sphere_decrease_factor = 0.999` factor by which the radius of the hypersphere is decreased
  (at each step the radius is multiplied by this number). Number closer to 1 means
  more refined accuracy
"""
Base.@kwdef struct MFSBruteForce
    initial_iterations::Int64 = 10000
    sphere_iterations::Int64 = 10000
    sphere_decrease_factor::Float64 = 0.999
end

function _mfs(algorithm::MFSBruteForce, mapper, u0, search_area, idchecker, metric)
    algorithm.sphere_decrease_factor ≥ 1 && error("Sphere decrease factor cannot be ≥ 1.")
    dim = dimension(mapper.ds)
    best_shock, best_dist = crude_initial_radius(
        mapper, u0, search_area, idchecker, metric, algorithm.initial_iterations
    )
    best_shock, best_dist = mfs_brute_force(
        mapper, u0, best_shock, best_dist, dim, idchecker, metric,
        algorithm.sphere_iterations, algorithm.sphere_decrease_factor
    )
    return best_shock
end

"""
This function generates a random perturbation of the initial point `u0` within
specified "search_area" and checks if it is in the same basin of attraction.
It does so by generating a random vector of length dim and then adding it to u0.
If the perturbation is not in the same basin of attraction, it calculates the norm
of the perturbation and compares it to the best perturbation found so far.
If the norm is smaller, it updates the best perturbation found so far.
It repeats this process total_iterations times and returns the best perturbation found.
"""
function crude_initial_radius(mapper::AttractorMapper, u0, search_area, idchecker, metric, total_iterations)
    best_dist = Inf
    region = StateSpaceSets.HRectangle([s[1] for s in search_area], [s[2] for s in search_area])
    generator, _ = statespace_sampler(region)
    best_shock = copy(generator())
    shock = copy(best_shock)

    for _ in 1:total_iterations
        perturbation = generator()
        @. shock = perturbation + u0
        new_id = mapper(shock)
        if idchecker(new_id)
            dist = metric(perturbation)
            if dist < best_dist
                best_dist = dist
                best_shock .= perturbation
            end
        end
    end

    return best_shock, best_dist
end


"""
This function works on the results obtained by `crude_initial_radius`. It starts from
the best shock found so far and tries to find a better one by continuously reducing
the radius of the sphere on the surface of which it generates random perturbations.
If perturbation with the same basin of attraction is found, it updates the best shock found
so far and reduces the radius of the sphere. It repeats this process total_iterations times
and returns the best perturbation found.
"""
function mfs_brute_force(mapper::AttractorMapper, u0,
        best_shock, best_dist, dim, idchecker, metric,
        total_iterations, sphere_decrease_factor
    )

    temp_dist = best_dist*sphere_decrease_factor
    region = HSphereSurface(temp_dist, dim)
    generator, = statespace_sampler(region)
    i = 0
    new_shock = zeros(dim)
    while i < total_iterations
        perturbation = generator()
        @. new_shock = perturbation + u0
        new_id = mapper(new_shock)
        # if perturbation leading to another basin:
        if idchecker(new_id)
            # record best
            best_shock .= perturbation
            best_dist = temp_dist
            # update radius
            temp_dist = temp_dist*sphere_decrease_factor
            region = HSphereSurface(temp_dist, dim)
            generator, = statespace_sampler(region)
        end
        i += 1
    end
    return best_shock, best_dist
end


"""
    MFSBlackBoxOptim(; kwargs...)

The black box derivative-free optimization algorithm used in [`minimal_fatal_shock`](@ref).

## Keyword arguments

- `guess = nothing`: a initial guess for the minimal fatal shock given to the
  optimization algorithm. If not `nothing`, `random_algo` below is ignored.
- `max_steps = 10000`: maximum number of steps for the optimization algorithm.
- `penalty = 1000.0`: penalty value for the objective function for perturbations that do
  not lead to a different basin of attraction. This value is added to the norm of the
  perturbation and its value should be much larger than the typical sizes of the basins of
  attraction.
- `print_info`: boolean value, if true, the optimization algorithm will print information on
  the evaluation steps of objective function, `default = false`.
- `random_algo = MFSBruteForce(100, 100, 0.99)`: an instance of [`MFSBruteForce`](@ref)
  that can be used to provide an initial guess.
- `bbkwargs = NamedTuple()`: additional keyword arguments propagated to
  `BlackBoxOptim.bboptimize` for selecting solver, accuracy, and more.

## Description

The algorithm uses BlackBoxOptim.jl and a penalized objective function to minimize.
y function used as a constraint function.
So, if we hit another basin during the search we encourage the algorithm otherwise we
punish it with some penalty. The function to minimize is (besides some details):
```julia
function mfs_objective(perturbation, u0, mapper, penalty)
    dist = norm(perturbation)
    if mapper(u0 + perturbation) == mapper(u0)
        # penalize if we stay in the same basin:
        return dist + penalty
    else
        return dist
    end
end
```
Using an initial guess can be beneficial to both performance and accuracy,
which is why the output of a crude [`MFSBruteForce`](@ref) is used to provide a guess.
This can be disabled by either passing a `guess` vector explicitly or
by giving `nothing` as `random_algo`.
"""
Base.@kwdef struct MFSBlackBoxOptim{G, RA, BB}
    guess::G = nothing
    max_steps::Int64 = 10_000
    penalty::Float64 = 0.999
    print_info::Bool = false
    random_algo::RA = MFSBruteForce(100, 100, 0.99)
    bbkwargs::BB = NamedTuple()
end

function _mfs(algorithm::MFSBlackBoxOptim, mapper, u0, search_area, idchecker, metric)
    function objective_function(perturbation)
        return mfs_objective(perturbation, u0, idchecker, metric, mapper, algorithm.penalty)
    end
    dim = dimension(mapper.ds)
    if algorithm.print_info == true
        TraceMode = :compact
    else
        TraceMode = :silent
    end

    if isnothing(algorithm.random_algo) && isnothing(algorithm.guess)
            result = bboptimize(objective_function;
            MaxSteps = algorithm.max_steps, SearchRange = search_area,
            NumDimensions = dim, TraceMode
        )
    else
        if !isnothing(algorithm.guess)
            guess = algorithm.guess
        else
            guess = minimal_fatal_shock(mapper, u0, search_area, algorithm.random_algo)
        end
        result = bboptimize(objective_function, guess;
            MaxSteps = algorithm.max_steps,
            SearchRange = search_area,
            NumDimensions = dim, TraceMode,
            algorithm.bbkwargs...
        )
    end
    best_shock = best_candidate(result)
    return best_shock
end

function mfs_objective(perturbation, u0, idchecker, metric, mapper::AttractorMapper, penalty)
    dist = metric(perturbation)
    if dist == 0
        return penalty
    end
    new_shock = perturbation + u0
    new_id = mapper(new_shock)
    if idchecker(new_id) # valid shock, it brings us to a desired basin
        return dist
    else
        return dist + penalty
    end
end
