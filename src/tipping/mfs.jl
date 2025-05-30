using BlackBoxOptim: bboptimize, best_candidate
using Random: GLOBAL_RNG
import LinearAlgebra
export minimal_critical_shock, MCSBruteForce, MCSBlackBoxOptim
export excitability_threshold

"""
    minimal_critical_shock(mapper::AttractorMapper, u0, search_area, algorithm; kw...)

Return the _minimal critical shock_ for the initial point `u0` according to the
specified `algorithm` given a `mapper` that satisfies the `id = mapper(u0)` interface
(see [`AttractorMapper`](@ref) if you are not sure which mappers do that).
The output `mfs` is a vector like `u0`.

The `mapper` contains a reference to a [`DynamicalSystem`](@ref).
The options for `algorithm` are: [`MCSBruteForce`](@ref) or [`MCSBlackBoxOptim`](@ref).
For high dimensional systems [`MCSBlackBoxOptim`](@ref) is likely more accurate.

The `search_area` dictates the state space range for the search of the `mfs`.
It can be a 2-tuple of (min, max) values, in which case the same values are used
for each dimension of the system in `mapper`. Otherwise, it can be a vector of 2-tuples,
each for each dimension of the system. The search area is defined w.r.t. to `u0`
(i.e., it is the search area for perturbations of `u0`).

An alias to `minimal_critical_shock` is `excitability_threshold`.
Other names for the concept are or _stability threshold_ or _minimal fatal shock_.

## Keyword arguments

- `metric = LinearAlgebra.norm`: a metric function that gives the norm of a perturbation vector.
  This keyword is ignored for the [`MCSBruteForce`](@ref) algorithm.
- `target_id = nothing`: when not `nothing`, it should be an integer or a vector of
  integers corresponding to target attractor label(s).
  Then, the MFS is estimated based only on perturbations that lead to the target
  attractor(s).

## Description

The minimal critical shock is defined as the smallest-norm perturbation of the initial
point `u0` that will lead it a different basin of attraction than the one it was originally in.
This alternative basin is not returned, do `mapper(u0 .+ mfs)` if you need the ID.

The minimal critical shock has many names. Many papers computed this quantity without explicitly
naming it, or naming it something simple like "distance to the threshold".
The first work that proposed the concept as a nonlocal stability quantifier
was by [Klinshov2015](@cite) with the name "stability threshold".
Here we use the name of [Halekotte2020](@cite).

Our implementation is generic and works for _any_ dynamical system,
using either black box optimization or brute force searching approaches
and the unique interface of Attractors.jl for mapping initial conditions to attractors.
In contrast to [Klinshov2015](@cite) or [Halekotte2020](@cite), our implementation
does not place any assumptions on the nature of the dynamical system, or whether
the basin boundaries are smooth.

The _excitability threshold_ is a concept nearly identical, however, instead of looking
for a perturbation that simply brings us out of the basin, we look for the smallest
perturbation that brings us into specified basin(s). This is enabled via the keyword
`target_id`.
"""
function minimal_critical_shock(mapper::AttractorMapper, u0, search_area, algorithm;
        metric = LinearAlgebra.norm, target_id = nothing
    )
    dim = length(u0)
    if typeof(search_area) <: Tuple{<:Real,<:Real}
        search_area = [search_area for _ in 1:dim]
    elseif length(search_area) != dim
        error("Input search area does not match the dimension of the system")
    end
    # generate a function that returns `true` for ids that that are in the target basin
    id_u0 = mapper(u0)
    idchecker = id_check_function(id_u0, target_id)
    return _mfs(algorithm, mapper, u0, search_area, idchecker, metric)
end
const excitability_threshold = minimal_critical_shock

id_check_function(id::Int, ::Nothing) = i -> i != id
function id_check_function(id::Integer, ti::Integer)
    id == ti && error("target id and attractor id of u0 are the same.")
    return i -> i == ti
end
function id_check_function(id::Integer, ti::AbstractVector{<:Integer})
    id ∈ ti && error("attractor id of u0 belongs in target ids.")
    return i -> i ∈ ti
end


"""
    MCSBruteForce(; kwargs...)

The brute force randomized search algorithm used in [`minimal_critical_shock`](@ref).

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

Because this algorithm is based on hyperspheres, it assumes the Euclidean norm as the metric.

## Keyword arguments

- `initial_iterations = 10000`: number of random perturbations to try in the first step of the
  algorithm.
- `sphere_iterations = 10000`: number of steps while initializing random points on hypersphere and
  decreasing its radius.
- `sphere_decrease_factor = 0.999`: factor by which the radius of the hypersphere is decreased
  (at each step the radius is multiplied by this number). Number closer to 1 means
  more refined accuracy.
- `seed = rand(1:10000)`: seed for the random number generator used when sampling
  random perturbations.
"""
Base.@kwdef struct MCSBruteForce
    initial_iterations::Int = 10000
    sphere_iterations::Int = 10000
    sphere_decrease_factor::Float64 = 0.999
    seed::Int = rand(1:10000)
end

function _mfs(algorithm::MCSBruteForce, mapper, u0, search_area, idchecker, _metric)
    metric = LinearAlgebra.norm
    algorithm.sphere_decrease_factor ≥ 1 && error("Sphere decrease factor cannot be ≥ 1.")
    dim = length(u0)
    best_shock, best_dist = crude_initial_radius(
        mapper, u0, search_area, idchecker, metric, algorithm.initial_iterations, algorithm.seed
    )
    best_shock, best_dist = mfs_brute_force(
        mapper, u0, best_shock, best_dist, dim, idchecker, metric,
        algorithm.sphere_iterations, algorithm.sphere_decrease_factor, algorithm.seed
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
function crude_initial_radius(mapper::AttractorMapper, u0, search_area, idchecker, metric, total_iterations, seed)
    best_dist = Inf
    region = StateSpaceSets.HRectangle([s[1] for s in search_area], [s[2] for s in search_area])
    generator, _ = statespace_sampler(region, seed)
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
        total_iterations, sphere_decrease_factor, seed
    )

    temp_dist = best_dist*sphere_decrease_factor
    region = HSphereSurface(temp_dist, dim)
    generator, = statespace_sampler(region, seed)
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
    MCSBlackBoxOptim(; kwargs...)

The black box derivative-free optimization algorithm used in [`minimal_critical_shock`](@ref).

## Keyword arguments

- `guess = nothing`: a initial guess for the minimal critical shock given to the
  optimization algorithm. If not `nothing`, `random_algo` below is ignored.
- `max_steps = 10000`: maximum number of steps for the optimization algorithm.
- `penalty = 1000.0`: penalty value for the objective function for perturbations that do
  not lead to a different basin of attraction. This value is added to the norm of the
  perturbation and its value should be much larger than the typical sizes of the basins of
  attraction.
- `print_info`: boolean value, if true, the optimization algorithm will print information on
  the evaluation steps of objective function, `default = false`.
- `random_algo = MCSBruteForce(100, 100, 0.99)`: an instance of [`MCSBruteForce`](@ref)
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
which is why the output of a crude [`MCSBruteForce`](@ref) is used to provide a guess.
This can be disabled by either passing a `guess` vector explicitly or
by giving `nothing` as `random_algo`.
"""
Base.@kwdef struct MCSBlackBoxOptim{G, RA, BB}
    guess::G = nothing
    max_steps::Int64 = 10_000
    penalty::Float64 = 0.999
    print_info::Bool = false
    random_algo::RA = MCSBruteForce(100, 100, 0.99, rand(1:10000))
    bbkwargs::BB = NamedTuple()
end

function _mfs(algorithm::MCSBlackBoxOptim, mapper, u0, search_area, idchecker, metric)
    function objective_function(perturbation)
        return mfs_objective(perturbation, u0, idchecker, metric, mapper, algorithm.penalty)
    end
    dim = length(u0)
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
            guess = minimal_critical_shock(mapper, u0, search_area, algorithm.random_algo)
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
