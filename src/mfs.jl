using BlackBoxOptim: bboptimize, best_candidate
using Random: GLOBAL_RNG
using StateSpaceSets: statespace_sampler
export minimal_fatal_shock, MFSBruteForce, MFSBlackBoxOptim


"""
    minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, algorithm) → mfs

Return the minimal fatal shock `mfs` for the initial point `u0` according to the
specified `algorithm` given a `mapper` that satisfies the `id = mapper(u0)` interface
(see [`AttractorMapper`](@ref) if you are not sure which mappers do that).
The `mapper` contains a reference to a [`DynamicalSystem`](@ref).
The options for `algorithm` are: [`MFSBruteForce`](@ref) or [`MFSBlackBoxOptim`](@ref).

The `search_area` dictates the state space range for the search of the `mfs`.
It can be a 2-tuple of (min, max) values, in which case the same values are used
for each dimension of the system in `mapper`. Otherwise, it can be a vector of 2-tuples,
each for each dimension of the system.

## Description

The minimal fatal shock is defined as the smallest perturbation of the initial
point `u0` that will lead it a different basin of attraction. It is inspired by the paper
"Minimal fatal shocks in multistable complex networks" [Halekotte2020](@cite),
however the implementation here is generic: it works for _any_ dynamical system.

While working with high dimensional systems, we recommend using [`MFSBlackBoxOptim`](@ref)
algorithm, as it will give more precise results. While you can use [`MFSBruteForce`](@ref)
algorithm to quickly obtain some guess which you can pass as the argument in initialization
of [`MFSBlackBoxOptim`](@ref) algorithm to optimize it and obtain exact results more
efficiently. Or you can also pass [`MFSBruteForce`](@ref) algorithm as the argument in
initialization of [`MFSBlackBoxOptim`](@ref) with parameter `random_algo`. We recommend to pay
attention to setup parameters of the algorithms, as default settings may not be sufficient
to obtain precise results in some cases. Parameters `MaxSteps` and `search_area` are crucial
in initialization of [`MFSBlackBoxOptim`](@ref) algorithm. The higher the `MaxSteps` value,
the more precise results may be obtained with the cost of longer computation time.
By manually decresing approximaton of `search_area`, you may significantly optimize
algorithm's performance.
"""
function minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, algorithm)
    id_u0 = mapper(u0)
    dim = dimension(mapper.ds)
    if typeof(search_area) <: Tuple{Any,Any}
        search_area  = [search_area for _ in 1:dim]
    elseif length(search_area) == 1
        search_area = [search_area[1] for _ in 1:dim]
    elseif length(search_area) != dim
        error("Input search area does not match the dimension of the system")
    end
    return _mfs(algorithm, mapper, u0, search_area, id_u0)
end


"""
    MFSBruteForce(; kwargs...)

The brute force randomized search algorithm used in [`minimal_fatal_shock`](@ref).

It consists of
two steps: random initialization and sphere radius reduction. On the first step,
the algorithm generates random pertubations of the initial point `u0` and records
the perturbation that leads to a different basin but with the smallest magnitude.
With this
obtained pertubation it proceeds to the second step. On the second step, the algorithm
generates random pertubations on the surface of the hypersphere with radius equal to the
norm of the pertubation found in the first step. It reduces the radius
of the hypersphere to continue searching for the better result within smaller radius.
It records the shock with smallest radius that leads to a different basin.
(a new radius is recorded only if a perturbation that leads to a different basin is found)

## Keyword arguments

- `initial_iterations = 10000`: number of random pertubations to try in the first step of the
  algorithm.
- `sphere_iterations = 10000`: number of steps while initializing random points on hypersphere and
  decreasing its radius.
- `sphere_decrease_factor = 0.999` factor by which the radius of the hypersphere is decreased
  (at each step the radius is multiplied by this number).
"""
Base.@kwdef struct MFSBruteForce
    initial_iterations::Int64 = 10000
    sphere_iterations::Int64 = 10000
    sphere_decrease_factor::Float64 = 0.999
end

function _mfs(algorithm::MFSBruteForce, mapper, u0, search_area, id_u0)
    dim = dimension(mapper.ds)
    best_shock, best_dist = crude_initial_radius(
        mapper, u0, search_area, id_u0, algorithm.initial_iterations
    )
    best_shock, best_dist = mfs_brute_force(
        mapper, u0, best_shock, best_dist, dim, id_u0,
        algorithm.sphere_iterations, algorithm.sphere_decrease_factor
    )
    return best_shock
end

"""
This function generates a random pertubation of the initial point `u0` within
specified "search_area" and checks if it is in the same basin of attraction.
It does so by generating a random vector of length dim and then adding it to u0.
If the pertubation is not in the same basin of attraction, it calculates the norm
of the pertubation and compares it to the best pertubation found so far.
If the norm is smaller, it updates the best pertubation found so far.
It repeats this process total_iterations times and returns the best pertubation found.
"""
function crude_initial_radius(mapper::AttractorMapper, u0, search_area, id_u0, total_iterations)
    best_dist = Inf
    best_shock = nothing
    generator, _ = statespace_sampler(; min_bounds = [s[1] for s in search_area],
        max_bounds = [s[2] for s in search_area]
    )

    for _ in 1:total_iterations
        perturbation = generator()
        shock = u0 + perturbation
        if !(id_u0 == mapper(shock))
            dist = norm(perturbation)

            if dist < best_dist
                best_dist = dist
                best_shock = perturbation

            end
        end
    end

    return best_shock, best_dist
end


"""
This function works on the results obtained by `crude_initial_radius`. It starts from
the best shock found so far and tries to find a better one by continuously reducing
the radius of the sphere on the surface of which it generates random pertubations.
If pertubation with the same basin of attraction is found, it updates the best shock found
so far and reduces the radius of the sphere. It repeats this process total_iterations times
and returns the best pertubation found.
"""
function mfs_brute_force(mapper::AttractorMapper, u0,
                        best_shock, best_dist, dim, id_u0,
                        total_iterations, sphere_decrease_factor)

    temp_dist = best_dist*sphere_decrease_factor
    perturbation = zeros(dim)

    for _ in 1:total_iterations
        generator, _ = statespace_sampler(; radius = temp_dist, spheredims = dim)
        perturbation = generator()
        new_shock = perturbation + u0

        if !(id_u0 == mapper(new_shock))
            best_dist = norm(perturbation)
            best_shock = perturbation
            temp_dist = best_dist*sphere_decrease_factor
        end
    end

    return best_shock, best_dist
end


"""
    MFSBlackBoxOptim(; kwargs...)

The black box derivative-free optimization algorithm used in [`minimal_fatal_shock`](@ref).


## Keyword arguments

- `guess = nothing` a initial guess for the minimal fatal shock given to the
  optimization algorithm. If not `nothing`, `random_algo` below is ignored.
- `MaxSteps = 10000` maximum number of steps for the optimization algorithm.
- `penalty` penalty value for the objective function, allows to adjust optimization algorithm
  to find the minimal fatal shock, `default = 1000.0`
- `PrintInfo` boolean value, if true, the optimization algorithm will print information on
  the evaluation steps of objective function, `default = false`.
- `random_algo = nothing` algorithm used to find the initial guess for the optimization algorithm,
  To activate it, you need to initialize and provide an instance of [`MFSBruteForce`](@ref).

## Description

Initialize the optimization algorithm used in `minimal_fatal_shock` function.
The algorithm uses BlackBoxOptim.jl and a penaltized objective function to minimize.
y function used as a constraint function.
So, if we hit another basin during the search we encourage the algorithm otherwise we
punish it with some penalty. The function to minimize is (besides some details):
```julia
function mfs_objective(perturbation, u0, mapper, penalty)
    dist = norm(perturbation)
    if mapper(u0 + perturbation) == mapper(u0)
        # penaltize if we stay in the same basin:
        return dist + penalty
    else
        return dist
    end
end
```
"""
Base.@kwdef struct MFSBlackBoxOptim{G, RA}
    guess::G = nothing
    MaxSteps::Int64 = 10_000
    penalty::Float64 = 0.999
    PrintInfo::Bool = false
    random_algo::RA = nothing
end

function _mfs(algorithm::MFSBlackBoxOptim, mapper, u0, search_are, id_u0)
    function objective_function(perturbation)
        return mfs_objective(perturbation, u0, id_u0, mapper, algorithm.penalty)
    end
    dim = dimension(mapper.ds)
    if algorithm.PrintInfo == true
        TraceMode = :compact
    else
        TraceMode = :silent
    end

    if isnothing(algorithm.random_algo) && isnothing(algorithm.guess)
            result = bboptimize(objective_function;
            MaxSteps = algorithm.MaxSteps, SearchRange = search_area,
            NumDimensions = dim, TraceMode
        )
    else
        if !isnothing(algorithm.guess)
            guess = algorithm.guess
        else
            guess = minimal_fatal_shock(mapper, u0, search_area, algorithm.random_algo)
        end
        result = bboptimize(objective_function, guess;
            MaxSteps = algorithm.MaxSteps,
            SearchRange = search_area,
            NumDimensions = dim, TraceMode
        )
    end
    best_shock = best_candidate(result)
    return best_shock
end

function mfs_objective(perturbation, u0, id_u0, mapper::AttractorMapper, penalty=1000.0)
    dist = norm(perturbation)
    if dist == 0
        return penalty
    end
    if mapper(u0 + perturbation) == id_u0
        return dist + penalty
    else
        return dist
    end
end
