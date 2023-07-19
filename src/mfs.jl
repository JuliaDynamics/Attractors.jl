using BlackBoxOptim: bboptimize, best_candidate
using Random: GLOBAL_RNG
using StateSpaceSets: statespace_sampler
export minimal_fatal_shock, MFSBruteForce, MFSBlackBoxOptim

"""
    MFSBruteForce(; kwargs...)

Initialize the randomised algorithm used in `minimal_fatal_shock` function. It consists of
two steps: random initialization and sphere radius reduction. On the first step,
the algorithm generates random pertubations of the initial point `u0` and with the best 
obtained pertubation it proceeds to the second step. On the second step, the algorithm 
generates random pertubations  on the surface of the hypersphere with radius equal to the 
norm of the best pertubation found in the first step. It reduces the radius 
of the hypersphere to continue searching for the better result within smaller radius. 

## Keyword arguments
`initial_iterations` number of random pertubations to try in the first step of the 
algorithm, default = 10000.
`sphere_iterations` number of steps while initializing random points on hypersphere and
decreasing its radius, default = 10000.
`sphere_decrease_step` factor by which the radius of the hypersphere is decreased, 
default = 1000.0, the higher the value, the more precise values may be obtained.

## Description
 
"""
mutable struct MFSBruteForce 
    initial_iterations::Int64
    sphere_iterations::Int64
    sphere_decrease_step::Float64
end

function MFSBruteForce(; initial_iterations = 10000, sphere_iterations = 10000, 
                                                    sphere_decrease_step = 100.0)

    MFSBruteForce(initial_iterations, sphere_iterations, sphere_decrease_step)
end




"""
    MFSBlackBoxOptim(; kwargs...)

Initialize the optimization algorithm used in `minimal_fatal_shock` function. It uses 
BlackBoxOptim.jl package to find the best shock. It is based on derivative free 
optimization and uses the objective function with penalties to find the minimal fatal shock.

## Keyword arguments
`guess` vector of initial guesses for the optimization algorithm, `default = []`.
`MaxSteps` maximum number of steps for the optimization algorithm, default = 10000.
`penalty` penalty value for the objective function, allows to adjust optimization algorithm 
to find the minimal fatal shock, `default = 1000.0`
`random_algo` algorithm used to find the initial guess for the optimization algorithm, 
by default it is initialized as `MFSBruteForce(0,0,0)` and not used. To activate it,
you need to initialize it with the parameters you want to use, 
e.g. `MFSBruteForce(1000,1000,100.0)` or `MFSBruteForce()` with default parameters.

"""
struct MFSBlackBoxOptim
    guess::Vector{Float64}
    MaxSteps::Int64
    penalty::Float64
    random_algo::MFSBruteForce
end

function MFSBlackBoxOptim(; guess = [], MaxSteps = 10000,  penalty = 1000.0, 
                                            random_algo = MFSBruteForce(0,0,0))

    MFSBlackBoxOptim(guess, MaxSteps, penalty, random_algo)
end


"""
    minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, algorithm)  

Run the minimal fatal shock algorithm on the initial point `u0` and output tuple 
`minimal fatal shock` and its `norm`. Two algorithms are available: [`MFSBruteForce`](@ref) 
and [`MFSBlackBoxOptim`](@ref). 

`mapper::AttractorMapper` one of available in Attractors.jl [`AttractorMapper`](@ref) 
constructed with respect to dynamical system.
`u0` vector containing coordinates of initial point to be tested.
`search_area` array of tuples or a tuple, by specifying single value 
e.g.  `(-1.5, 1.5)` or `[(-1.5, 1.5)]` you specify this value for each of the variables, 
or you may specify separetely 
search area for each of your states by, for example in 2 dimensions, `[(-2, 2), (-1, 1)]`.
`algorithm` one of the two algorithms available: [`MFSBruteForce`](@ref) 
or [`MFSBlackBoxOptim`](@ref).

## Setup 
While working with high dimensional systems, we recommend using [`MFSBlackBoxOptim`](@ref) 
algorithm, as it will give more precise results. While you can use [`MFSBruteForce`](@ref) 
algorithm to quickly obtain some guess which you can pass as the argument in initialization
of [`MFSBlackBoxOptim`](@ref) algorithm to optimize it and obtain exact results more 
efficiently. Or you can also pass [`MFSBruteForce`](@ref) algorithm as the argument in 
initialization of [`MFSBlackBoxOptim`](@ref) with parameter random_algo. We recommend to pay
attention to setup parameters of the algorithms, as default settings may not be sufficient 
to obtain precise results in some cases. Parameters `MaxSteps` and `search_area` are crucial 
in initialization of [`MFSBlackBoxOptim`](@ref) algorithm. The higher the `MaxSteps` value, 
the more precise results may be obtained with the cost of longer computation time. 
By manually decresing approximaton of `search_area`, you may significantly optimize 
algorithm's performance.

## Description
The minimal fatal shock algorithm is used to find the smallest perturbation of the initial 
point `u0` that will lead to a different basin of attraction. It is inspired by the paper 
"Minimal fatal shocks in multistable complex networks"[^MFS2020] by L. Halekotte and U. Feudel.

[^MFS2020]:
    Halekotte, L., Feudel, U. Minimal fatal shocks in multistable complex networks. 
    [Sci Rep 10, 11783 (2020).](https://doi.org/10.1038/s41598-020-68805-6)

"""
function minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, 
                                                algorithm::MFSBruteForce) 

    id_u0 = mapper(u0)
    dim = dimension(mapper.ds)
    if typeof(search_area) <: Tuple{Any,Any}
        search_area  = [search_area for _ in 1:dim]
    elseif length(search_area) == 1
        search_area = [search_area[1] for _ in 1:dim]
    elseif length(search_area) != dim
        error("Input search area does not match the dimension of the system")
    end

    best_shock, best_dist = crude_initial_radius(mapper, u0, search_area, dim, id_u0; 
      total_iterations = algorithm.initial_iterations)
    best_shock, best_dist = mfs_brute_force(mapper, u0, best_shock, best_dist, dim, id_u0, 
     algorithm.sphere_iterations, algorithm.sphere_decrease_step)
    return best_shock, best_dist

end

function minimal_fatal_shock(mapper::AttractorMapper, u0, search_area, 
                                                algorithm::MFSBlackBoxOptim)
    id_u0 = mapper(u0)
    dim = dimension(mapper.ds)
    if typeof(search_area) <: Tuple{Any,Any}
        search_area  = [search_area for _ in 1:dim]
    elseif length(search_area) == 1
        search_area = [search_area[1] for _ in 1:dim]
    elseif length(search_area) != dim
        error("Input search area does not match the dimension of the system")
    end

    function objective_function(perturbation)
        return mfs_objective(perturbation, u0, id_u0, mapper, algorithm.penalty)
    end
    
    rand_guess = []
    if algorithm.random_algo.initial_iterations != 0
        rand_guess = minimal_fatal_shock(mapper, u0, search_area, algorithm.random_algo)[1]
    end

    if (algorithm.guess != [] || rand_guess != [])  && dim == length(algorithm.guess)
        result = bboptimize(objective_function, algorithm.guess; 
                                            MaxSteps = algorithm.MaxSteps, 
                                            SearchRange = search_area, 
                                            NumDimensions = dim)
    else
        result = bboptimize(objective_function; MaxSteps = algorithm.MaxSteps, 
                                                SearchRange = search_area, 
                                                NumDimensions = dim)
    end

    best_shock = best_candidate(result)
    best_dist = norm(best_shock)
    return best_shock, best_dist

end





""" 
    crude_initial_radius(mapper::AttractorMapper, u0, search_area, dim,  
                                                    total_iterations=10000)

This function generates a random pertubation of the initial point `u0` within 
specified "search_area" and checks if it is in the same basin of attraction.
It does so by generating a random vector of length dim and then adding it to u0.
If the pertubation is not in the same basin of attraction, it calculates the norm 
of the pertubation and compares it to the best pertubation found so far.
If the norm is smaller, it updates the best pertubation found so far.
It repeats this process total_iterations times and returns the best pertubation found.

`mapper`  one of available in Attractors.jl [`AttractorMapper`](@ref) constructed with 
respect to dynamical system
`u0` initial point to be tested
`search_area` array of two points defining the search area
`dim` dimension of the system
`total_iterations` number of random pertubations to try, `default = 10000`

"""
function crude_initial_radius(mapper::AttractorMapper, u0, search_area, dim, id_u0;  
                                                                 total_iterations=10000)
    best_dist = Inf
    best_shock = nothing
    for _ in 1:total_iterations
        perturbation = rand(Uniform(search_area[1][1],search_area[1][2]), dim)
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
    mfs_brute_force(mapper, u0, best_shock, best_dist, dim,  total_iterations=10000)

This function works on the results obtained by `crude_initial_radius`. It starts from 
the best shock found so far and tries to find a better one by continuously reducing 
the radius of the sphere on the surface of which it generates random pertubations. 
If pertubation with the same basin of attraction is found, it updates the best shock found 
so far and reduces the radius of the sphere. It repeats this process total_iterations times
and returns the best pertubation found.

`mapper` one of available in Attractors.jl [`AttractorMapper`](@ref) constructed 
with respect to dynamical system.
`u0` initial point to be tested.
`best_shock` best shock found by `crude_initial_radius`.
`best_dist` norm of the best shock found by `crude_initial_radius`.
`dim` dimension of the system.
`total_iterations` number of random pertubations to try, default = 10000.
`sphere_decrease_step` factor by which the radius of the hypersphere is decreased, 
default = 1000.0, the higher the value the more precise values may be obtained.
"""
function mfs_brute_force(mapper::AttractorMapper, u0, 
                        best_shock, best_dist, dim, id_u0,  
                        total_iterations=10000, sphere_decrease_step = 1000.0)

    temp_dist = best_dist-best_dist/sphere_decrease_step
    perturbation = zeros(dim)

    for _ in 1:total_iterations
        generator, _ = statespace_sampler(GLOBAL_RNG;
                                         radius = 1.0, spheredims = dim, 
                                         center = zeros(dim) ) 
        perturbation = generator() * temp_dist
        new_shock = perturbation + u0
        
        if !(id_u0 == mapper(new_shock))
            best_dist = norm(perturbation)
            best_shock = perturbation
            temp_dist = best_dist-best_dist/sphere_decrease_step
        end
    end

    return best_shock, best_dist
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

