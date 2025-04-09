export stagger_and_step
using LinearAlgebra: norm
using Random
using ProgressMeter

"""
    stagger_and_step(ds::DynamicalSystem, x0, N::Int, 
        isinside::Function; kwargs...) -> trajectory

Implement the stagger-and-step method 
[Sweet2001](@cite) and [Sala2016](@cite) to approximate the invariant 
non-attracting set governing the chaotic transient dynamics 
of a system, namely the stable manifold of a chaotic saddle. 

Given the dynamical system `ds` and a initial guess `x0` in a 
region *with no attractors* defined by the membership function 
`isinside`, the algorithm returns a `StateSpaceState` named `trajectory` 
of `N` points from the phase space close to the stable manifold. 
If we set one of this point as an initial condition of `ds`, 
the trajectory escape from the region after at least `Tm` 
steps of `ds`. The search is stochastic and depends on the 
parameter `δ` defining a (small) neighborhood of search. 

The function `isinside(x)` must return `true` if the point `x` is 
inside the chosen bounded region and `false` otherwise. See 
[`statespace_sampler`](@ref) as a helper to construct this 
function.

## Keyword arguments
* `δ = 1e-10`: A small number constraining the random
  search around a particular point. The interpretation of this 
  number will depend on the probability distribution chosen 
  for the sampling (see `stagger_mode`).

* `Tm = 30`: The minimum number of iterations of `ds` before 
  the trajectory escapes from the bounding box defined by
  `isinside`. 

* `max_steps = 10^5`: The search for a new candidate point may 
  fail at some point. If the search fails after `max_steps`, 
  a new initial point is set and the method starts from a new
  point.

* `max_escape_time = 10000`: If the trajectory stays in the 
  defined region after `max_escape_time` steps, there is probably
  an attractor in the region and the algorithm will throw an error. 

* `stagger_mode = :exp`: There are several ways to produce 
  candidate points `x` to fulfill the condition `T(x) > Tm`. 
  The available methods are: 

    * `:exp`: A candidate is sampled from a truncated exponential 
      distribution in a random direction `u` around the current
      `x` such that `x_c = x + u*r`. `r = 10^-s` with `s` taken
      from a uniform distribution in [-15, δ]. 

    * `:unif`: The next candidate is `x_c = x + u*r` with `r` 
      taken from a uniform distribution [0,δ] and `u` a random 
      direction around `x`. 

    * `:adaptive`: The next candidate is `x_c = x + u*r` with 
      `r` drawn from a gaussian distribution with variance  δ 
      and mean zero.
      The variance is adapted according to a free parameter `γ`
      such that: `δ = δ/γ` if no candidate is found and `δ = δ*γ`
      when it succeeds. 

* `γ = 1.1`: It is free parameter for the adaptive stagger method `:adaptive`. 

* `δ₀ = 1.0`: This is the radius for the first stagger 
  trajectory search. The algorithm looks for a point 
  sufficiently close to the saddle before switching to the 
  stagger-and-step routine. The search radius must be large 
  enough to find a suitable initial candidate. 
  To type `δ₀` use `\delta<TAB>\_0<TAB>`.

* `rng::AbstractRNG = Xoshiro()`:  Random number generator. Use this for
  reproducibility.

## Description 

The method relies on the stagger-and-step algorithm that 
search initial conditions close to the saddle with escapes time 
`T(x_n) > Tm`. The function `T` represents the iteration number 
at which the trajectory with initial condition `x_n` steps out 
from a region defined by the user (see the argument `isinside`).

Given the dynamical mapping `F`, if the iteration `x_{n+1} = 
F(x_n)` respects the condition `T(x_{n+1}) > Tm` we accept 
this next point, this is the _step_ part of the method. If 
not, the method search randomly the next point in a 
neighborhood following a given probability distribution, this 
is the _stagger_ part. This part sometimes fails to find a new 
candidate and a new starting point of the trajectory is chosen 
within the defined region. See the keyword argument 
`stagger_mode` for the different available methods.   

The method produces a pseudo-trajectory of `N` points δ-close
to the stable manifold of the chaotic saddle. 

"""
function stagger_and_step(ds::DynamicalSystem, x0, N::Int, isinside::Function; δ = 1e-10, Tm = 30, 
    γ = 1.1, max_steps = Int64(1e5), max_escape_time = 10000, stagger_mode = :exp, δ₀ = 1.0, 
    show_progress = true,  rng::AbstractRNG = Xoshiro())

    progress = ProgressMeter.Progress(
        N; desc = "Saddle estimation: ", dt = 1.0
    )
    xi = stagger_trajectory(ds, x0, Tm, isinside; δ₀, stagger_mode = :unif, 
                            max_steps, γ, max_escape_time, rng) 
    if isnothing(xi)
        error("Cannot find a stagger trajectory. Choose a different starting point or 
                search radius δ₀.")
    end

    v = Vector{typeof(current_state(ds))}(undef, N)
    v[1] = xi
    for n in 1:N
        show_progress && ProgressMeter.update!(progress, n)
        if escape_time!(ds, xi, isinside; max_escape_time) > Tm
            reinit!(ds, xi)
        else
            xp, Tp = stagger!(ds, xi, δ, Tm, isinside; stagger_mode, max_steps, γ, max_escape_time, rng)
            # The stagger step may fail. We reinitiate the algorithm from a new initial condition.
            if Tp < 0
                xp = stagger_trajectory(ds, x0, Tm, isinside; δ₀, stagger_mode = :exp,max_steps, γ, max_escape_time, rng) 
                if isnothing(xp)
                    error("Cannot find a stagger trajectory. Choose a different starting 
                          point or search radius δ₀.")
                end
                δ = 0.1
            end
            reinit!(ds, xp)
        end 
        step!(ds)
        xi = copy(current_state(ds))
        v[n] = xi
    end
    return StateSpaceSet(v)
end


"""
    stagger_trajectory(ds, x0, Tm, isinside; kwargs...) -> xi

On success, the function returns a point `xi` with the property `T(xi) > 
Tm` with a random walk search around the initial coordinates 
`x0`. `T(xi)` is the escape time of the initial condition `xi` from 
the bounding box defined by the function `isinside`.  
In the case where the algorithm cannot find a suitable point, the algorithm 
returns nothing. 

This is an auxiliary function for [`stagger_and_step`](@ref). 
Keyword arguments and definitions are identical for both functions. 
"""
function stagger_trajectory(ds, x0, Tm, isinside; δ₀ = 1., stagger_mode = :exp,
        max_steps = Int(1e5), γ = 1.1, max_escape_time = 10000, rng::AbstractRNG = Xoshiro())
    T = escape_time!(ds, x0, isinside; max_escape_time)
    xi = copy(x0) 
    while !(T > Tm)  # we must have T > Tm at each step 
        xi, T = stagger!(ds, xi, δ₀, T, isinside; γ, stagger_mode, max_steps, max_escape_time, rng)
        if T < 0
            return nothing
        end 
            
    end
    return xi
end
    


function escape_time!(ds, x0, isinside; max_escape_time = 10000) 
    set_state!(ds,x0)
    reinit!(ds, x0)
    k = 1; 
    while isinside(current_state(ds)) 
        if k > max_escape_time 
            error("The trajectory did not escape for ", k, " steps, you probably have \
an attractor in the defined region. Last point evaluated: ", current_state(ds))
        end
        step!(ds)
        k += 1
    end
    return current_time(ds)
end

function rand_u!(u, δ, n, stagger_mode, rng)
    if stagger_mode == :exp 
        a = -log10(δ)
        s = (15-a)*rand(rng) + a
        u .= randn(rng,n)
        u .*= (10.0^(-s))/norm(u)
        return
    elseif stagger_mode == :unif
        s = δ*rand(rng)
        u .= randn(rng,n)
        u .*= s/norm(u)
        return
    elseif stagger_mode == :adaptive
        s = δ*randn(rng)
        u .= randn(rng,n)
        u .*= s/norm(u)
        return
    else
        error("Invalid stagger_mode: $stagger_mode")
    end
end

"""
    stagger!(ds, x0, δ, Tm, isinside; kwargs...) -> stagger_point

This function searches a new candidate in a neighborhood of x0 with a random search 
depending on some distribution. If the search fails it returns a negative time.
"""
function stagger!(ds, x0, δ, Tm, isinside; max_steps = Int(1e6), γ = 1.1, stagger_mode = :exp, verbose = false, max_escape_time = 10000, rng::AbstractRNG)
    Tp = 0; xp = zeros(length(x0)); k = 1; 
    u = zeros(length(x0))
    T0 = escape_time!(ds, x0, isinside; max_escape_time)
    if !isinside(x0)
        error("x0 must be in grid")
    end
    while Tp ≤ Tm 
        rand_u!(u, δ, length(x0), stagger_mode, rng)
        xp .= x0 .+ u 

        if k > max_steps 
           if verbose 
               @warn "Stagger search fails, δ is too small or T is too large. 
                We reinitiate the algorithm
               "
           end
           return xp, -1
        end
        Tp = escape_time!(ds, xp, isinside; max_escape_time)
        if stagger_mode == :adaptive
            # We adapt the variance of the search
            # if the alg. can't find a candidate
            if Tp < T0
                δ = δ/γ
            elseif Tp == T0
                δ = δ*γ
            end
            if Tp == Tm
                # The adaptive alg. accepts T == Tp
                return xp, Tp
            end
        end
        k = k + 1
    end
    return xp, Tp
end

