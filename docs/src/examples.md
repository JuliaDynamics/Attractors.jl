# Examples for Attractors.jl

Note that the examples utilize some convenience plotting functions offered by Attractors.jl which come into scope when using `Makie` (or any of its backends such as `CairoMakie`), see the [visualization utilities](@ref) for more.

## Newton's fractal (basins of 2D map)

```@example MAIN
using Attractors
function newton_map(z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
xg = yg = range(-1.5, 1.5; length = 400)
# Use non-sparse for using `basins_of_attraction`
mapper_newton = AttractorsViaRecurrences(ds, (xg, yg);
    sparse = false, consecutive_lost_steps = 1000
)
basins, attractors = basins_of_attraction(mapper_newton; show_progress = false)
basins
```
```@example MAIN
attractors
```

Now let's plot this as a heatmap, and on top of the heatmap, let's scatter plot the attractors. We do this in one step by utilizing one of the pre-defined plotting functions offered by Attractors.jl

```@example MAIN
using CairoMakie
xg = yg = range(-1.5, 1.5; length = 400)
grid = (xg, yg)
fig = heatmap_basins_attractors(grid, basins, attractors)
```

Instead of computing the full basins, we could get only the fractions of the basins of attractions using [`basins_fractions`](@ref), which is typically the more useful thing to do in a high dimensional system.
In such cases it is also typically more useful to define a sampler that generates initial conditions on the fly instead of pre-defining some initial conditions (as is done in [`basins_of_attraction`](@ref). This is simple to do:

```@example MAIN
sampler, = statespace_sampler(grid)

basins = basins_fractions(mapper_newton, sampler)
```

in this case, to also get the attractors we simply extract them from the underlying storage of the mapper:
```@example MAIN
attractors = extract_attractors(mapper_newton)
```


## Minimal Fatal Shock
Here we find the Minimal Fatal Shock (MFS, see [`minimal_fatal_shock`](@ref)) for the attractors (i.e., fixed points) of Newton's fractal
```@example MAIN
shocks = Dict()
algo_bb = Attractors.MFSBlackBoxOptim()
for atr in values(attractors)
    u0 = atr[1]
    shocks[u0] = minimal_fatal_shock(mapper_newton, u0, (-1.5,1.5), algo_bb)
end
shocks
```
To visualize results we can make use of previously defined heatmap
```@example MAIN
ax =  content(fig[1,1])
for (atr, shock) in shocks
    lines!(ax, [atr, atr + shock]; color = :orange, linewidth = 3)
end
fig
```

## Fractality of 2D basins of the (4D) magnetic pendulum
In this section we will calculate the basins of attraction of the four-dimensional magnetic pendulum. We know that the attractors of this system are all individual fixed points on the (x, y) plane so we will only compute the basins there. We can also use this opportunity to highlight a different method, the [`AttractorsViaProximity`](@ref) which works when we already know where the attractors are. Furthermore we will also use a `ProjectedDynamicalSystem` to project the 4D system onto a 2D plane, saving a lot of computational time!

### Computing the basins

First we need to load in the magnetic pendulum from the predefined dynamical systems library
```@example MAIN
using Attractors, CairoMakie
using PredefinedDynamicalSystems
ds = PredefinedDynamicalSystems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
```

Then, we create a projected system on the x-y plane
```@example MAIN
psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])
```

For this systems we know the attractors are close to the magnet positions. The positions can be obtained from the equations of the system, provided that one has seen the source code (not displayed here), like so:
```@example MAIN
attractors = Dict(i => StateSpaceSet([dynamic_rule(ds).magnets[i]]) for i in 1:3)
```

and then create a
```@example MAIN
mapper = AttractorsViaProximity(psys, attractors)
```

and as before, get the basins of attraction
```@example MAIN
xg = yg = range(-4, 4; length = 201)
grid = (xg, yg)
basins, = basins_of_attraction(mapper, grid; show_progress = false)

heatmap_basins_attractors(grid, basins, attractors)
```

### Computing the uncertainty exponent
Let's now calculate the [`uncertainty_exponent`](@ref) for this system as well.
The calculation is straightforward:
```@example MAIN
using CairoMakie
ε, f_ε, α = uncertainty_exponent(basins)
fig, ax = lines(log.(ε), log.(f_ε))
ax.title = "α = $(round(α; digits=3))"
fig
```
The actual uncertainty exponent is the slope of the curve (α) and indeed we get an exponent near 0 as we know a-priory the basins have fractal boundaries for the magnetic pendulum.

### Computing the tipping probabilities
We will compute the tipping probabilities using the magnetic pendulum's example
as the "before" state. For the "after" state we will change the `γ` parameter of the
third magnet to be so small, its basin of attraction will virtually disappear.
As we don't know _when_ the basin of the third magnet will disappear, we switch the attractor finding algorithm back to [`AttractorsViaRecurrences`](@ref).

```@example MAIN
set_parameter!(psys, :γs, [1.0, 1.0, 0.1])
mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 1)
basins_after, attractors_after = basins_of_attraction(
    mapper, (xg, yg); show_progress = false
)
# matching attractors is important!
rmap = match_statespacesets!(attractors_after, attractors)
# Don't forget to update the labels of the basins as well!
replace!(basins_after, rmap...)

# now plot
heatmap_basins_attractors(grid, basins_after, attractors_after)
```

And let's compute the tipping "probabilities":

```@example MAIN
P = tipping_probabilities(basins, basins_after)
```
As you can see `P` has size 3×2, as after the change only 2 attractors have been identified
in the system (3 still exist but our state space discretization isn't fine enough to
find the 3rd because it has such a small basin).
Also, the first row of `P` is 50% probability to each other magnet, as it should be due to
the system's symmetry.


## 3D basins via recurrences

To showcase the true power of [`AttractorsViaRecurrences`](@ref) we need to use a system whose attractors span higher-dimensional space. An example is
```@example MAIN
using Attractors
using PredefinedDynamicalSystems
ds = PredefinedDynamicalSystems.thomas_cyclical(b = 0.1665)
```
which, for this parameter, contains 3 coexisting attractors which are entangled
periodic orbits that span across all three dimensions.

To compute the basins we define a three-dimensional grid and call on it
[`basins_of_attraction`](@ref).

```julia
# This computation takes about an hour
xg = yg = zg = range(-6.0, 6.0; length = 251)
mapper = AttractorsViaRecurrences(ds, (xg, yg, zg); sparse = false)
basins, attractors = basins_of_attraction(mapper)
attractors
```
```
Dict{Int16, StateSpaceSet{3, Float64}} with 5 entries:
  5 => 3-dimensional StateSpaceSet{Float64} with 1 points
  4 => 3-dimensional StateSpaceSet{Float64} with 379 points
  6 => 3-dimensional StateSpaceSet{Float64} with 1 points
  2 => 3-dimensional StateSpaceSet{Float64} with 538 points
  3 => 3-dimensional StateSpaceSet{Float64} with 537 points
  1 => 3-dimensional StateSpaceSet{Float64} with 1 points
```

_Note: the reason we have 6 attractors here is because the algorithm also finds 3 unstable fixed points and labels them as attractors. This happens because we have provided initial conditions on the grid `xg, yg, zg` that start exactly on the unstable fixed points, and hence stay there forever, and hence are perceived as attractors by the recurrence algorithm. As you will see in the video below, they don't have any basin fractions_

The basins of attraction are very complicated. We can try to visualize them by animating the 2D slices at each z value, to obtain:

```@raw html
<video width="75%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/cyclical_basins.mp4?raw=true" type="video/mp4">
</video>
```

Then, we visualize the attractors to obtain:

```@raw html
<video width="75%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/cyclical_attractors.mp4?raw=true" type="video/mp4">
</video>
```

In the animation above, the scattered points are the attractor values the function [`AttractorsViaRecurrences`](@ref) found by itself. Of course, for the periodic orbits these points are incomplete. Once the function's logic understood we are on an attractor, it stops computing. However, we also simulated lines, by evolving initial conditions colored appropriately with the basins output.

The animation was produced with the code:
```julia
using GLMakie
fig = Figure()
display(fig)
ax = fig[1,1] = Axis3(fig; title = "found attractors")
cmap = cgrad(:dense, 6; categorical = true)

for i in keys(attractors)
    tr = attractors[i]
    markersize = length(attractors[i]) > 10 ? 2000 : 6000
    marker = length(attractors[i]) > 10 ? :circle : :rect
    scatter!(ax, columns(tr)...; markersize, marker, transparency = true, color = cmap[i])
    j = findfirst(isequal(i), bsn)
    x = xg[j[1]]
    y = yg[j[2]]
    z = zg[j[3]]
    tr = trajectory(ds, 100, SVector(x,y,z); Ttr = 100)
    lines!(ax, columns(tr)...; linewidth = 1.0, color = cmap[i])
end

a = range(0, 2π; length = 200) .+ π/4

record(fig, "cyclical_attractors.mp4", 1:length(a)) do i
    ax.azimuth = a[i]
end
```

## Basins of attraction of a Poincaré map

[`PoincareMap`](@ref) is just another discrete time dynamical system within the DynamicalSystems.jl ecosystem. With respect to Attractors.jl functionality, there is nothing special about Poincaré maps. You simply initialize one use it like any other type of system. Let's continue from the above example  of the Thomas cyclical system
```@example MAIN
using Attractors
using PredefinedDynamicalSystems
ds = PredefinedDynamicalSystems.thomas_cyclical(b = 0.1665);
```
The three limit cycles attractors we have above become fixed points in the Poincaré map (for appropriately chosen hyperplanes). Since we already know the 3D structure of the basins, we can see that an appropriately chosen hyperplane is just the plane `z = 0`. Hence, we define a Poincaré map on this plane:

```@example MAIN
plane = (3, 0.0)
pmap = PoincareMap(ds, plane)
```

We define the same grid as before, but now only we only use the x-y coordinates. This is because we can utilize the special `reinit!` method of the [`PoincareMap`](@ref), that allows us to initialize a new state directly on the hyperplane (and then the remaining variable of the dynamical system takes its value from the hyperplane itself).
```@example MAIN
xg = yg = range(-6.0, 6.0; length = 251)
grid = (xg, yg)
mapper = AttractorsViaRecurrences(pmap, grid; sparse = false)
```
All that is left to do is to call [`basins_of_attraction`](@ref):

```@example MAIN
basins, attractors = basins_of_attraction(mapper; show_progress = false);
```

```@example MAIN
heatmap_basins_attractors(grid, basins, attractors)
```
_just like in the example above, there is a fourth attractor with 0 basin fraction. This is an unstable fixed point, and exists exactly because we provided a grid with the unstable fixed point exactly on this grid_


## Irregular grid for `AttractorsViaRecurrences`
It is possible to provide an irregularly spaced grid to `AttractorsViaRecurrences`. This can make algorithm performance better for continuous time systems where the state space flow has significantly different speed in some state space regions versus others.

In the following example the dynamical system has only one attractor: a limit cycle. However, near the origin (0, 0) the timescale of the dynamics becomes very slow. As the trajectory is stuck there for quite a while, the recurrences algorithm may identify this region as an "attractor" (incorrectly). The solutions vary and can be to increase drastically the max time checks for finding attractors, or making the grid much more fine. Alternatively, one can provide a grid that is only more fine near the origin and not fine elsewhere.

The example below highlights that for rather coarse settings of grid and convergence thresholds, using a grid that is finer near (0, 0) gives correct results:

```@example MAIN
using Attractors, CairoMakie

function predator_prey_fastslow(u, p, t)
    α, γ, ϵ, ν, h, K, m = p
    N, P = u
    du1 = α*N*(1 - N/K) - γ*N*P / (N+h)
    du2 = ϵ*(ν*γ*N*P/(N+h) - m*P)
    return SVector(du1, du2)
end
γ = 2.5
h = 1
ν = 0.5
m = 0.4
ϵ = 1.0
α = 0.8
K = 15
u0 = rand(2)
p0 = [α, γ, ϵ, ν, h, K, m]
ds = CoupledODEs(predator_prey_fastslow, u0, p0)

fig = Figure()
ax = Axis(fig[1,1])

# when pow > 1, the grid is finer close to zero
for pow in (1, 2)
    xg = yg = range(0, 18.0^(1/pow); length = 200).^pow
    mapper = AttractorsViaRecurrences(ds, (xg, yg);
        Dt = 0.1, sparse = true,
        consecutive_recurrences = 10, attractor_locate_steps = 10,
        maximum_iterations = 1000,
    )

    # Find attractor and its fraction (fraction is always 1 here)
    sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
    fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
    attractors = extract_attractors(mapper)
    scatter!(ax, vec(attractors[1]); markersize = 16/pow, label = "pow = $(pow)")
end

axislegend(ax)

fig
```

## Subdivision Based Grid for `AttractorsViaRecurrences`

To achieve even better results for this kind of problematic systems than with previuosly introduced `Irregular Grids`  we provide a functionality to construct `Subdivision Based Grids` in which
one can obtain more coarse or dense structure not only along some axis but for a specific regions where the state space flow has
significantly different speed. [`subdivided_based_grid`](@ref) enables automatic evaluation of velocity vectors for regions of originally user specified
grid to further treat those areas as having more dense or coarse structure than others.

```@example MAIN
using Attractors, CairoMakie

function predator_prey_fastslow(u, p, t)
    α, γ, ϵ, ν, h, K, m = p
    N, P = u
    du1 = α*N*(1 - N/K) - γ*N*P / (N+h)
    du2 = ϵ*(ν*γ*N*P/(N+h) - m*P)
return SVector(du1, du2)
end
γ = 2.5
h = 1
ν = 0.5
m = 0.4
ϵ = 1.0
α = 0.8
K = 15
u0 = rand(2)
p0 = [α, γ, ϵ, ν, h, K, m]
ds = CoupledODEs(predator_prey_fastslow, u0, p0)

xg = yg = range(0, 18, length = 30)
# Construct `Subdivision Based Grid`
grid = subdivision_based_grid(ds, (xg, yg))
grid.lvl_array
```
The constructed array corresponds to levels of discretization for specific regions of the grid as a powers of 2,
meaning that if area index is assigned to be `3`, for example, the algorithm will treat the region as one being
`2^3 = 8` times more dense than originally user provided grid `(xg, yg)`.

Now upon the construction of this structure, one can simply pass it into mapper function as usual.

```@example MAIN
fig = Figure()
ax = Axis(fig[1,1])
# passing SubdivisionBasedGrid into mapper
mapper = AttractorsViaRecurrences(ds, grid;
        Dt = 0.1, sparse = true,
        consecutive_recurrences = 10, attractor_locate_steps = 10,
        maximum_iterations = 1000,
    )

# Find attractor and its fraction (fraction is always 1 here)
sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
attractors_SBD = extract_attractors(mapper)
scatter!(ax, vec(attractors_SBD[1]); label = "SubdivisionBasedGrid")


# to compare the results we also construct RegularGrid of same length here
xg = yg = range(0, 18, length = 30)
mapper = AttractorsViaRecurrences(ds, (xg, yg);
        Dt = 0.1, sparse = true,
        consecutive_recurrences = 10, attractor_locate_steps = 10,
        maximum_iterations = 1000,
    )

sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
attractors_reg = extract_attractors(mapper)
scatter!(ax, vec(attractors_reg[1]); label = "RegularGrid")

axislegend(ax)
fig

```

## Basin fractions continuation in the magnetic pendulum

Perhaps the simplest application of [`continuation`](@ref) is to produce a plot of how the fractions of attractors change as we continuously change the parameter we changed above to calculate tipping probabilities.

### Computing the fractions

This is what the following code does:

```@example MAIN
# initialize projected magnetic pendulum
using Attractors, PredefinedDynamicalSystems
using Random: Xoshiro
ds = Systems.magnetic_pendulum(; d = 0.3, α = 0.2, ω = 0.5)
xg = yg = range(-3, 3; length = 101)
ds = ProjectedDynamicalSystem(ds, 1:2, [0.0, 0.0])
# Choose a mapper via recurrences
mapper = AttractorsViaRecurrences(ds, (xg, yg); Δt = 1.0)
# What parameter to change, over what range
γγ = range(1, 0; length = 101)
prange = [[1, 1, γ] for γ in γγ]
pidx = :γs
# important to make a sampler that respects the symmetry of the system
region = HSphere(3.0, 2)
sampler, = statespace_sampler(region, 1234)
# continue attractors and basins:
# `Inf` threshold fits here, as attractors move smoothly in parameter space
rsc = RecurrencesFindAndMatch(mapper; threshold = Inf)
fractions_curves, attractors_info = continuation(
    rsc, prange, pidx, sampler;
    show_progress = false, samples_per_parameter = 100
)
# Show some characteristic fractions:
fractions_curves[[1, 50, 101]]
```


### Plotting the fractions
We visualize them using a predefined function that you can find in `docs/basins_plotting.jl`

```@example MAIN
# careful; `prange` isn't a vector of reals!
plot_basins_curves(fractions_curves, γγ)
```


### Fixed point curves

A by-product of the analysis is that we can obtain the curves of the position of fixed points for free. However, only the stable branches can be obtained!

```@example MAIN
using CairoMakie
fig = Figure()
ax = Axis(fig[1,1]; xlabel = L"\gamma_3", ylabel = "fixed point")
# choose how to go from attractor to real number representation
function real_number_repr(attractor)
    p = attractor[1]
    return (p[1] + p[2])/2
end

for (i, γ) in enumerate(γγ)
    for (k, attractor) in attractors_info[i]
        scatter!(ax, γ, real_number_repr(attractor); color = Cycled(k))
    end
end
fig
```

as you can see, two of the three fixed points, and their stability, do not depend at all on the parameter value, since this parameter value tunes the magnetic strength of only the third magnet. Nevertheless, the **fractions of basin of attraction** of all attractors depend strongly on the parameter. This is a simple example that highlights excellently how this new approach we propose here should be used even if one has already done a standard linearized bifurcation analysis.


## Extinction of a species in a multistable competition model

In this advanced example we utilize both [`RecurrencesFindAndMatch`](@ref) and [`aggregate_attractor_fractions`](@ref) in analyzing species extinction in a dynamical model of competition between multiple species.
The final goal is to show the percentage of how much of the state space leads to the extinction or not of a pre-determined species, as we vary a parameter. The model however displays extreme multistability, a feature we want to measure and preserve before aggregating information into "extinct or not".

To measure and preserve this we will apply [`RecurrencesFindAndMatch`](@ref) as-is first. Then we can aggregate information. First we have
```julia
using Attractors, OrdinaryDiffEq
using PredefinedDynamicalSystems
using Random: Xoshiro
# arguments to algorithms
samples_per_parameter = 1000
total_parameter_values = 101
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9, maxiters = Inf)
recurrences_kwargs = (; Δt= 1.0, consecutive_recurrences=9, diffeq);
# initialize dynamical system and sampler
ds = PredefinedDynamicalSystems.multispecies_competition() # 8-dimensional
ds = CoupledODEs(ODEProblem(ds), diffeq)
# define grid in state space
xg = range(0, 60; length = 300)
grid = ntuple(x -> xg, 8)
prange = range(0.2, 0.3; length = total_parameter_values)
pidx = :D
sampler, = statespace_sampler(grid, 1234)
# initialize mapper
mapper = AttractorsViaRecurrences(ds, grid; recurrences_kwargs...)
# perform continuation of attractors and their basins
continuation = RecurrencesFindAndMatch(mapper; threshold = Inf)
fractions_curves, attractors_info = continuation(
    continuation, prange, pidx, sampler;
    show_progress = true, samples_per_parameter
);
plot_basins_curves(fractions_curves, prange; separatorwidth = 1)
```

![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/multispecies_competition_fractions.png)

_this example is not actually run when building the docs, because it takes about 60 minutes to complete depending on the computer; we load precomputed results instead_

As you can see, the system has extreme multistability with 64 unique attractors
(according to the default matching behavior in [`RecurrencesFindAndMatch`](@ref); a stricter matching with less than `Inf` threshold would generate more "distinct" attractors).
One could also isolate a specific parameter slice, and do the same as what we do in
the [Fractality of 2D basins of the (4D) magnetic pendulum](@ref) example, to prove that the basin boundaries are fractal, thereby indeed confirming the paper title "Fundamental Unpredictability".

Regardless, we now want to continue our analysis to provide a figure similar to the
above but only with two colors: fractions of attractors where a species is extinct or not. Here's how:

```julia
species = 3 # species we care about its existence

featurizer = (A) -> begin
    i = isextinct(A, species)
    return SVector(Int32(i))
end
isextinct(A, idx = unitidxs) = all(a -> a <= 1e-2, A[:, idx])

# `minneighbors = 1` is crucial for grouping single attractors
groupingconfig = GroupViaClustering(; min_neighbors=1, optimal_radius_method=0.5)

aggregated_fractions, aggregated_info = aggregate_attractor_fractions(
    fractions_curves, attractors_info, featurizer, groupingconfig
)

plot_basins_curves(aggregated_fractions, prange;
    separatorwidth = 1, colors = ["green", "black"],
    labels = Dict(1 => "extinct", 2 => "alive"),
)
```

![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/multispecies_competition_fractions_aggr.png)

(in hindsight, the labels are reversed; attractor 1 is the alive one, but oh well)

## Trivial featurizing and grouping for basins fractions

This is a rather trivial example showcasing the usage of [`AttractorsViaFeaturizing`](@ref). Let us use once again the magnetic pendulum example. For it, we have a really good idea of what features will uniquely describe each attractor: the last points of a trajectory (which should be very close to the magnetic the trajectory converged to). To provide this information to the [`AttractorsviaFeaturizing`](@ref) we just create a julia function that returns this last point

```@example MAIN
using Attractors
using PredefinedDynamicalSystems

ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])

function featurizer(X, t)
    return X[end]
end

mapper = AttractorsViaFeaturizing(psys, featurizer; Ttr = 200, T = 1)

xg = yg = range(-4, 4; length = 101)

region = HRectangle([-4, 4], [4, 4])
sampler, = statespace_sampler(region)

fs = basins_fractions(mapper, sampler; show_progress = false)
```
As expected, the fractions are each about 1/3 due to the system symmetry.


## Featurizing and grouping across parameters (MCBB)
Here we showcase the example of the Monte Carlo Basin Bifurcation publication.
For this, we will use [`FeaturizeGroupAcrossParameter`](@ref) while also providing a `par_weight = 1` keyword.
However, we will not use a network of 2nd order Kuramoto oscillators (as done in the paper by Gelbrecht et al.) because it is too costly to run on CI.
Instead, we will use "dummy" system which we know analytically the attractors and how they behave versus a parameter.

 the Henon map and try to group attractors into period 1 (fixed point), period 3, and divergence to infinity. We will also use a pre-determined optimal radius for clustering, as we know a-priory the expected distances of features in feature space (due to the contrived form of the `featurizer` function below).

```@example MAIN
using Attractors, Random

function dumb_map(dz, z, p, n)
    x, y = z
    r = p[1]
    if r < 0.5
        dz[1] = dz[2] = 0.0
    else
        if x > 0
            dz[1] = r
            dz[2] = r
        else
            dz[1] = -r
            dz[2] = -r
        end
    end
    return
end

r = 3.833
ds = DiscreteDynamicalSystem(dumb_map, [0., 0.], [r])
```


```@example MAIN
sampler, = statespace_sampler(HRectangle([-3.0, -3.0], [3.0, 3.0]), 1234)

rrange = range(0, 2; length = 21)
ridx = 1

featurizer(a, t) = a[end]
clusterspecs = GroupViaClustering(optimal_radius_method = "silhouettes", max_used_features = 200)
mapper = AttractorsViaFeaturizing(ds, featurizer, clusterspecs; T = 20, threaded = true)
gap = FeaturizeGroupAcrossParameter(mapper; par_weight = 1.0)
fractions_curves, clusters_info = continuation(
    gap, rrange, ridx, sampler; show_progress = false
)
fractions_curves
```

Looking at the information of the "attractors" (here the clusters of the grouping procedure) does not make it clear which label corresponds to which kind of attractor, but we can look at the:

```@example MAIN
clusters_info
```

## Using histograms and histogram distances as features

One of the aspects discussed in the original MCBB paper and implementation was the usage of histograms of the means of the variables of a dynamical system as the feature vector. This is useful in very high dimensional systems, such as oscillator networks, where the histogram of the means is significantly different in synchronized or unsychronized states.

This is possible to do with current interface without any modifications, by using two more packages: ComplexityMeasures.jl to compute histograms, and Distances.jl for the Kullback-Leibler divergence (or any other measure of distance in the space of probability distributions you fancy).

The only code we need to write to achieve this feature is a custom featurizer and providing an alternative distance to `GroupViaClustering`. The code would look like this:

```julia
using Distances: KLDivergence
using ComplexityMeasures: ValueHistogram, FixedRectangularBinning, probabilities

# you decide the binning for the histogram, but for a valid estimation of
# distances, all histograms must have exactly the same bins, and hence be
# computed with fixed ranges, i.e., using the `FixedRectangularBinning`
const binning = FixedRectangularBinning(range(-5, 5; length = 11))

function histogram_featurizer(A, t)
    ms = mean.(columns(A)) # vector of mean of each variable
    p = probabilities(ValueHistogram(binning), ms) # this is the histogram
    return vec(p) # because Distances.jl doesn't know `Probabilities`
end

gconfig = GroupViaClustering(;
    clust_distance_metric = KLDivergence(), # or any other PDF distance
)
```

You can then pass the `histogram_featurizer` and `gconfig` to an [`AttractorsViaFeaturizing`](@ref) and use the rest of the library as usual.


## Animation illustrating `AttractorsViaRecurrences`

The following Julia script inputs a 2D continuous time dynamical system and animates its time evolution while illustrating how [`AttractorsViaRecurrences`](@ref) works.

```julia
using Attractors, CairoMakie
using PredefinedDynamicalSystems
using OrdinaryDiffEq

# Set up dynamical system: bi-stable predator pray
function predator_prey_rule(u, p, t)
    r, c, μ, ν, α, β, χ, δ = p
    N, P = u
    common = α*N*P/(β+N)
    dN = r*N*(1 - (c/r)*N)*((N-μ)/(N+ν)) - common
    dP = χ*common - δ*P
    return SVector(dN, dP)
end

u0 = SVector(8.0, 0.01)
r = 2.0
# r, c, μ, ν, α, β, χ, δ = p
p = [r, 0.19, 0.03, 0.003, 800, 1.5, 0.004, 2.2]

diffeq = (alg = Rodas5P(), abstol = 1e-9, rtol = 1e-9)
ds = CoupledODEs(predator_prey_rule, u0, p; diffeq)

u0s = [ # animation will start from these initial conditions
    [10, 0.012],
    [15, 0.02],
    [12, 0.01],
    [13, 0.015],
    [5, 0.02],
]

density = 31
xg = range(-0.1, 20; length = density)
yg = range(-0.001, 0.03; length = density)
Δt = 0.1
grid = (xg, yg)
mapper = AttractorsViaRecurrences(ds, grid;
    Δt, consecutive_attractor_steps = 10, consecutive_basin_steps = 10, sparse = false,
    consecutive_recurrences = 100, attractor_locate_steps = 100,
)

##########################################################################

function animate_attractors_via_recurrences(
        mapper::AttractorsViaRecurrences, u0s;
        colors = ["#FFFFFF", "#7143E0","#0A9A84","#AF9327","#791457", "#6C768C", "#4287f5",],
        filename = "recurrence_algorithm.mp4",
    )

    grid_nfo = mapper.bsn_nfo.grid_nfo

    fig = Figure()
    ax = Axis(fig[1,1])

    # Populate the grid with poly! rectangle plots. However! The rectangles
    # correspond to the same "cells" of the grid. Additionally, all
    # rectangles are colored with an _observable_, that can be accessed
    # later using the `basin_cell_index` function. The observable
    # holds the face color of the rectangle!

    # Only 6 colors; need 3 for base, and extra 2 for each attractor.
    # will choose initial conditions that are only in the first 2 attractors
    COLORS = map(c -> Makie.RGBA(Makie.RGB(to_color(c)), 0.9), colors)

    function initialize_cells2!(ax, grid; kwargs...)
        # These are all possible outputs of the `basin_cell_index` function
        idxs = all_cartesian_idxs(grid)
        color_obs = Matrix{Any}(undef, size(idxs)...)
        # We now need to reverse-engineer
        for i in idxs
            rect = cell_index_to_rect(i, grid)
            color = Observable(COLORS[1])
            color_obs[i] = color
            poly!(ax, rect; color = color, strokecolor = :black, strokewidth = 0.5)
        end
        # Set the axis limits better
        mini, maxi = Attractors.minmax_grid_extent(grid)
        xlims!(ax, mini[1], maxi[1])
        ylims!(ax, mini[2], maxi[2])
        return color_obs
    end

    all_cartesian_idxs(grid::Attractors.RegularGrid) = CartesianIndices(length.(grid.grid))

    # Given a cartesian index, the output of `basin_cell_index`, create
    # a `Rect` object that corresponds to that grid cell!
    function cell_index_to_rect(n::CartesianIndex, grid::Attractors.RegularGrid)
        x = grid.grid[1][n[1]]
        y = grid.grid[2][n[2]]
        dx = grid.grid_steps[1]
        dy = grid.grid_steps[2]
        rect = Rect(x - dx/2, y - dy/2, dx, dy)
        return rect
    end

    color_obs = initialize_cells2!(ax, grid_nfo)

    # plot the trajectory
    state2marker = Dict(
        :att_search => :circle,
        :att_found => :dtriangle,
        :att_hit => :rect,
        :lost => :star5,
        :bas_hit => :xcross,
    )

    # This function gives correct color to search, recurrence, and
    # the individual attractors. Ignores the lost state.
    function update_current_cell_color!(cellcolor, bsn_nfo)
        # We only alter the cell color at specific situations
        state = bsn_nfo.state
        if state == :att_search
            if cellcolor[] == COLORS[1] # empty
                cellcolor[] = COLORS[2] # visited
            elseif cellcolor[] == COLORS[2] # visited
                cellcolor[] = COLORS[3] # recurrence
            end
        elseif state == :att_found
            attidx = (bsn_nfo.current_att_label ÷ 2)
            attlabel = (attidx - 1)*2 + 1
            cellcolor[] = COLORS[3+attlabel]
        end
        return
    end

    # Iteration and labelling
    ds = mapper.ds
    bsn_nfo = mapper.bsn_nfo
    u0 = current_state(ds)

    traj = Observable(SVector{2, Float64}[u0])
    point = Observable([u0])

    marker = Observable(:circle)
    lines!(ax, traj; color = :black, linewidth = 1)
    scatter!(ax, point; color = (:black, 0.5), markersize = 20, marker, strokewidth = 1.0, strokecolor = :black)

    stateobs = Observable(:att_search)
    consecutiveobs = Observable(0)
    labeltext = @lift("state: $($(stateobs))\nconsecutive: $($(consecutiveobs))")

    Label(fig[0, 1][1,1], labeltext; justification = :left, halign = :left, tellwidth = false)
    # add text  with options
    kwargstext = prod("$(p[1])=$(p[2])\n" for p in mapper.kwargs)
    Label(fig[0, 1][1, 2], kwargstext; justification = :right, halign = :right, tellwidth = false)


    # make legend
    entries = [PolyElement(color = c) for c in COLORS[2:end]]
    labels = ["visited", "recurrence", "attr. 1", "basin 1", "attr. 2", "basin 2"]

    Legend(fig[:, 2][1, 1], entries, labels)


    # %% loop
    # The following code is similar to the source code of `recurrences_map_to_label!`

    cell_label = 0
    record(fig, filename) do io
        for u0 in u0s
            reinit!(ds, copy(u0))
            traj[] = [copy(u0)]

            while cell_label == 0
                step!(ds, bsn_nfo.Δt)
                u = current_state(ds)

                # update FSM
                n = Attractors.basin_cell_index(u, bsn_nfo.grid_nfo)
                cell_label = Attractors.finite_state_machine!(bsn_nfo, n, u; mapper.kwargs...)

                state = bsn_nfo.state

                if cell_label ≠ 0 # FSM terminated; we assume no lost/divergence in the system
                    stateobs[] = :terminated

                    # color-code initial condition if we converged to attractor
                    # or to basin (even or odd cell label)
                    u0n = Attractors.basin_cell_index(u0, bsn_nfo.grid_nfo)

                    basidx = (cell_label - 1)
                    color_obs[u0n][] = COLORS[3+basidx]

                    # Clean up: all "visited" cells become white again
                    visited_idxs = findall(v -> (v[] == COLORS[2] || v[] == COLORS[3]), color_obs)
                    for n in visited_idxs
                        color_obs[n][] = COLORS[1] # empty
                    end
                    # clean up trajectory line
                    traj[] = []

                    for i in 1:15; recordframe!(io); end
                    cell_label = 0
                    break
                end

                # update visuals:
                point[] = [u]
                push!(traj[], u)
                notify(traj)
                marker[] = state2marker[state]
                stateobs[] = state
                consecutiveobs[] = bsn_nfo.consecutive_match

                update_current_cell_color!(color_obs[n], bsn_nfo)

                recordframe!(io)
            end
        end
    end
end

animate_attractors_via_recurrences(mapper, u0s)
```

```@raw html
<video width="75%" height="auto" controls autoplay loop>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/recurrence_algorithm.mp4?raw=true" type="video/mp4">
</video>
```
