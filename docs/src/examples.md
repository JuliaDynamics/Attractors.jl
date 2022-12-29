# Examples for Attractors.jl

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
mapper = AttractorsViaRecurrences(ds, (xg, yg); sparse = false)
basins, attractors = basins_of_attraction(mapper; show_progress = false)
basins
```
```@example MAIN
attractors
```

Now let's plot this as a heatmap
```@example MAIN
using CairoMakie
# Set up some code for plotting attractors
function scatter_attractors!(ax, attractors)
    for k ∈ keys(attractors)
        x, y = columns(attractors[k])
        scatter!(ax, vec(attractors[k]);
            color = Cycled(k),
            strokewidth = 3, strokecolor = :white
        )
    end
end

generate_cmap(n) = cgrad(Main.COLORS[1:n], n; categorical = true)
ids = sort!(unique(basins))
cmap = generate_cmap(length(ids))

fig, ax = heatmap(xg, yg, basins;
    colormap = cmap, colorrange = (ids[1] - 0.5, ids[end]+0.5),
)
scatter_attractors!(ax, attractors)
fig
```


## Fractality of 2D basins of the (4D) magnetic pendulum
In this section we will calculate the basins of attraction of the four-dimensional magnetic pendulum. We know that the attractors of this system are all individual fixed points on the (x, y) plane so we will only compute the basins there. We can also use this opportunity to highlight a different method, the [`AttractorsViaProximity`](@ref) which works when we already know where the attractors are. Furthermore we will also use a `projected_integrator` to project the 4D system onto a 2D plane, saving a lot of computational time!

### Computing the basins

```@example MAIN
ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)
psys = projected_integrator(ds, [1, 2], [0.0, 0.0])
```

For this systems we know the attractors are close to the magnet positions, so we can just do
```@example MAIN
attractors = Dict(i => Dataset([ds.f.magnets[i]]) for i in 1:3)
mapper = AttractorsViaProximity(psys, attractors; sparse = false)
```

and as before, get the basins of attraction
```@example MAIN
xg = yg = range(-4, 4; length=150)
grid = (xg, yg)
basins, = basins_of_attraction(mapper, grid; show_progress = false)
ids = sort!(unique(basins))
cmap = generate_cmap(length(ids))
fig, ax = heatmap(xg, yg, basins;
    colormap = cmap, colorrange = (ids[1] - 0.5, ids[end]+0.5),
)
scatter_attractors!(ax, attractors)
fig
```

### Computing the uncertainty exponent
Let's now calculate the [`uncertainty_exponent`](@ref) for this system as well.
The calculation is straightforward:
```@example MAIN
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
ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3, γs = [1.0, 1.0, 0.1])
psys = projected_integrator(ds, [1, 2], [0.0, 0.0]; diffeq = (reltol = 1e-9,))
mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 1)
basins_after, attractors_after = basins_of_attraction(
    mapper, (xg, yg); show_progress = false
)
# matching attractors is important!
match_attractor_ids!(attractors_after, attractors)

# now plot
ids = sort!(unique(basins_after))
cmap = generate_cmap(length(ids))
fig, ax = heatmap(xg, yg, basins_after;
    colormap = cmap, colorrange = (ids[1] - 0.5, ids[end]+0.5),
)
scatter_attractors!(ax, attractors_after)
fig
```


```@example MAIN
P = tipping_probabilities(basins, basins_after)
```
As you can see `P` has size 3×2, as after the change only 2 attractors have been identified
in the system (3 still exist but our state space discretization isn't accurate enough to
find the 3rd because it has such a small basin).
Also, the first row of `P` is 50% probability to each other magnet, as it should be due to
the system's symmetry.


## Basin fractions continuation in the magnetic pendulum
Perhaps the simplest application of [`basins_fractions_continuation`](@ref) is to produce a plot of how the fractions of attractors change as we continuously change the parameter we changed above to calculate tipping probabilities.



TODO: Write it.


## 3D basins via recurrences
To showcase the true power of [`AttractorsViaRecurrences`](@ref) we need to use a system whose attractors span higher-dimensional space. An example is
```@example MAIN
ds = Systems.thomas_cyclical(b = 0.1665)
```
which, for this parameter, contains 5 coexisting attractors. 3 of them are entangled
periodic orbits that span across all three dimensions, and the remaining 2 are fixed points.

To compute the basins we define a three-dimensional grid and call on it
[`basins_of_attraction`](@ref).

```julia
# This computation takes about an hour
xg = yg = zg = range(-6.0, 6.0; length = 251)
mapper = AttractorsViaRecurrences(ds, (xg, yg, zg))
basins, attractors = basins_of_attraction(mapper)
attractors
```
```
Dict{Int16, Dataset{3, Float64}} with 5 entries:
  5 => 3-dimensional Dataset{Float64} with 1 points
  4 => 3-dimensional Dataset{Float64} with 379 points
  6 => 3-dimensional Dataset{Float64} with 1 points
  2 => 3-dimensional Dataset{Float64} with 538 points
  3 => 3-dimensional Dataset{Float64} with 537 points
  1 => 3-dimensional Dataset{Float64} with 1 points
```

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



## Extinction of a species in a multistable competition model

In this advanced example we utilize both [`RecurrencesSeededContinuation`](@ref) and [`aggregate_attractor_fractions`](@ref) in analyzing species extinction in a dynamical model of competition between multiple species.
The final goal is to show the percentage of how much of the state space leads to the extinction or not of a pre-determined species, as we vary a parameter. The model however displays extreme multistability, a feature we want to measure and preserve before aggregating information into "extinct or not".

To measure and preserve this we will apply [`RecurrencesSeededContinuation`](@ref) as-is first. Then we can aggregate information. First we have
```julia
using Attractors, OrdinaryDiffEq
using Random: Xoshiro
# arguments to algorithms
samples_per_parameter = 1000
total_parameter_values = 101
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9, maxiters = Inf)
recurrences_kwargs = (; Δt= 1.0, mx_chk_fnd_att=9, diffeq);
# initialize dynamical systerm and sampler
ds = Systems.multispecies_competition() # 8-dimensional by default
xg = range(0, 60; length = 300)
grid = ntuple(x -> xg, 8)
prange = range(0.2, 0.3; length = total_parameter_values)
pidx = :D
sampler, = statespace_sampler(Xoshiro(1234);
    min_bounds = minimum.(grid), max_bounds = maximum.(grid)
)
# initialize mapper
mapper = AttractorsViaRecurrences(ds, grid; recurrences_kwargs...)
# perform continuation of attractors and their basins
# NOTE: in a realistic applicaiton the threshold should not be infinite,
# but some sensible distance in state space units
continuation = RecurrencesSeedingContinuation(mapper; threshold = Inf)
fractions_curves, attractors_info = basins_fractions_continuation(
    continuation, prange, pidx, sampler;
    show_progress = true, samples_per_parameter
);
Main.basins_fractions_plot(fractions_curves, prange; separatorwidth = 1)
```

![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/multispecies_competition_fractions.png)

_this  example is not actually run when building the docs, because it takes about 30 minutes to complete depending on the computer; we load precomputed results instead_

As you can see, the system has extreme multistability with 64 unique attractors
(according to the default matching behavior in [`RecurrencesSeededContinuation`](@ref); a stricter matching with less than `Inf` threshold would generate more "distinct" attractors).
One could also isolate a specific parameter slice, and to the same as what we do in
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

basins_fractions_plot(aggregated_fractions, prange; separatorwidth = 1,
labels = Dict(1 => "extinct", 2 => "alive"), colors = ["green", "black"])
```

![](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/attractors/multispecies_competition_fractions_aggr.png)

(in hindsight, the labels are reversed; attractor 1 is the alive one)


## An example with attractors via featurizing

TODO: example of using AttractorsViaFeaturizing and GroupAcrossParameterContinuation.
