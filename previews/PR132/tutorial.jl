# # [Attractors.jl Tutorial](@id tutorial)

# [`Attractors`](@ref) is a submodule of the **DynamicalSystems.jl** library.
# This tutorial will walk you through its main functionality.
# That is, given a `DynamicalSystem` instance, find all its attractors and their basins
# of attraction. Then,
# continue these attractors, and their stability properties, across a parameter value.
# It also offers various functions that compute nonlocal stability properties for an
# attractor, any of which can be used in the continuation to quantify stability.

# Besides this main functionality, there are plenty of other stuff,
# like for example [`edgestate`](@ref) or [`basins_fractal_dimension`](@ref),
# but we won't cover anything else in this introductory tutorial.
# See the [examples](@ref examples) page instead.

# ## Input: a `DynamicalSystem`

# The key input for most functionality of Attractors.jl is an instance of
# a `DynamicalSystem`. If you don't know how to make
# a `DynamicalSystem`, you need to consult the main tutorial of the
# [DynamicalSystems.jl library](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/tutorial/).
# For this tutorial we will use a modified Lorenz-like system with equations

# ```math
# \begin{align*}
# \dot{x} & = y - x \\
# \dot{y}  &= -x*z + b*|z| \\
# \dot{z}  &= x*y - a \\
# \end{align*}
# ```

# which we define in code as
using Attractors # part of `DynamicalSystems`, so it re-exports functionality for making them!
using OrdinaryDiffEq # for accessing advanced ODE Solvers

function modified_lorenz_rule(u, p, t)
    x, y, z = u; a, b = p
    dx = y - x
    dy = - x*z + b*abs(z)
    dz = x*y - a
    return SVector(dx, dy, dz)
end

p0 = [5.0, 0.1] # parameters
u0 = [-4.0, 5, 0] # state
ds = CoupledODEs(modified_lorenz_rule, u0, p0;
    diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9, dt = 0.01)
)

# ## Finding attractors

# There are two major methods for finding attractors in dynamical systems.
# Explanation of how they work is in their respective docs.

# 1. [`AttractorsViaRecurrences`](@ref).
# 2. [`AttractorsViaFeaturizing`](@ref).

# You can consult [Datseris2023](@cite) for a comparison between the two.

# As far as the user is concerned, both algorithms are part of the same interface,
# and can be used in the same way. The interface is extendable as well,
# and works as follows.

# First, we create an instance of such an "attractor finding algorithm",
# which we call `AttractorMapper`. For example, [`AttractorsViaRecurrences`](@ref)
# requires a tesselated grid of the state space to search for attractors in.
# It also allows the user to tune some meta parameters, but in our example
# they are already tuned for the dynamical system at hand. So we initialize

grid = (
    range(-10.0, 10.0; length = 100), # x
    range(-15.0, 15.0; length = 100), # y
    range(-15.0, 15.0; length = 100), # z
)

mapper = AttractorsViaRecurrences(ds, grid;
    consecutive_recurrences = 1000, attractor_locate_steps = 1000,
    consecutive_lost_steps = 100,
)

# This `mapper` can map any initial condition to the corresponding
# attractor ID, for example

mapper([-4.0, 5, 0])

# while

mapper([4.0, 2, 0])

# the fact that these two different initial conditions got assigned different IDs means
# that they converged to a different attractor.
# The attractors are stored in the mapper internally, to obtain them we
# use the function

attractors = extract_attractors(mapper)

# In Attractors.jl, all information regarding attractors is always a standard Julia
# `Dict`, which maps attractor IDs (positive integers) to the corresponding quantity.
# Here the quantity are the attractors themselves, represented as `StateSpaceSet`.

# Let's visualize them
using CairoMakie
fig = Figure()
ax = Axis(fig[1,1]; title = "bistable lorenz-like")
for (k, A) in attractors
    scatter!(ax, A[:, 1], A[:, 2]; label = "ID = $(k)")
end
axislegend(ax; position = :lt)
fig

# In our example system we see that for the chosen parameters there are two coexisting attractors:
# a limit cycle and a chaotic attractor.
# There may be more attractors though! We've only checked two initial conditions,
# so we could have found at most two attractors!
# However, it can get tedious to manually iterate over initial conditions, which is why
# this `mapper` is typically given to higher level functions for finding attractors
# and their basins of attraction. The simplest one
# is [`basins_fractions`](@ref). Using the `mapper`,
# it finds "all" attractors of the dynamical system and reports the state space fraction
# each attractors attracts. The search is probabilistic, so "all" attractors means those
# that at least one initial condition converged to.

# We can provide explicitly initial conditions to [`basins_fraction`](@ref),
# however it is typically simpler to provide it with with a state space sampler instead:
# a function that generates random initial conditions in the region of the
# state space that we are interested in. Here this region coincides with `grid`,
# so we can simply do:

sampler, = statespace_sampler(grid)

sampler() # random i.c.

#

sampler() # another random i.c.

# and finally call

fs = basins_fractions(mapper, sampler)

# The returned `fs` is a dictionary mapping each attractor ID to
# the fraction of the state space the corresponding basin occupies.
# With this we can confirm that there are (likely) only two attractors
# and that both attractors are robust as both have sufficiently large basin fractions.

# To obtain the full basins, which is computationally much more expensive,
# use [`basins_of_attraction`](@ref).

# You can use alternative algorithms in [`basins_fractions`](@ref), see
# the documentation of [`AttractorMapper`](@ref) for possible subtypes.
# [`AttractorMapper`](@ref) defines an extendable interface and can be enriched
# with other methods in the future!

# ## Global continuation

# If you have heard before the word "continuation", then you are likely aware of the
# **traditional continuation-based bifurcation analysis (CBA)** offered by many software,
# such as AUTO, MatCont, and in Julia [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).
# Here we offer a completely different kind of continuation called **global continuation**.

# A direct comparison of the two approaches is not truly possible, because they do different things.
# The traditional continuation analysis continues the curves of individual fixed
# points across the joint state-parameter space and tracks their _local (linear) stability_.
# The global continuation in Attractors.jl finds all attractors, including chaotic ones,
# in the whole of the state space (that it searches in), and continues all of these attractors
# concurrently along a parameter axis.
# Additionally, this global continuation tracks a _nonlocal_ stability property which by
# default is the basin fraction.

# This is a fundamental difference. Because all attractors are simultaneously
# tracked across the parameter axis, the user may arbitrarily estimate _any_
# property of the attractors and how it varies as the parameter varies.
# A more detailed comparison between these two approaches can be found in [Datseris2023](@cite).

# To perform the continuation is extremely simple. First, we decide what parameter,
# and what range, to continue over:

prange = 4.7:0.02:6
pidx = 1 # index of the parameter

# Then, we may call the [`global_continuation`](@ref) function.
# We have to provide a continuation algorithm, which itself references an [`AttractorMapper`](@ref).
# In this example we will re-use the `mapper` to create a [`RecurrencesFindAndMatch`](@ref) continuation algorithm.
# This algorithm uses the `mapper` to find all attractors at each parameter value.
# Then, it performs a "matching" step, ensuring a "continuity" of the attractor
# label across the parameter axis. You can read the docstring for more details,
# as this algorithm is quite sophisticated!

# For now we can use all of its default options which are reliable most of the time

rafm = RecurrencesFindAndMatch(mapper)

# and call

fractions_cont, attractors_cont = global_continuation(
	rafm, prange, pidx, sampler; samples_per_parameter = 1_000
)

# the output is given as two vectors. Each vector is a dictionary
# mapping attractor IDs to their basin fractions, or their state space sets, respectively.
# Both vectors have the same size as the parameter range.
# For example, the attractors at the 34-th parameter value are:

attractors_cont[34]

# There is a fantastic convenience function for animating
# the attractors evolution, that utilizes things we have
# already defined:

animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, prange, pidx;
);

# ```@raw html
# <video width="auto" controls autoplay loop>
# <source src="../attracont.mp4" type="video/mp4">
# </video>
# ```

# Hah, how cool is that! The attractors pop in and out of existence like out of nowhere!
# It would be incredibly difficult to find these attractors in traditional continuation software
# where a rough estimate of the period is required! (It would also be too hard due to the presence
# of chaos for most of the parameter values, but that's another issue!)

# Now typically a continuation is visualized in a 2D plot where the x axis is the
# parameter axis. We can do this with the convenience function:

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), prange,
)

# In the top panel are the basin fractions, by default plotted as stacked bars.
# Bottom panel is a visualization of the tracked attractors.
# The argument `A -> minimum(A[:, 1])` is simply a function that maps
# an attractor into a real number for plotting.
# We can provide more functions to visualize other aspects of the attractors:

a2rs = [
    A -> minimum(A[:, 1]),
    A -> log(length(A)), # proxy for "complexity"
]

fig = plot_basins_attractors_curves(
	fractions_cont, attractors_cont, a2rs, prange; add_legend = false
)

ax1, ax2 = content.((fig[2,1], fig[3,1]))

ax1.ylabel = "min(A₁)"
ax2.ylabel = "log(len(A))"

fig

# ## Enhancing the continuation

# The biggest strength of Attractors.jl is that it is not an isolated software.
# It is part of **DynamicalSystems.jl**. Here, we will use the full power of
# **DynamicalSystems.jl** and enrich the above continuation with various other
# measures of nonlocal stability, in particular Lyapunov exponents and
# the minimal fatal shock.

# First, let's estimate the maximum Lyapunov exponent (MLE) for all attractors,
# using the `lyapunovspectrum` function that comes from the ChaosTools.jl submodule.

using ChaosTools: lyapunov

lis = map(enumerate(prange)) do (i, p) # loop over parameters
    set_parameter!(ds, pidx, p) # important! We use the dynamical system!
    attractors = attractors_cont[i]
    Dict(k => lyapunov(ds, 2000.0; u0 = A[1]) for (k, A) in attractors)
end

# The above `map` loop may be intimidating if you are a beginner, but it is
# really just a shorter way to write a `for` loop for our example.
# We iterate over all parameters, and for each we first update the dynamical
# system with the correct parameter, and then extract the MLE
# for each attractor. `map` just means that we don't have to pre-allocate a
# new vector before the loop; it creates it for us.

# We can visualize the LE with the other convenience function [`plot_continuation_curves!`](@ref),

ax3 = Axis(fig[4, 1]; ylabel = "MLE")
plot_continuation_curves!(ax3, lis, prange; add_legend = false)

fig

# This reveals crucial information for tha attractors, whether they are chaotic or not, that we would otherwise obtain only by visualizing the system dynamics at every single parameter.
# The story we can see now is that the dynamics start with a limit cycle (0 Lyapunov exponent), go into bi-stability of chaos and limit cycle, then there is only one limit cycle again, and then a chaotic attractor appears again, for a second bistable regime.

# The last piece of information to add is yet another measure of nonlocal stability: the minimal fatal shock (MFS), which is provided by [`minimal_fatal_shock`](@ref).
# The code to estimate this is similar with the `map` block for the MLE.
# Here however we re-use the created `mapper`, but now we must not forget to reset it inbetween parameter increments:

using LinearAlgebra: norm
search_area = collect(extrema.(grid ./ 2)) # smaller search = faster results
search_algorithm = MFSBlackBoxOptim(max_steps = 1000, guess = ones(3))

mfss = map(enumerate(prange)) do (i, p)
    set_parameter!(ds, pidx, p)
    reset_mapper!(mapper) # reset so that we don't have to re-initialize
    ## We need a special clause here: if there is only 1 attractor,
    ## then there is no MFS. It is undefined. We set it to `NaN`,
    ## which conveniently, will result to nothing being plotted by Makie.
    attractors = attractors_cont[i]
    if length(attractors) == 1
        return Dict(k => NaN for (k, A) in attractors)
    end
    ## otherwise, compute the actual MFS from the first point of each attractor
    Dict(k =>
        norm(minimal_fatal_shock(mapper, A[1], search_area, search_algorithm))
        for (k, A) in attractors
    )
end

# In a real application we wouldn't use the first point of each attractor,
# as the first point is completely random on the attractor (at least, for the
# [`AttractorsViaRecurrences`] mapper we use here).
# We would do this by examining the whole `A` object in the above block
# instead of just using `A[1]`. But this is a tutorial so we don't care!

# Right, so now we can visualize the MFS with the rest of the other quantities:

ax4 = Axis(fig[5, 1]; ylabel = "MFS", xlabel = "parameter")
plot_continuation_curves!(ax4, mfss, prange; add_legend = false)

## make the figure prettier
for ax in (ax1, ax2, ax3); hidexdecorations!(ax; grid = false); end
resize!(fig, 500, 600)
fig

# And that's the end of the tutorial! See the [examples](@ref examples) for
# more runnable code, and see the [API](@ref) for a list of all functions and algorithms!