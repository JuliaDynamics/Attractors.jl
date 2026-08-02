# # [Attractors.jl Tutorial](@id tutorial)

# ```@raw html
# <video width="auto" controls loop>
# <source src="../attracont.mp4" type="video/mp4">
# </video>
# ```

# [`Attractors`](@ref) is a component of the **DynamicalSystems.jl** library.
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

# ### Package versions used

import Pkg

#nb # Activate an environment in the folder containing the notebook
#nb Pkg.activate(dirname(@__DIR__))
#nb Pkg.add(["DynamicalSystems", "CairoMakie", "GLMakie", "OrdinaryDiffEqVerner"])
Pkg.status(["Attractors", "CairoMakie", "OrdinaryDiffEqVerner"])

#nb # ## Attractors.jl summary

#nb using Attractors # re=exported by `DynamicalSystems`
#nb @doc Attractors

# ## Tutorial - copy-pasteable version

# _Gotta go fast!_

# ```julia
# using Attractors, CairoMakie, OrdinaryDiffEqVerner
# ## Define key input: a `DynamicalSystem`
# function modified_lorenz_rule(u, p, t)
#     x, y, z = u; a, b = p
#     dx = y - x
#     dy = - x*z + b*abs(z)
#     dz = x*y - a
#     return SVector(dx, dy, dz)
# end
# p0 = [5.0, 0.1] # parameters
# u0 = [-4.0, 5, 0] # state
# diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9, dt = 0.01) # solver options
# ds = CoupledODEs(modified_lorenz_rule, u0, p0; diffeq)

# ## Define key input: an `AttractorMaper` that finds
# ## attractors of a `DynamicalSystem`
# grid = (
#     range(-15.0, 15.0; length = 150), # x
#     range(-20.0, 20.0; length = 150), # y
#     range(-20.0, 20.0; length = 150), # z
# )
# mapper = AttractorsViaRecurrences(ds, grid;
#     consecutive_recurrences = 1000,
#     consecutive_lost_steps = 100,
# )

# ## Find attractors and their basins of attraction state space fraction
# ## by randomly sampling initial conditions in state sapce
# sampler, = statespace_sampler(grid)
# algo = AttractorSeedContinueMatch(mapper)
# fs = basins_fractions(mapper, sampler)
# attractors = extract_attractors(mapper)

# ## found two attractors: one is a limit cycle, the other is chaotic
# ## visualize them
# plot_attractors(attractors)

# ## continue all attractors and their basin fractions across any arbigrary
# ## curve in parameter space using a global continuation algorithm
# algo = AttractorSeedContinueMatch(mapper)
# params(θ) = [1 => 5 + 0.5cos(θ), 2 => 0.1 + 0.01sin(θ)]
# angles = range(0, 2π; length = 101)
# pcurve = params.(angles)
# fractions_cont, attractors_cont = global_continuation(
# 	algo, pcurve, sampler; samples_per_parameter = 1_000
# )

# ## and visualize the results
# fig = plot_basins_attractors_curves(
# 	fractions_cont, attractors_cont, A -> minimum(A[:, 1]), angles; add_legend = false
# )
# ```


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
using OrdinaryDiffEqVerner # for accessing advanced ODE Solvers

function modified_lorenz_rule(u, p, t)
    x, y, z = u; a, b = p
    dx = y - x
    dy = - x*z + b*abs(z)
    dz = x*y - a
    return SVector(dx, dy, dz)
end

p0 = [5.0, 0.1] # parameters
u0 = [-4.0, 5, 0] # state
diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9, dt = 0.01) # solver options
ds = CoupledODEs(modified_lorenz_rule, u0, p0; diffeq)

# ## Finding attractors

# In this tutorial we will utilize two methods for finding attractors in dynamical systems.
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
    range(-10.0, 10.0; length = 150), # x
    range(-15.0, 15.0; length = 150), # y
    range(-15.0, 15.0; length = 150), # z
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
# that they converged to a different attractor. Indeed,

mapper([1.0, 3, 2])

# gets the same ID as the first initial condition.

# This functionality is already incredibly powerful!
# To our knowledge the DynamicalSystems.jl library is the only
# dynamical systems software (in any language) that provides such an
# infrastructure for mapping initial conditions of any arbitrary dynamical system to its
# unique attractors. And this is only the tip of this iceberg!
# The rest of the functionality of Attractors.jl is all full of brand new cutting edge
# progress in dynamical systems research.

# Okay, back to the tutorial now!
# The found attractors are stored in the mapper internally, to obtain them we
# use the function

attractors = extract_attractors(mapper)

# In Attractors.jl, all information regarding attractors is always a standard Julia
# `Dict`, which maps attractor IDs (positive integers) to the corresponding quantity.
# Here the quantity are the attractors themselves, represented as `StateSpaceSet`.

# We can visualize them with the convenience plotting function
using CairoMakie
plot_attractors(attractors)

# (this convenience function is a simple loop over scattering the values of
# the `attractors` dictionary)

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

# ## Different Attractor Mapper

# Attractors.jl utilizes composable interfaces throughout its functionality.
# In the above example we used one particular method to find attractors,
# via recurrences in the state space. An alternative is [`AttractorsViaFeaturizing`](@ref).

# For this method, we need to provide a "featurizing" function that given an
# trajectory (which is likely an attractor), it returns some features that will
# hopefully distinguish different attractors in a subsequent grouping step.
# Finding good features is typically a trial-and-error process, but for our system
# we already have some good features:

using Statistics: mean

function featurizer(A, t) # t is the time vector associated with trajectory A
    xmin = minimum(A[:, 1])
    ycen = mean(A[:, 2])
    return SVector(xmin, ycen)
end

# from which we initialize

mapper2 = AttractorsViaFeaturizing(ds, featurizer; Δt = 0.1)

# [`AttractorsViaFeaturizing`](@ref) allows for a third input, which is a
# "grouping configuration", that dictates how features will be grouped into
# attractors, as features are extracted from (randomly) sampled state space trajectories.
# In this tutorial we leave it at its default value, which is clustering using the DBSCAN
# algorithm. The keyword arguments are meta parameters which control how long
# to integrate each initial condition for, and what sampling time, to produce
# a trajectory `A` given to the `featurizer` function. Because one of the two attractors
# is chaotic, we need denser sampling time than the default.

# We can use `mapper2` exactly as `mapper`:

fs2 = basins_fractions(mapper2, sampler)

attractors2 = extract_attractors(mapper2)

plot_attractors(attractors2)

# This mapper also found the attractors, but we should warn you: this mapper is less
# robust than [`AttractorsViaRecurrences`](@ref). One of the reasons for this is
# that [`AttractorsViaFeaturizing`](@ref) is not auto-terminating. For example, if we do not
# have enough transient integration time, the two attractors will get confused into one:

mapper3 = AttractorsViaFeaturizing(ds, featurizer; Ttr = 10, Δt = 0.1)
fs3 = basins_fractions(mapper3, sampler)
attractors3 = extract_attractors(mapper3)
plot_attractors(attractors3)

# On the other hand, the downside of [`AttractorsViaRecurrences`](@ref) is that
# it can take quite a while to converge for chaotic high dimensional systems.

# ## [Global continuation](@id global_cont_tutorial)

# If you have heard before the word "continuation", then you are likely aware of the
# **traditional continuation-based bifurcation analysis (CBA)** offered by many software,
# such as AUTO, MatCont, and in Julia [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).
# These software perform **local continuation**.
# Here we offer a completely different kind of continuation called **global continuation**.

# The local continuation continues the curves of individual _fixed
# points (and under some conditions limit cycles)_ across the joint state-parameter space and
# tracks their _local (linear) stability_.
# This approach needs to manually be "re-run" for every individual branch of fixed points
# or limit cycles.
# The global continuation in Attractors.jl finds _all_ attractors, _including chaotic
# or quasiperiodic ones_,
# in the whole of the state space (that it searches in), without manual intervention.
# It then continues all of these attractors concurrently along a parameter axis.
# Additionally, the global continuation tracks a _nonlocal_ stability property which by
# default is the basin fraction.

# Because all attractors are simultaneously
# tracked across the parameter axis, the user may arbitrarily estimate _any_
# property of the attractors and how it varies as the parameter varies.
# A more detailed comparison between these two approaches can be found in [Datseris2023](@cite).
# See also the [comparison page](@ref bfkit_comparison) in our docs
# that attempts to do the same analysis of our Tutorial with traditional continuation software.

# To perform a global continuation is surprisingly simple. First, we decide what parameter,
# and what range of that parameter, to continue over:

prange = 4.5:0.01:6
pidx = 1 # index of the parameter

# Then, we may call the [`global_continuation`](@ref) function.
# We have to provide a continuation algorithm, which itself references an [`AttractorMapper`](@ref).
# In this example we will re-use the `mapper` to create the "flagship product" of Attractors.jl
# which is the general [`AttractorSeedContinueMatch`](@ref).
# This algorithm uses the `mapper` to find all attractors at each parameter value
# and from the found attractors it continues them along a parameter axis
# using a seeding process (see its documentation string).
# Then, it performs a "matching" step, ensuring a "continuity" of the attractor
# label across the parameter axis. For now we ignore the matching step, leaving it to the
# default value. We'll use the `mapper` we created above and define

ascm = AttractorSeedContinueMatch(mapper)

# and call

fractions_cont, attractors_cont = global_continuation(
	ascm, prange, pidx, sampler; samples_per_parameter = 1_000
)

# the output is given as two vectors. Each vector is a dictionary
# mapping attractor IDs to their basin fractions, or their state space sets, respectively.
# Both vectors have the same size as the parameter range.
# For example, the attractors at the 34-th parameter value are:

attractors_cont[34]

# If you want to transform the output to the alternative format of
# a dictionary of vectors, use [`continuation_series`](@ref).
# You typically don't have to though, because there is a fantastic convenience function
# for animating the attractors evolution, that utilizes things we have
# already defined:

animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, prange, pidx;
);

# ```@raw html
# <video width="auto" controls loop>
# <source src="../attracont.mp4" type="video/mp4">
# </video>
# ```

# Hah, how cool is that! The attractors pop in and out of existence like out of nowhere!
# It can be difficult to find these attractors in traditional continuation software
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

# ## Different matching procedures

# By default attractors are matched by their distance in state space.
# The default matcher is [`MatchBySSSetDistance`](@ref), and is given implicitly
# as a default 2nd argument when creating [`AttractorSeedContinueMatch`](@ref).
# But like anything else in Attractors.jl, "matchers" also follow a well-defined
# and extendable interface, see [`IDMatchers`](@ref) for that.

# Let's say that the default matching that we chose above isn't desirable.
# For example, one may argue that the attractor that pops up
# at the end of the continuation should have been assigned the same ID
# as attractor 1, because they are both to the left (see the video above).
# In reality one wouldn't really request that, because looking
# the video of attractors above shows that the attractors labelled "1", "2", and "3"
# are all completely different. But we argue here for example that "3" should have been
# the same as "1".

# Thankfully during a global continuation the "matching" step is completely
# separated from the "finding and continuing" step. If we don't like the
# initial matching, we can call [`match_sequentially!`](@ref) with a new
# instance of a matcher, and match again, without having to recompute
# the attractors and their basin fractions.
# For example, using this matcher:

matcher = MatchBySSSetDistance(use_vanished = true)

# will compare a new attractor with the latest instance of attractors
# with a given ID that have ever existed, irrespectively if they exist in the
# current parameter or not. This means, that the attractor "3" would in fact be compared
# with both attractor "2" and "1", even if "1" doesn't exist in the parameter "3"
# started existing at. And because "3" is closer to "1" than to "2", it will get
# matched to attractor "1" and get the same ID.

# Let's see this in action:

attractors_cont2 = deepcopy(attractors_cont)

match_sequentially!(attractors_cont2, matcher)

fig = plot_attractors_curves(
	attractors_cont2, A -> minimum(A[:, 1]), prange,
)

# and as we can see, the new attractor at the end of the parameter range got
# assigned the same ID as the original attractor "1".
# For more ways of matching attractors see [`IDMatcher`](@ref).

# ## Continuation along arbitrary parameter curves

# One of the many advantages of the global continuation is that we can choose
# what parameters to continue over. We can provide any arbitrary curve
# in parameter space. This is possible because (1) finding and matching attractors
# are two completely orthogonal steps, and (2) it is completely fine for
# attractors to dissapear (and perhaps re-appear) during a global continuation.

#For example, we can probe an elipsoid defined as

params(θ) = [1 => 5 + 0.5cos(θ), 2 => 0.1 + 0.01sin(θ)]
θs = range(0, 2π; length = 101)
pcurve = params.(θs)

# here each component maps the parameter index to its value.
# We can just give this `pcurve` to the global continuation,
# using the same mapper and continuation algorithm,
# but adjusting the matching process so that vanished attractors
# are kept in "memory"

matcher = MatchBySSSetDistance(use_vanished = true)

ascm = AttractorSeedContinueMatch(mapper, matcher)

fractions_cont, attractors_cont = global_continuation(
	ascm, pcurve, sampler; samples_per_parameter = 1_000
)

# and animate the result
animate_attractors_continuation(
    ds, attractors_cont, fractions_cont, pcurve;
    savename = "curvecont.mp4"
);

# ```@raw html
# <video width="auto" controls loop>
# <source src="../curvecont.mp4" type="video/mp4">
# </video>
# ```

# %% #src
# ## Enhancing the continuation with more stability measures or other quantities

# The standard stability measure that is reported during a global continuation
# is the basin fractions. This is primarily because it is computed automatically
# as we find the different attractors. There are many more stability measures
# that could be more useful in different contexts. Attractors.jl offers
# the unique possibility of estimating _almost all_ known measures of stability in the
# literature of dynamical systems during a _single_ global continuation pass.
# This is done with the [`StabilityMeasuresAccumulator`](@ref) data structure.
# You can visit its documentation string to learn about all different stability measures.
# If you find some stability measures not included in the [`StabilityMeasuresAccumulator`](@ref)
# then either open an Issue and tell us about it or even better make a Pull Request and
# contribute it yourself!

# Using [`StabilityMeasuresAccumulator`](@ref) is very easy. If we have already performed
# a global continuation then we can utilize the function [`stability_measures_along_continuation`](@ref)
# to run through it again and estimate now all stability measures.

result = stability_measures_along_continuation(
    ds, attractors_cont, pcurve, sampler; ε = 0.1
)
keys(result)

# The result is a dictionary mapping the stability measure name (as a string) to
# the continuation of the measure.
# We can see there are quite a lot of measures that have been estimated!
# So the values of `result` are the same type as
# the `fractions_cont` we have computed before. This means it is straightforward to
# visualize these new stability measures. For example, let's say we want to visualize

chosen = ["median_convergence_pace", "minimal_critical_shock_magnitude"]
abbrev = ["MCP", "MCS"]

# then we just use the visualization function `plot_continuation_curves`

ukeys = unique_keys(attractors_cont) # so that we ignore key -1 (divergent orbits)
fig = plot_attractors_curves(attractors_cont, A -> minimum(A[:, 1]), θs)
for (i, c) in enumerate(chosen)
    ax = Axis(fig[1 - i, 1]; ylabel = abbrev[i])
    measure_cont = result[c]
    plot_continuation_curves!(ax, measure_cont, θs; ukeys, add_legend = false)
    hidexdecorations!(ax; grid = false)
end
resize!(fig, 600, 500)
fig


# For more specialization on estimating these stability measures,
# see the documentation of [`StabilityMeasuresAccumulator`](@ref).


# One of the biggest strengths of Attractors.jl is that it is not an isolated software.
# It is part of **DynamicalSystems.jl**. We can straightforwardly use any other
# functionality of the library to enhance this continuation, even beyond stability measures.
# Let's also estimate and visualize the maximum Lyapunov exponent for each attractor.
# We first import the function that estimates the MLE

using ChaosTools: lyapunov

# and then, the estimation itself is rather simple:

lis = map(enumerate(pcurve)) do (i, p) # loop over parameters
    set_parameters!(ds, p) # important! We use the dynamical system!
    attractors = attractors_cont[i]
    ## Return a dictionary mapping attractor IDs to their MLE
    Dict(k => lyapunov(ds, 10000.0; u0 = A[1]) for (k, A) in attractors)
end

# The above `map` loop may be intimidating if you are a beginner, but it is
# really just a shorter way to write a `for` loop for our example.
# We iterate over all parameters, and for each we first update the dynamical
# system with the correct parameter, and then extract the MLE
# for each attractor. `map` just means that we don't have to pre-allocate a
# new vector before the loop; it creates it for us.

# We now visualize the MLE with the same way as any other quantity over the continuation:

axλ = Axis(fig[-length(chosen), 1]; ylabel = "MLE")
hidexdecorations!(axλ)
plot_continuation_curves!(axλ, lis, θs; add_legend = false)
fig

# This reveals crucial information for tha attractors, whether they are chaotic or not,
# that we would otherwise obtain only by visualizing the system dynamics at every single
# point in the continuation parameter.

# ## Conclusion and comparison with traditional local continuation

# We've reached the end of the tutorial! Some aspects we haven't highlighted is
# how most of the infrastructure of Attractors.jl is fully extendable.
# You will see this when reading the documentation strings of key structures
# like [`AttractorMapper`](@ref). All documentation strings are in the [API](@ref) page.
# See the [examples](@ref examples) page for more varied applications.
# And lastly, see the [comparison page](@ref bfkit_comparison) in our docs
# that attempts to do the same analysis of our Tutorial with traditional local continuation and bifurcation analysis software
# showing that (at least for this example) using Attractors.jl is clearly beneficial
# over the alternatives.