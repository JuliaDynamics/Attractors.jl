# # [Comparison with traditional local continuation and bifurcation analysis software](@id bfkit_comparison)

# !!! note "Continues from tutorial"
#     This page continues after the end of the main [Tutorial](@ref)!
#     Ensure you have gon through it first!

# As we discussed in the subsection on [global continuation](@ref global_cont_tutorial),
# the approach of Attractors.jl is fundamentally different from traditional local continuation
# and bifurcation analysis software like AUTO, MatCont, or BifurcationKit.jl.
# Nevertheless in this page we will compare using BifurcationKit.jl
# to attempt to find and continue the limit cycle of the tutorial modified Lorenz-like system.
# We forfeit looking for the chaotic attractors, as to our knowledge there exists no
# software on dynamical systems beyond Attractors.jl can find chaotic attractors.

# The goal of this comparison is to highlight
# the differences in usage complexity and overall capability
# when using Attractors.jl or traditional continuation tools to study
# **complex** dynamical systems.

# ## BifurcationKit.jl version

# To use BifurcationKit.jl (BK) for periodic orbits (POs) we need to choose one of its
# several Newton-based algorithms for POs, and in addition supply it with
# both an initial guess for the location of the periodic orbit, as well as
# a guess for the period.
# In this example we translate almost verbatim the example of the [Periodic predator prey model](https://bifurcationkit.github.io/BifurcationKitDocs.jl/v0.3/tutorials/ode/tutorialsCodim2PO/#Periodic-predator-prey-model)).
# from the BK docs. Finding a periodic orbit this way is already considered
# an advanced use case in BK documentation,
# requiring "high level of knowledge of (numerical) bifurcation theory".
# For Attractors.jl on the other hand, this is as basic of a use-case as it can get,
# which highlights the simplicity of our computational approach.

# To use BK we need to import it and initialize
# various continuation-related structures.
# The entire input BK requires to find a periodic orbit is:
# 1. a periodic orbit problem like `BK.ShootingProblem` or `BK.PeriodicOrbitTrapProblem`
#    (and its meta parameters)
# 1. a `BK.BifurcationProblem`
# 1. a `DifferentialEquations.Solution`
# 1. an estimate of the period
# 1. a `BK.ContinuationPar` parameter container
# 1. a predictor for the continuation
# 1. arguments for what aspect of the periodic orbit to record.

# Some of this input isn't particularly important,
# but some of it is can be crucial and the values of the meta-parameters matter
# for whether the continuation will succeed or not.

# Let's start with the bifurcation problem.
# This is basically the same thing as a `DynamicalSystem`, but BK does not
# support efficient `StaticVector`-based out of place format for low dimensional
# systems (see main tutorial of DynamicalSystems.jl if you don't understand what this means).
# So we have to re-create

## BK requires to modify the rule with t=0 as well
function modified_lorenz_rule!(du, u, p, t = 0)
    x, y, z = u; a, b = p
    du[1] = y - x
    du[2] = - x*z + b*abs(z)
    du[3] = x*y - a
    return du
end
prange = 4.7:0.02:6
pidx = 1
p0 = [5.0, 0.1] # parameters
u0 = [-4.0, 5, 0] # state

# %% #src
# Now we can create the bifurcation problem

import BifurcationKit as BK
using OrdinaryDiffEq
using CairoMakie

bf_prob = BK.BifurcationProblem(
    modified_lorenz_rule!, u0, p0, (BK.@lens _[pidx])
)

# and then a full solution structure from DifferentialEquations.jl, that
# **must** start on the periodic orbit. Requiring that the solution
# starts on the periodic orbit defeats the purpose of "wanting to find it",
# but oh well, we do as we must.

## This guess comes from the Attractors.jl main tutorial
point_on_lc = [
    -1.622956992666447,
    -4.527917984019188,
    -5.178825669659272,
]

ode_prob = ODEProblem(modified_lorenz_rule!, point_on_lc, (0.0, 100.0), p0)
sol = OrdinaryDiffEq.solve(ode_prob; alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
j = length(sol)÷2
fig, ax = lines(sol.t[j:end], sol[1, j:end])
lines!(ax, sol.t[j:end], sol[2, j:end])
lines!(ax, sol.t[j:end], sol[3, j:end])
fig

# We need an estimate of the period besides providing the full DifferentialEquations.jl
# solution. From the figure this appears to be around 20.0 (note: the periodic
# orbit wraps around four times before repeating exactly).

# Right, and lastly we need a continuation parameter container,
# which has some options regarding convergence and stability that one
# would need to fine tune to the problem at hand.
opts_br = BK.ContinuationPar(
    p_min = prange[1], p_max = prange[end],
    ds = 0.002, dsmax = 0.01, n_inversion = 6,
    detect_bifurcation = 3, max_bisection_steps = 25, nev = 4,
    max_steps = 2000, tol_stability = 1e-3,
)

# We now create a periodic orbit problem type, by choosing a periodic
# orbit finding algorithm

periodic_orbit_algo = BK.PeriodicOrbitOCollProblem(40, 4)

# and creating the problem type giving the period guess 19.0

probpo, cish = BK.generate_ci_problem(
    periodic_orbit_algo, bf_prob, sol, 19.0
)

# To call the continuation we need to also tell it what aspects of the
# periodic orbit to record, so we define

argspo = (record_from_solution = (x, p) -> begin
		xtt = BK.get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = BK.getperiod(p.prob, x, p.p),
                p = p.p,)
	end,
)

# we also define the predictor
predictor = BK.PALC(tangent = BK.Bordered())

# and _finally_ call the continuation from BK

@time branch = BK.continuation(probpo, cish, predictor, opts_br;
    verbosity = 0, plot = false, argspo...
)

stability = branch.stable
color = [s ? "black" : "red" for s in stability]
marker = [s ? :circle : :x for s in stability]
scatter(branch.branch.p, branch.branch.min; color, marker)

# The above code takes about 5 seconds to run.
# Thankfully it works. Or rather, almost.
# The code finds no stable limit cycle for parameter less than 5.0,
# even though we know (see Tutorial final plot) that there is one.
# Or maybe, it is an extremely weakly chaotic attractor with MLE almost 0.
# Or maybe it is a quasiperiodic attractor. One needs to analyze further,
# but Attractors.jl finds everything without much difficulty.
# Altering the above code to start the continuation at 4.7 finds a limit cycle
# there, but only manages to continue it only up to parameter 4.74,
# instead of well into parameter 5.5 or more that Attractors.jl shows.

# ## Attractors.jl version

# We have already seen the code for this version in the main [Tutorial](@ref tutorial),
# but here we copy it again using exactly the same input as that given to BK.
# To make the comparison explicit, let's also make an enumerated list for
# the Attractors.jl info:

# 1. A `DynamicalSystem`,
# 1. an `AttractorMapper` instance (and its meta parameters).
#    For the mapper used here, `AttractorsViaRecurrences`, the meta parameters are:
#    1. A state space tesselation
#    1. A recurrence threshold
#    1. A lost iterations threshold
# 1. a global continuation algorithm, and optionally a matcher for it.
# 1. a sampler to sample initial conditions in the state space.


# %% #src
using Attractors
ds = CoupledODEs(modified_lorenz_rule!, u0, p0;
    diffeq = (alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
)

grid = (
    range(-15.0, 15.0; length = 100), # x
    range(-20.0, 20.0; length = 100), # y
    range(-20.0, 20.0; length = 100), # z
)

mapper = AttractorsViaRecurrences(ds, grid;
    consecutive_recurrences = 1000,
    consecutive_lost_steps = 100,
)

sampler, = statespace_sampler(grid)

algo = AttractorSeedContinueMatch(mapper)

fractions_cont, attractors_cont = global_continuation(
	algo, prange, pidx, sampler; samples_per_parameter = 1_000
)

plot_attractors_curves(
    attractors_cont,  A -> minimum(A[:, 1]), prange,
)

# This code takes about 15 seconds to run.
# This number however is for 1000 initial conditions, not one (i.e., the one
# branch generated during the traditional continuation).

# ## Discussion and comparison

# On the one hand, Attractors.jl found not only a single limit cycle,
# but also all system attractors, including chaotic ones.
# It didn't require any initial guess regarding the limit cycle or its period,
# but only a state space box that may contain attractors.
# Attractors.jl is extremely robust w.r.t. to its input parameters:
# the state space box can be arbitrarily large, as long as it is large enough.
# Similarly, all meta parameters of `AttractorsViaRecurrences` only need to be
# large enough; the larger, the more accurate the result.
# These algorithms are also robust in the sense of working well for
# many different types of dynamical systems, including discrete ones,
# see [Datseris2022](@cite) for a demonstration.
# And finally, Attractors.jl estimates a more general nonlocal measure of stability,
# in the sense that if a set is nonlocally stable, it is guaranteed to be locally stable,
# however the other way around isn't guaranteed.

# On the other hand, traditional local continuation can track the unstable
# branches, and automatically detect and label local bifurcations
# (note we didn't bother to plot any of the detected bifurcations - if any are found).
# In our experience having the local bifurcations is always useful.
# Now, whether the unstable branches are of a limit cycle are useful or not,
# depends on the research question and whether the analysis is done for some
# sort of real world understanding (unstable limit cycles / fixed points don't
# actually exist in the real world).
# Beyond this however, BifurcationKit.jl is also optimised for PDE systems,
# while Attractors.jl isn't.

# Now, we have to be a bit more transparent here.
# Besides this absence of stable orbits for parameter less than 5.0
# that we discussed in the end of the BifurcationKit.jl version,
# we need to say that it took **a lot of effort** to make the
# BifurcationKit.jl code work. At least, a lot of effort compared with the effort
# it took to make the Attractors.jl version work, which was simply "increase the recurrences
# threshold", which is standard practice when dealing with chaotic systems
# [Datseris2022](@cite). For example, any other of the
# periodic orbit algorithms of BifurcationKit.jl (such as shooting or trapezoid) fails.
# Using a slightly incorrect initial period guess of 20.0 instead of 19.0 also fails.
# We imagine that this sensitivity would apply also to some other of the
# several meta-parameters that enter a traditional continuation routine,
# but we didn't check further. This is exactly what we were alluding to
# in the comparison we did in [Datseris2023](@cite), that traditional
# local continuation "requires expertise and constant interventions".
# See [Datseris2023](@cite) for a more thorough comparison
# (that was based on and older and less powerful version of Attractors.jl).