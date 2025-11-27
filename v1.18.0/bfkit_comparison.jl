# # [Comparison with traditional continuation & bifurcation software](@id bfkit_comparison)

# !!! note "Continues from tutorial"
#     This page continues after the end of the main [Tutorial](@ref)!
#     Ensure you have gon through it first!

# As we discussed in the subsection on [global continuation](@ref global_cont_tutorial),
# the approach of Attractors.jl is fundamentally different from traditional continuation
# and bifurcation software like AUTO, MatCont, BifurcationKit.jl.
# Nevertheless in this page we will compare using BifurcationKit.jl
# to attempt to find and continue the limit cycle of the tutorial modified Lorenz-like system.
# We forfeit looking for the chaotic attractors, as to our knowledge there exists no
# software on dynamical systems# beyond Attractors.jl can find chaotic attractors.

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
# The key structure is the `ShootingProblem`, which is a "problem type"
# that utilizes a multiple shooting method to find periodic orbits
# (see BK docs). The entire input BK requires to find a periodic orbit is:
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

bf_prob = BK.BifurcationProblem(
    modified_lorenz_rule!, u0, p0, (BK.@lens _[pidx])
)

# Then the algorithm particularly for periodic orbits
periodic_orbit_algo = BK.ShootingProblem(M = 5)

# and then a full solution structure from DifferentialEquations.jl, that
# **must** start on the periodic orbit. Requiring that the solution
# starts on the periodic orbit defeats the purpose of "wanting to find it",
# but oh well, we do as we must.

## point_on_lc = Vector(attractors_cont[1][1][end]) # we copy paste the guess below
point_on_lc = [
    -1.622956992666447,
    -4.527917984019188,
    -5.178825669659272,
]

ode_prob = ODEProblem(modified_lorenz_rule!, point_on_lc, (0.0, 100.0), p0)
sol = OrdinaryDiffEq.solve(ode_prob; alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
lines(sol.t, sol[1, :])

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
    max_steps = 2000, tol_stability = 1e-3
)

# Now we put everything together in

probsh, cish = BK.generate_ci_problem(
    periodic_orbit_algo, bf_prob, ode_prob, sol, 20.0; alg = Vern9(), maxiters = Int(1e6),
)

# To call the continuation we need to also tell it what aspects of the
# periodic orbit to record, so we define

argspo = (record_from_solution = (x, p) -> begin
		xtt = get_periodic_orbit(p.prob, x, p.p)
		return (max = maximum(xtt[1,:]),
				min = minimum(xtt[1,:]),
				period = getperiod(p.prob, x, p.p))
	end,
)

# we also define the predictor
predictor = BK.PALC(tangent = BK.Bordered())

# and _finally_ call the continuation from BK

br_fold_sh = BK.continuation(probsh, cish, predictor, opts_br;
    verbosity = 3, plot = false,
    argspo...
)

# The above code takes a good 30 seconds to run, but unfortunately it fails.
# During the evaluation we get lots of warnings of the type
# ```julia-repl
# ┌ Warning: Interrupted. Larger maxiters is needed.
# If you are using an integrator for non-stiff ODEs or an automatic switching
# algorithm (the default), you may want to consider using a method for stiff equations.
# See the solver pages for more details (e.g. https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Stiff-Problems).
# ```
# which actually, have nothing to do with our solver not being stiff. What happens here
# is that the Newton algorithm used to find periodic orbits never converges.
# This is almost surely because these Newton-based periodic orbit algorithms
# fail to converge in the presence of chaos, and the periodic orbit we try to find
# here starts without chaos, but continues to exist while chaos exists,
# hence the continuation fails.

# Since BifurcationKit.jl doesn't offer a (public) interface to find
# "just the limit cycle at a single parameter" like Attractors.jl does, we can't
# use it to detect when chaos starts by utilizing a `for` loop and stopping at the first
# parameter where the Newton algorithm fails to converge. We have to manually
# scan the parameter axis and evolve various initial conditions with random sampling
# and then "look" at the trajectories to see if they are chaotic or not
# (or, use a tool from DynamicalSystems.jl like `lyapunov`). Then, we could
# limit the continuation at the parameter that chaos starts.
# Alternatively one could try to fine-tune the parameters in `BK.ContinuationPar`
# or in `BK.ShootingProblem`, but this is unlikely to work.

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

# This code takes about 15 seconds to run.
# Of course, the above code didn't find and continue just a single limit cycle.
# It found all system attractors, and it didn't require a specific initial condition
# as a guess, or a period, but rather an arbitrarily large box that may contain attractors.
# Not only only that, but some of the found attractors are _chaotic_!

# The code is relatively agnostic to the dynamical system
# as well, meaning, the algorithms are robust: the recurrences mapper would work for
# most systems with these parameters; it doesn't need as much tuning
# (see [Datseris2022](@cite) for a proof of that).
# And furthermore, it estimates a more general nonlocal measure of stability,
# in the sense that if a set is nonlocally stable, it is guaranteed to be locally stable,
# however the other way around isn't guaranteed.

# Now, this example is one that definitely we Attractors.jl is much more suitable for.
# But to be transparent and discuss the silver linings,
# BK can do some things not possible (yet) with Attractors.jl, such find
# the unstable branches of fixed points and limit cycles (in simpler example systems than the one here).
# Now, whether the unstable branches are useful or not, depends on the research question.
# Beyond this however, BK is also optimised for PDE systems, while Attractors.jl isn't.
# A more thorough comparison of the two approaches (that was based on and older and
# less powerful version of Attractors.jl) is discussed in [Datseris2023](@cite).