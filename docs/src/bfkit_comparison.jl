
# ## Comparison with traditional continuation software

# As we discussed in subsection on [global continuation](@ref global_cont_tutorial),
# the approach of Attractors.jl is fundamentally different from traditional continuation
# software. In this section we use BifurcationKit.jl to highlight the differences,
# and attempt to find and continue the limit cycle of our example
# (to our knowledge there exists no software on dynamical systems
# beyond Attractors.jl can find _chaotic_ attractors).

# To use BifurcationKit.jl (BK) for periodic orbits (POs) we need to choose one of its
# several Newton-based algorithms for POs, and in addition supply it with
# both an initial guess for the location of the periodic orbit, as well as
# a guess for the period. Finding a periodic orbit this way is already considered
# an advanced use case in BK documentation,
# requiring "high level of knowledge of (numerical) bifurcation theory"
# (following the example of the [Periodic predator prey model](https://bifurcationkit.github.io/BifurcationKitDocs.jl/v0.3/tutorials/ode/tutorialsCodim2PO/#Periodic-predator-prey-model)).
# For Attractors.jl on the other hand, this is as basic of a use-case as it can get,
# which highlights the simplicity of our computational approach.
# Notice also that, at least from the public interface,
# it is not possible to use BK outside of a continuation context, to e.g.,
# just find a limit cycle at a particular parameter. We have to do the full continuation
# or limit it to an infinitesimally small parameter range.

# To use BK we need to import it and initialize
# various continuation-related structures.
# The key structure is the `ShootingProblem`, which is a "problem type"
# that utilizes a multiple shooting method to find periodic orbits
# (see BK docs). The entire input BK requires to find a periodic orbit is:
# 1. a `BK.ShootingProblem` shooting problem
# 1. a `BK.BifurcationProblem`
# 1. a `DifferentialEquations.Solution`
# 1. an estimate of the period
# 1. a `BK.ContinuationPar` parameter container
# 1. a predictor for the continuation
# 1. arguments for what aspect of the periodic orbit to record.
# That's already lots of input, and each one of this can have dozens of options.

# Let's start with the bifurcation problem.
# This is basically the same thing as a `DynamicalSystem`, but BK does not
# support efficient `StaticVector`-based out of place format for low dimensional
# systems (see main tutorial of DynamicalSystems.jl if you don't understand what this means).
# So we have to re-create

# %% #src
import BifurcationKit as BK
using OrdinaryDiffEq

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

bf_prob = BK.BifurcationProblem(
    modified_lorenz_rule!, u0, p0, (BK.@lens _[pidx])
)

# Then the shooting problem
shooting_prob = BK.ShootingProblem(M=5)

# and then a full solution structure from DifferentialEquations.jl, that
# **must** start exactly on the periodic orbit. Requiring that the solution
# starts on the periodic orbit defeats the purpose of "wanting to find it",
# but oh well, we do as we must.

## point_on_lc = Vector(attractors_cont[1][1][end]) # we copy paste the guess below
point_on_lc = [
    -1.622956992666447,
    -4.527917984019188,
    -5.178825669659272,
]

ode_prob = ODEProblem(modified_lorenz_rule!, point_on_lc, p0, (0.0, 40.0))
sol = OrdinaryDiffEq.solve(ode_prob; alg = Tsit5())
sol[end]

# We need an estimate of the period besides providing the full DifferentialEquations.jl
# solution. From the figure this appears to be around 20.0 (note: the periodic
# orbit wraps around four times before repeating exactly).

# Right, and lastly we need a continuation parameter container,
# which has some options regarding convergence and stability that one
# would need to fine tune to the problem at hand.
opts_br = BK.ContinuationPar(
    p_min = prange[1], p_max = prange[1],
    ds = 0.002, dsmax = 0.01, n_inversion = 6,
    detect_bifurcation = 3, max_bisection_steps = 25, nev = 4,
    max_steps = 2000, tol_stability = 1e-3
)

# Now we put everything together in
probsh, cish = BK.generate_ci_problem(
    shooting_prob, bf_prob, ode_prob, sol, 20.0; alg = Vern9()
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

try
    br_fold_sh = BK.continuation(probsh, cish, predictor, opts_br;
        verbosity = 3, plot = false,
        argspo...
    )
catch err
    return string(err)
end

# Unfortunately, the above errors. So we would need to go back to the drawing board,
# and fine tune literally every single input to the final call.

# **But this is about the point we give up**. Attractors.jl found the limit cycle
# without any problem, and with much more robust algorithms (with respect to input parameters).
# Given a dynamical system and parameter range (input to both Attractors.jl global
# continuation or traditional bifurcation-based continuation), the Attractors.jl
# full code is:

ds = CoupledODEs(modified_lorenz_rule!, u0, p0)

grid = (
    range(-10.0, 10.0; length = 100), # x
    range(-15.0, 15.0; length = 100), # y
    range(-15.0, 15.0; length = 100), # z
)

mapper = AttractorsViaRecurrences(ds, grid;
    consecutive_recurrences = 1000,
    attractor_locate_steps = 1000,
    consecutive_lost_steps = 100,
)

sampler, = statespace_sampler(grid)

matcher = MatchBySSSetDistance()

algo = AttractorSeedContinueMatch(mapper, matcher)

fractions_cont, attractors_cont = global_continuation(
	algo, prange, pidx, sampler; samples_per_parameter = 1_000
)

# Of course, the above code didn't find and continue just a single limit cycle.
# It found all system attractors, and it didn't require a specific initial condition
# as a guess, but rather a box that may contain attractors.
# The code is relatively agnostic to the dynamical system
# as well, meaning, the algorithms are robust: the recurrences mapper would work for
# most systems with this parameters; it doesn't need as much tuning
# (see [Datseris2022](@cite) for a proof of that).
# Not only only does the code find all attractors, it also finds chaotic attractors.
# And furthermore, it estimates a more powerful nonlocal measure of stability.