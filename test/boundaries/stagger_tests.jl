using Attractors
using Test
using LinearAlgebra:norm
using Random

@testset "Stagger Tests" begin

# Dynamical system with a saddle at 0
function F!(du, u ,p, n)
    x,y = u
    du[1] = x + y 
    du[2] = x - y
    return 
end


R_min = [-1.; -1.]; R_max = [1.; 1.]

# Initial condition.
sampler, isinside = statespace_sampler(HRectangle(R_min,R_max))
x0 = sampler()
df = DeterministicIteratedMap(F!, x0) 

xi = Attractors.stagger_trajectory(df, x0, 30, isinside, :exp, 2., 1.1, Int(1e5), 10000, Xoshiro(123)) 

@test isinside(xi)
@test Attractors.escape_time!(df, xi, isinside, 10000) > 30

# Test if all the points have escape time ≥ Tm 

for m in [:exp, :unif, :adaptive]
    v = stagger_and_step(df, xi, 10, isinside; δ₀ = 1e-3, stagger_mode = m) 
    for u in v
        @test Attractors.escape_time!(df, u, isinside, 10000) ≥ 30 - 1
    end
end


# Dynamical system with a stable fixed point at 0
function F_stable!(du, u ,p, n)
    x,y = u
    du[1] = -x 
    du[2] = -y 
    return 
end

df = DeterministicIteratedMap(F_stable!, x0) 
# Test if the program throws an error when the 
# trajectory stays inside the square
@test_throws ErrorException Attractors.escape_time!(df, sampler(), isinside, 10000)

# Test with a real continuous dynamical system. 
# There is a saddle point at (0,0,0)
function thomas_rule(u, p, t)
    x,y,z = u
    b = p
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector{3}(xdot, ydot, zdot)
end
ds =  CoupledODEs(thomas_rule, rand(3), 0.1665)
R_min = -1*ones(3); R_max = 1*ones(3)
sampler, isinside = statespace_sampler(HRectangle(R_min,R_max))
xi = Attractors.stagger_trajectory(ds, sampler(), 30, isinside, :exp, 2., 1.1, Int(1e5), 10000, Xoshiro(123)) 
 
@test isinside(xi)
@test Attractors.escape_time!(ds, xi, isinside, 10000) > 30

v = stagger_and_step(ds, xi, 1000, isinside; stagger_mode = :adaptive, δ = 1e-4, Tm = 10, max_steps = Int(1e5), δ₀ = 1., rng = Xoshiro(123)) 
# The last point of the series should be close to the saddle 
@test norm(v[end]) < 0.1 

end


