using Attractors
using Test

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
xi = stagger_trajectory!(df, x0, 30, isinside; δ₀ = 2.) 
@test isinside(xi)
@test Attractors.escape_time!(df, xi, isinside) > 30

# Test if all the points have escape time ≥ Tm 
# :exp mode
v = stagger_and_step!(df, xi, 10, isinside; δ₀ = 1e-3, stagger_mode = :exp) 
for u in v
    @test Attractors.escape_time!(df, u, isinside) ≥ 30
end

# Test if all the points have escape time ≥ Tm 
# :unif mode
v = stagger_and_step!(df, xi, 10, isinside; δ₀ = 1e-3, stagger_mode = :unif) 
for u in v
    @test Attractors.escape_time!(df, u, isinside) ≥ 30
end

# Test if all the points have escape time ≥ Tm 
# :adaptive mode
v = stagger_and_step!(df, xi, 10, isinside; δ₀ = 1e-3, stagger_mode = :adaptive) 
for u in v
    @test Attractors.escape_time!(df, u, isinside) ≥ 30 - 1
end
end


