using Attractors

# Coupled Hénon maps
function F!(du, u ,p, n)
    x,y,u,v = u
    A = 3; B = 0.3; C = 5.; D = 0.3; k = 0.4; 
    du[1] = A - x^2 + B*y + k*(x-u)
    du[2] = x
    du[3] = C - u^2 + D*v + k*(u-x)
    du[4] = u
    return 
end

# The region should not contain any attractors. 
R_min = [-4; -4.; -4.; -4.] 
R_max = [4.; 4.; 4.; 4.]


# Initial box and initial condition.
sampler, isinside = statespace_sampler(HRectangle(R_min,R_max))
x0 = sampler()
df = DeterministicIteratedMap(F!, x0) 
xi = stagger_trajectory!(df, x0, 20, isinside; δ₀= 2.) 
if isnothing(xi) 
    error("couldn't find a suitable point")
end
@show Tp = Attractors.escape_time!(df, xi, isinside)

v = stagger_and_step!(df, x0, 10000, isinside; stagger_mode = :adaptive, δ = 1e-4, Tm = 10, max_steps = Int(1e5), δ₀ = 2.) 
v = hcat(v...)'

# using CairoMakie
# scatter(v[:,1], v[:,3]; markersize = 3)

