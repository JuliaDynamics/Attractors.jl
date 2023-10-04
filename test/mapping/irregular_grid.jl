using Attractors
using Test

function newton_map(z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
yg = collect(range(-1.5, 1.5; length = 400))
xg = log.(collect(range(exp(-1.5), exp(1.5), length=400)))


newton = AttractorsViaRecurrences(ds, (xg, yg);
    sparse = false, mx_chk_lost = 1000, Dt =1,
)


@test ((newton([-0.5, 0.86]) != newton([-0.5, -0.86]))& (newton([-0.5, 0.86]) != newton([1.0, 0.0])) & (newton([-0.5, -0.86]) != newton([1.0, 0.0])))

basins, attractors = basins_of_attraction(newton)

@test length(attractors) == 3


xg = vcat(range(-1.5, 0, length = 300), range(0, 1.5, length = 100))
yg = vcat(range(-1.5, 0, length = 100), range(0, 1.5, length = 300))
newton = AttractorsViaRecurrences(ds, (xg, yg);
    sparse = false, mx_chk_lost = 1000, Dt = 1,
)
@test ((newton([-0.5, 0.86]) != newton([-0.5, -0.86]))& (newton([-0.5, 0.86]) != newton([1.0, 0.0])) & (newton([-0.5, -0.86]) != newton([1.0, 0.0])))
basins, attractors = basins_of_attraction(newton)

@test length(attractors) == 3

xg = yg = range(-1.5, 1.5, length = 400)
newton = AttractorsViaRecurrences(ds, (xg, yg);
    sparse = false, mx_chk_lost = 1000, Dt =1,
)
@test ((newton([-0.5, 0.86]) != newton([-0.5, -0.86]))& (newton([-0.5, 0.86]) != newton([1.0, 0.0])) & (newton([-0.5, -0.86]) != newton([1.0, 0.0])))
basins, attractors = basins_of_attraction(newton)

@test length(attractors) == 3

#########################################################################
###                 Basin cell index test                             ###
#########################################################################



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


xg = yg = range(0, 5, length = 6)

grid = subdivision_based_grid(ds, (xg, yg); maxlevel = 3)
lvl_array = grid.lvl_array

@test (Attractors.basin_cell_index((2, 1.1) , grid) == CartesianIndex((17,9)))
@test (Attractors.basin_cell_index((3.6, 2.1) , grid) == CartesianIndex((29,17)))
@test (Attractors.basin_cell_index((4.9, 4.9) , grid) == CartesianIndex((33,33)))
@test (Attractors.basin_cell_index((0.5, 0.5) , grid) == CartesianIndex((5,5)))
@test (Attractors.basin_cell_index((3.4, 1.6) , grid) == CartesianIndex((27,13)))
@test (Attractors.basin_cell_index((4.7, 2.2) , grid) == CartesianIndex((37,17)))




############################################################
####  SubdivisionBasedGrid tests                        ####
############################################################
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



#####################
## IrregularGrid  ###
#####################
attractors = []
for pow in (1, 2)
    xg = yg = range(0, 18.0^(1/pow); length = 200).^pow
    mapper = AttractorsViaRecurrences(ds, (xg, yg);
        Dt = 0.1, sparse = true,
        mx_chk_fnd_att = 10, mx_chk_loc_att = 10,
        mx_chk_safety = 1000,
    )
    sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
    fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
    push!(attractors, extract_attractors(mapper))

end

@test length(vec(values(attractors)[1][1])) *10 < length(vec(values(attractors)[2][1]))



######################
### New approach  ####
######################


xg = yg = range(0, 18, length = 30)

grid = subdivision_based_grid(ds, (xg, yg))

mapper = AttractorsViaRecurrences(ds, grid;
        Dt = 0.1, sparse = true,
        mx_chk_fnd_att = 10, mx_chk_loc_att = 10,
        mx_chk_safety = 1000,
    )

sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
attractors_SBD = extract_attractors(mapper)

###############################
## same setup, regular grid  ##
###############################

xg = yg = range(0, 18, length = 30)
mapper = AttractorsViaRecurrences(ds, (xg, yg);
        Dt = 0.1, sparse = true,
        mx_chk_fnd_att = 10, mx_chk_loc_att = 10,
        mx_chk_safety = 1000,
    )


sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
attractors_reg = extract_attractors(mapper)


@test (length(vec(attractors_SBD[1]))/10) > length(vec(attractors_reg[1]))


############################################################
####  Automatic Δt works                                ####
############################################################
reinit!(ds)

xg = yg = range(0, 18, length = 30)
grid0 = (xg, yg)
xg = yg = range(0, 18.0^(1/2); length = 20).^2
grid1 = (xg, yg)
grid2 = SubdivisionBasedGrid(grid0, rand([0, 1, 2], (30, 30)))

using Attractors: RegularGrid, IrregularGrid

Dt0 = automatic_Δt_basins(ds, RegularGrid(grid0))
Dt1 = automatic_Δt_basins(ds, IrregularGrid(grid1))
Dt2 = automatic_Δt_basins(ds, grid2)

@test Dt0 > 0
@test Dt1 > 0
@test Dt2 > 0