
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


function newton_map(z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])

grid = (range(0, 5, length = 6),range(0, 5, length = 6))
lvl_array = subdivision_based_grid(ds, grid)
grid_steps = Dict{Int, Vector{Int}}(i => [6*2^i, 6*2^i] for i in 0:4)
grid_nfo = Attractors.SubdivisionBasedGrid(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), Array([]), grid)



# the function is modified, as current approach doesn't allow us to test it directly
# the math behind absolutely the same, just input was changed and expanded to allow quick testing
# Usage: you need to provide manually lvl of discretization of the cell with `cell_area` and
# `max_level` of discretization along the whole grid
# `cell_area` <= `max_level` in any scenraio

function basin_cell_index(y_grid_state, grid_nfo, cell_area, max_level)
    
    grid_maxima = grid_nfo.grid_maxima
    grid_minima = grid_nfo.grid_minima
    grid_steps = (grid_maxima - grid_minima .+1)./ grid_nfo.grid_steps[cell_area]
    iswithingrid = true
    @inbounds for i in eachindex(grid_minima)

        if !(grid_minima[i] ≤ y_grid_state[i] ≤ grid_maxima[i])
            iswithingrid = false
            break
        end 
    end
    if iswithingrid

        indices = @. round(Int,(y_grid_state - grid_minima)/grid_steps, RoundDown) * (2^(max_level-cell_area)) + 1

        return CartesianIndex(indices...) 

    else
        return CartesianIndex(-1)
    end
end
@testset "Basin cell index for SubdivisionBasedGrid" begin

@test basin_cell_index((1.4, 3.4), grid_nfo, 1 , 1) == CartesianIndex((3,7))
@test basin_cell_index((1.4, 3.4), grid_nfo, 1 , 3) == CartesianIndex((9,25))
@test basin_cell_index((1.4, 3.4), grid_nfo, 2 , 3) == CartesianIndex((11,27))

@test basin_cell_index((2.2, 3.7), grid_nfo, 0 , 0) == CartesianIndex((3,4))
@test basin_cell_index((2.2, 3.7), grid_nfo, 1 , 2) == CartesianIndex((9,15))
@test basin_cell_index((2.2, 3.7), grid_nfo, 3 , 3) == CartesianIndex((18,30))

@test basin_cell_index((0.35, 2.2), grid_nfo, 1 , 1) == CartesianIndex((1,5))
@test basin_cell_index((0.35, 2.2), grid_nfo, 2 , 2) == CartesianIndex((2,9))
@test basin_cell_index((0.35, 2.2), grid_nfo, 0 , 3) == CartesianIndex((1,17))

@test basin_cell_index((2.65, 3.3), grid_nfo, 0 , 1) == CartesianIndex((5,7))
@test basin_cell_index((2.65, 3.3), grid_nfo, 0 , 2) == CartesianIndex((9,13))
@test basin_cell_index((2.65, 3.3), grid_nfo, 1 , 3) == CartesianIndex((21,25))
end


############################################################
####  SubdivisionBasedGrid tests                        ####
############################################################

function newton_map(z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
xg = yg = range(-1.5, 1.5, length = 10)


grid_nfo = subdivision_based_grid(ds, (xg,yg))

newton = AttractorsViaRecurrences(ds, grid_nfo;
    sparse = false, mx_chk_lost = 1000, Dt =1,
)

@test ((newton([-0.5, 0.86]) != newton([-0.5, -0.86]))& (newton([-0.5, 0.86]) != newton([1.0, 0.0])) & (newton([-0.5, -0.86]) != newton([1.0, 0.0])))




#using CairoMakie

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

# fig = Figure()
# ax = Axis(fig[1,1])

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

    # Find attractor and its fraction (fraction is always 1 here)
    sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
    fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
    push!(attractors, extract_attractors(mapper))
    #scatter!(ax, vec(attractors[1]); markersize = 16/pow, label = "pow = $(pow)")
end

@test length(vec(values(attractors)[1][1])) *10 < length(vec(values(attractors)[2][1]))
#display(fig)



######################
### New approach  ####
######################

# fig = Figure()
# ax = Axis(fig[1,1])

xg = yg = range(0, 18, length = 30)

grid = subdivision_based_grid(ds, (xg, yg))

# constructed lvl_array
#grid.lvl_array

# passing SubdivisionBasedGrid into mapper
mapper = AttractorsViaRecurrences(ds, grid;
        Dt = 0.1, sparse = true,
        mx_chk_fnd_att = 10, mx_chk_loc_att = 10,
        mx_chk_safety = 1000,
    )

    # Find attractor and its fraction (fraction is always 1 here)
sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
attractors_SBD = extract_attractors(mapper)
#scatter!(ax, vec(attractors_SBD[1]); label = "SubdivisionBasedGrid")
#println(length(vec(attractors[1])))
#display(fig)

###############################
## same setup, regular grid  ##
###############################

xg = yg = range(0, 18, length = 30)
mapper = AttractorsViaRecurrences(ds, (xg, yg);
        Dt = 0.1, sparse = true,
        mx_chk_fnd_att = 10, mx_chk_loc_att = 10,
        mx_chk_safety = 1000,
    )

    # Find attractor and its fraction (fraction is always 1 here)
sampler, _ = statespace_sampler(HRectangle(zeros(2), fill(18.0, 2)), 42)
fractions = basins_fractions(mapper, sampler; N = 100, show_progress = false)
attractors_reg = extract_attractors(mapper)
#scatter!(ax, vec(attractors_reg[1]); label = "RegularGrid")
#println(length(vec(attractors[1])))
#axislegend(ax)
#display(fig)


@test (length(vec(attractors_SBD[1]))/10) > length(vec(attractors_reg[1]))

