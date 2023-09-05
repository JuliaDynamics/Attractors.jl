
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



matrices = []
for _ in 1:3
    matrix = zeros(Int, rand((2,3,4)), rand((2,3,4)))
    for i in eachindex(matrix)
        matrix[i] = rand((0,1,2,3))
    end

    function newton_map(z, p, n)
        z1 = z[1] + im*z[2]
        dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
        z1 = z1 - dz1
        return SVector(real(z1), imag(z1))
    end
    newton_f(x, p) = x^p - 1
    newton_df(x, p)= p*x^(p-1)

    ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
    xg = yg = range(-1.5, 1.5, length = 100)

    newton = AttractorsViaRecurrences(ds, (xg, yg);
        sparse = false, mx_chk_lost = 1000, Dt =1, density_matrix = matrix,
    )

    @test ((newton([-0.5, 0.86]) != newton([-0.5, -0.86]))& (newton([-0.5, 0.86]) != newton([1.0, 0.0])) & (newton([-0.5, -0.86]) != newton([1.0, 0.0])))

    basins, attractors = basins_of_attraction(newton)

    @test length(attractors) == 3
end





# function get_dimension(matrix)
#     dim = [size(matrix)[1]]
#     while !(matrix[1] isa Int)  
#         println(dim)
#         matrix = matrix[1]
#         push!(dim, size(matrix)[1])
#     end


#     return dim
# end



# function collect_all_grids(grid, matrix)
#     unique = Set([0,1,2])
#     grid_steps = Dict()
#     for i in unique
#         grid_steps[i] = [length(axis)*2^i for axis in grid]
#         #grids[i] = [range(first(axis), last(axis), length = length(axis)*2^i) for axis in grid]
#     end

#     return grid_steps
# end

# function point_to_index(grid_nfo, matrix,  point)
#     dims = get_dimension(matrix)
#     grid_minima = grid_nfo.grid_minima
#     grid_maxima = grid_nfo.grid_maxima
#     ratios = []
#     for i in 1:length(dims)
#         if i == 1
#             ratio = abs(grid_minima[i] - point[i])/abs(grid_minima[i] - grid_maxima[i])
#         else
#             ratio = abs(grid_maxima[i] - point[i])/abs(grid_minima[i] - grid_maxima[i])
#         end
#         push!(ratios, ratio )

#     end
    
#     index = zeros(length(dims))

#     for i in 1:length(dims)
#         for j in 1:dims[i]
#             index[i] += 1/dims[i]
#             if index[i] >= ratios[i]
#                 index[i] = Int(j)
#                 break
#             end
#         end
#     end
#     println(index)
#     return matrix[Int(index[2])][Int(index[1])]
# end

# matrix = [[1,2], [1, 0]]

# grid = (range(-1.5, 1.5, length = 3),range(-1.5, 1.5, length = 3))
# grid_steps = collect_all_grids(grid, matrix)
# collect(range(-1.5, 1.5, length = 3))
# grid_info = Attractors.IrregularGridViaMatrix(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), matrix, (range(-1.5, 1.5, length = 80),range(-1.5, 1.5, length = 80)))








point_to_index(grid_info, matrix, (1,1))

matrix = [[1,2], [1, 0]]


using Attractors
using Test

function collect_all_grids(grid, matrix)
    unique = Set(matrix)
    grid_steps = Dict()
    for i in unique
        grid_steps[i] = [length(axis)*2^i for axis in grid]
        #grids[i] = [range(first(axis), last(axis), length = length(axis)*2^i) for axis in grid]
    end

    return grid_steps
end

function point_to_index(grid_nfo, matrix,  point)
    size_y, size_x = size(matrix)
    grid_minima = grid_nfo.grid_minima
    grid_maxima = grid_nfo.grid_maxima
    step_x = abs(grid_minima[1] - point[1])
    step_y = abs(grid_maxima[2]- point[2])

    ratio_x = step_x/abs(grid_minima[1] - abs(grid_maxima[1]))
    ratio_y = step_y/abs(grid_minima[2] - grid_maxima[2])
    temp_x = 0
    temp_y = 0
    for i in 1:size_x
        temp_x += 1/size_x
        if temp_x >= ratio_x
            temp_x = i
            break
        end
    end 
    for j in 1:size_y
        temp_y += 1/size_y
        if temp_y >= ratio_y
            temp_y = j
            break
        end
    end 
    return matrix[Int(temp_y), Int(temp_x)]
end

function basin_cell_index(y_grid_state, grid_nfo)
    
    cell_area = point_to_index(grid_nfo, grid_nfo.matrix, y_grid_state) ## correct first
    
    grid_maxima = grid_nfo.grid_maxima
    grid_minima = grid_nfo.grid_minima
    max_level = maximum(grid_nfo.matrix)
    println(grid_nfo.grid_steps)
    grid_steps = abs.(grid_minima - grid_maxima)./ grid_nfo.grid_steps[cell_area]
    iswithingrid = true
    @inbounds for i in eachindex(grid_minima)

        if !(grid_minima[i] ≤ y_grid_state[i] ≤ grid_maxima[i])
            iswithingrid = false
            break
        end
    end
    if iswithingrid
        # Snap point to grid
        # rough_ind = @. (y_grid_state - grid_minima)/grid_steps 
        # ind = []
        # for i in rough_ind
        #         (i % 1 <= 0.5) ? push!(ind, 1) : push!(ind, 0)
        # end
        # println(ind)
        # indices = @. round(Int, rough_ind ) *(2^(max_level-cell_area)) + ind


        indices = @. round(Int,(y_grid_state - grid_minima)/grid_steps, RoundDown) * (2^(max_level-cell_area)) + 1

        # indices = []
        
        # for ind in rough_indices
        #     (ind % 1 < 0.5) ? push!(indices, round(Int, ind+1)) : push!(indices, round(Int, ind))
        # end
        # println(indices)
        return CartesianIndex(indices...) 

    else
        return CartesianIndex(-1)
    end
end


matrix = [0 3; 3 0]
grid = (range(0, 10, length = 11),range(0, 10, length = 11))
grid_steps = collect_all_grids(grid, matrix)
# collect(range(0, 9, length = 10))
# collect(range(0, 9, length = 20))[2:8]
# collect(range(0, 9, length = 40))[1:6]
# collect(range(0, 9, length = 80))[10:13]

#grid_info = Attractors.IrregularGridViaMatrix(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), matrix, grid)

matrix = [0 3; 3 0]
grid = Attractors.IrregularGridViaMatrix(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), matrix, grid)
@test basin_cell_index((1.25, 1.25), grid) == CartesianIndex(12,12)

grid = (range(0, 10, length = 11),range(0, 10, length = 11))
grid_steps = collect_all_grids(grid, matrix)
matrix = [0 3; 2 0]
grid = Attractors.IrregularGridViaMatrix(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), matrix, grid)
@test basin_cell_index((1.25, 1.25), grid) == CartesianIndex(11, 11)

grid = (range(0, 10, length = 11),range(0, 10, length = 11))
grid_steps = collect_all_grids(grid, matrix)
matrix = [0 3; 0 0]
grid = Attractors.IrregularGridViaMatrix(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), matrix, grid)
@test basin_cell_index((1.25, 1.25), grid) == CartesianIndex(9, 9)













matrix = [0 0; 1 0]
grid = (range(-10, 10, length = 10),range(-10, 10, length = 10))
grid_steps = collect_all_grids(grid, matrix)
CartesianIndices(collect(range(0, 10, length = 10)))
collect(range(0, 10, length = 20))
collect(range(0, 10, length = 40))[1:10]
grid_info = Attractors.IrregularGridViaMatrix(grid_steps, SVector{2, Float64}(minimum.(grid)), SVector{2, Float64}(maximum.(grid)), matrix, grid)
grid_info
basin_cell_index((-10,-10), grid_info)



vect = [ [1, 1], [1, 1]]

A = rand(3,2)
CartesianIndices(A)

a = [-1 0 1; -1 0 1]
CartesianIndices(a)

round(1.4)
typeof(-1)

