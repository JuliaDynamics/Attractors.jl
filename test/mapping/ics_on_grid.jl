using Attractors, Test

@testset "Generating ics" begin
    function _construct_grid_varying_2_dims(idx_dim_1, idx_dim_2, ics_dim_1, ics_dim_2, ics_base_dims)
        _grid = [ic:ic for ic in ics_base_dims]
        _grid[idx_dim_1] = ics_dim_1
        _grid[idx_dim_2] = ics_dim_2
        return Tuple(_grid)
    end
    
    function construct_grid_varying_2_dims(size_grid, num_ics_1, num_ics_2, idx_dim_1, idx_dim_2) #silly function just to avoid repeating code
        ics_dim_1 = range(1, 10; length=num_ics_1)
        ics_dim_2 = range(1, 10; length=num_ics_2)
        ics_base = Float64.(collect(1:size_grid))
        grid = _construct_grid_varying_2_dims(idx_dim_1, idx_dim_2, ics_dim_1, ics_dim_2, ics_base)
        return grid, ics_dim_1, ics_dim_2
    end
    
   @testset "Varying in 10 dimensions" begin
        xg = range(1, 10; length=3)
        grid = ntuple(x->xg, 10)
        A = @time Attractors._generate_ics_on_grid_vary_few_dims(grid) 
        basins = zeros(Int32, map(length, grid))
        A_previous = @time Attractors._generate_ics_on_grid(grid, basins)
        @test A == A_previous
    end
    
    @testset "Varying in 2 dimensions - small" begin
        grid, _, _ = construct_grid_varying_2_dims(10, 3, 4, 1, 2)
        A = @time Attractors._generate_ics_on_grid_vary_few_dims(grid)
        basins = zeros(Int32, map(length, grid))
        A_previous = @time Attractors._generate_ics_on_grid(grid, basins) 
    end

    @testset "Varying in 2 dimensions - big" begin
        num_ics_1 = 100; num_ics_2 = 10; idx_dim_1 = 1; idx_dim_2 = 2
        grid, ics_dim_1, ics_dim_2 = construct_grid_varying_2_dims(100, num_ics_1, num_ics_2, idx_dim_1, idx_dim_2)
        A = @time Attractors._generate_ics_on_grid_vary_few_dims(grid)
        @test length(A) == num_ics_1 * num_ics_2
        @test all(ics_dim_1 .== [A[idx][idx_dim_1] for idx in 1:num_ics_1])
        @test all(ics_dim_2[1] .== [A[idx][idx_dim_2] for idx in 1:num_ics_1])
        @test all(ics_dim_2[2] .== [A[idx][idx_dim_2] for idx in (num_ics_1+1):2*num_ics_1])
    end

end

