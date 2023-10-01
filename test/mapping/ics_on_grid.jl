using Attractors, Test

@testset "vary 2 dimensions in higher dimension system" begin
    function _construct_grid(idx_unit_1, idx_unit_2, ics_unit_1, ics_unit_2, ics_base_units)
        _grid = [ic:ic for ic in ics_base_units]
        _grid[idx_unit_1] = ics_unit_1
        _grid[idx_unit_2] = ics_unit_2
        return Tuple(_grid)
    end

    size_grid = 100
    num_ics_1 = 100
    num_ics_2 = 10

    idx_dim_1 = 1 
    idx_dim_2 = 2
    ics_dim_1 = range(1, 10; length=num_ics_1)
    ics_dim_2 = range(1, 10; length=num_ics_2)
    ics_base = Float64.(collect(1:size_grid))

    grid = _construct_grid(idx_dim_1, idx_dim_2, ics_dim_1, ics_dim_2, ics_base)

    A = Attractors.generate_ics_on_grid(grid)

    @test length(A) == num_ics_1 * num_ics_2
    @test all(ics_dim_1 .== [A[idx][idx_dim_1] for idx in 1:num_ics_1])
    @test all(ics_dim_2[1] .== [A[idx][idx_dim_2] for idx in 1:num_ics_1])
    @test all(ics_dim_2[2] .== [A[idx][idx_dim_2] for idx in (num_ics_1+1):2*num_ics_1])

end
