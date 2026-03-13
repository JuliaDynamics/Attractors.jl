using Test
@testset "Wada property test" begin
    # Artificial smooth boundary
    M = [1 1 1 1 2 2 2 2  3 3 3 3]' * ones(Int, 1, 12)
    @test test_wada_merge(M, 1) == 12 * 2 / 12^2
    @test test_wada_merge(M, 2) == 12 * 2 / 12^2

    # Artificial fractal boundary
    M = [1 2 3 1 2 3 1 2  3 1 2 3]' * ones(Int, 1, 12)
    @test test_wada_merge(M, 1) == 0.0
end
