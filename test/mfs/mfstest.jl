using Test
using Attractors
using Distributions

function newton_map(z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
xg = yg = range(-1.5, 1.5; length = 400)

newton = AttractorsViaRecurrences(ds, (xg, yg);
    sparse = false, mx_chk_lost = 1000
)


attractors = [[1.0, 0.0], [-0.5, 0.8660254037844386], [-0.5, -0.8660254037844386]]
algo_r = Attractors.MFSBruteForce()
randomised = Dict([atr => Attractors.minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_r) for atr in attractors])


random_seed = [rand(Uniform(-0.5, 0.5), 2) for _ in 1:20]
randomised_r = Dict([atr => Attractors.minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_r) for atr in random_seed])
blackbox_r = Dict([atr => Attractors.minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_bb) for atr in random_seed])



@testset "Newton 2d" begin
    @testset begin
        test = true
        for i in (keys(randomised)) 
            
            if randomised[i][2] >= 0.63 || randomised[i][2] <= 0.62 || newton(randomised[i][1] + i) == newton(i)
                test = false
            end
        end
        @test test
    end
        
        
    algo_bb = Attractors.MFSBlackBoxOptim(dimension = 2)
    blackbox = Dict([atr => Attractors.minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_bb) for atr in attractors])

    @testset begin
        test = true
        for i in (keys(blackbox)) 
            
            if blackbox[i][2] >= 0.629 || blackbox[i][2] <= 0.62009 || newton(blackbox[i][1] + i) == newton(i)
                test = false
            end
        end
        @test test
    end



    @testset begin
        test = true
        for i in (keys(randomised_r)) 
            println(randomised_r[i][2] <= 0.5)
            if randomised_r[i][2] >= 0.5 || newton(randomised_r[i][1] + i) == newton(i)
                test = false
            end
        end
        @test test
    end


    @testset begin
        test = true
        for i in (keys(blackbox_r)) 
            
            if blackbox_r[i][2] >= 0.5 || newton(blackbox_r[i][1] + i) == newton(i)
                test = false
            end
        end
        @test test
    end
end