using Test
using Attractors
using Distributions



###############################################
#           Newton 2D fractal setup           #
###############################################
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

algo_bb = Attractors.MFSBlackBoxOptim()
blackbox = Dict([atr => Attractors.minimal_fatal_shock(newton, atr, [(-1.5, 1.5)], algo_bb) for atr in attractors])

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



###############################################
#           Magnetic 2D                       #
###############################################




ds = Systems.magnetic_pendulum(d=0.2, α=0.2, ω=0.8, N=3)

psys = ProjectedDynamicalSystem(ds, [1, 2], [0.0, 0.0])

attractors_m = Dict(i => StateSpaceSet([dynamic_rule(ds).magnets[i]]) for i in 1:3)

mapper_m = AttractorsViaProximity(psys, attractors_m)

xg = yg = range(-4, 4; length = 201)
grid = (xg, yg)
basins, attractors_m = basins_of_attraction(mapper_m, grid; show_progress = false)

attractor3 = vec((collect(values(attractors_m)))[3])
attractor2 = vec((collect(values(attractors_m)))[2])
attractor1 = vec((collect(values(attractors_m)))[1])


randomised_r = Dict([atr => Attractors.minimal_fatal_shock(mapper_m, atr, [(-4, 4), (-4, 4)], algo_r) 
                                         for atr in [attractor1[1], attractor2[1], attractor3[1]]])


@testset "Magnetic 2D" begin
    
    @test map(x -> (x[2] <= 0.4) && (x[2]) > 0.39, values(randomised_r)) |> all


    blackbox_r = Dict([atr => Attractors.minimal_fatal_shock(mapper, atr, [(-4, 4), (-4, 4)], algo_bb) 
                                            for atr in [attractor1[1], attractor2[1], attractor3[1]]])


    @test map(x -> (x[2] <= 0.395) && (x[2]) > 0.39, values(blackbox_r)) |> all
end





