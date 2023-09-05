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