using Attractors
using CairoMakie




function newton_map(z, p, n)
    z1 = z[1] + im*z[2]
    dz1 = newton_f(z1, p[1])/newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
newton_f(x, p) = x^p - 1
newton_df(x, p)= p*x^(p-1)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
xg = yg = collect(range(-1.5, 1.5; length = 400))

xg = vcat(range(-1.5, 0, length = 300), range(-1.5, 0, length = 100))
yg = vcat(range(-1.5, 0, length = 300), range(-1.5, 0, length = 100))


newton = AttractorsViaRecurrences(ds, (xg, yg);
    sparse = false, mx_chk_lost = 1000, Dt =1,
)


basins, attractors = basins_of_attraction(newton)
attractors

grid = (xg, yg)
fig = heatmap_basins_attractors(grid, basins, attractors)

# keys(attractors)
# attractors
newton.bsn_nfo.grid_nfo.grid
newton([-0.5, -0.8])






function thomas_rule(u, p, t)
    x,y,z = u
    b = p[1]
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector{3}(xdot, ydot, zdot)
end

thomas_cyclical(u0 = [1.0, 0, 0]; b = 0.2) = CoupledODEs(thomas_rule, u0, [b])

ds = thomas_cyclical(b = 0.1665)
xg = yg = zg = range(-6.0, 6.0; length = 251)
mapper_3d = AttractorsViaRecurrences(ds, (xg, yg, zg); sparse = false)

#_, attractors = basins_of_attraction(mapper_3d)
attractors
collect(values(attractors))[2].data[1]
mapper_3d(collect(values(attractors))[2].data[1])

xg = yg = zg = collect(range(-6.0, 6.0; length = 251))
mapper = AttractorsViaRecurrences(ds, (xg, yg, zg); sparse = false, Dt = 1)
mapper(collect(values(attractors))[6].data[1])
