using Attractors, Test

# This is a fake bistable map that has two equilibrium points
# for r > 0.5. It has predictable fractions.
function dumb_map(z, p, n)
    x, y = z
    r = p[1]
    if r < 0.5
        return SVector(0.0, 0.0)
    else
        if x ≥ 0
            return SVector(r, r)
        else
            return SVector(-r, -r)
        end
    end
end

r = 1.0
ds = DeterministicIteratedMap(dumb_map, [0., 0.], [r])
# deterministic grid, we know exactly the array layout
xg = yg = range(-1.5, 2.5; length=3)
grid = (xg, yg)
attrs = Dict(1 => StateSpaceSet([SVector(r, r)]), 2 => StateSpaceSet([SVector(-r, -r)]))
mapper = AttractorsViaProximity(ds, attrs; Ttr = 0)

basins, atts = basins_of_attraction(mapper, grid; show_progress=false)

@test basins[1, :] == fill(2, 3)
@test basins[2, :] == fill(1, 3)
@test basins[3, :] == fill(1, 3)



basins, atts, iterations = iterates_of_basins_of_attraction(mapper, grid; show_progress=false)

@test basins[1, :] == fill(2, 3)
@test basins[2, :] == fill(1, 3)
@test basins[3, :] == fill(1, 3)
@test length(atts) == 2
@test iterations == fill(1,3,3) 



