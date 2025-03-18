DO_EXTENSIVE_TESTS = get(ENV, "ATTRACTORS_EXTENSIVE_TESTS", "false") == "true"

using Test, Attractors
using Random

# use this dumb map, map initial conditiosn and then call finalize.
# you should know analytically the value for all nonlocal stability measures.
function dumb_map(dz, z, p, n)
    x, y = z
    r = p[1]
    if r < 0.5
        dz[1] = dz[2] = 0.0
    else
        if x > 0
            dz[1] = r
            dz[2] = r
        else
            dz[1] = -r
            dz[2] = -r
        end
    end
    return
end