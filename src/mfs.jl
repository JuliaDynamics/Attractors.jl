
""" 
    check_basin(mapper, X, Y)

Input:
    mapper: AttractorMapper
    X: Point 1
    Y: Point 2
Output: 
    true if X and Y are in the same basin of attraction, false otherwise

Description:
    This function checks if two points are in the same basin of attraction. 
    It does so by checking if the points are mapped to the same point by the mapper.

"""

function check_basin(mapper, X, Y)
    idx = mapper(X)
    idy = mapper(Y)
    if idx == idy
        return true
    else
        return false
    end
end


""" 
    random_point_on_sphere(radius, n)

Input:
    radius: radius of the sphere
    n: dimension of the sphere
Output:
    coordinates of a random point on the n-sphere

Description:
    This function generates a random point on the n-sphere of radius radius.
    It does so by generating n-1 random angles and then calculating the last coordinate
    using the formula for the n-sphere.
"""


function random_point_on_sphere(radius, n)
    angles = rand(n-1) * 2π  

    coordinates = [radius * prod(sin.(angles[1:i-1])) * cos(angles[i]) for i in 1:n-1]
    push!(coordinates, radius * prod(sin.(angles)))  
    
    return coordinates
end

""" 
    random_pertubation(mapper, X, search_area, dim,  Ntries=10000)

Input:
    mapper: AttractorMapper
    X: initial point to test
    search_area: array of two points defining the search area
    dim: dimension of the system
    Ntries: number of random pertubations to try, default = 10000

Output:
    best_shock: vector of the best shock found
    best_dist: norm of the best shock found

Description:
    This function generates a random pertubation of the initial point "X" within specified "search_area" and checks if it is in the same basin of attraction.
    It does so by generating a random vector of length dim and then adding it to X.
    If the pertubation is not in the same basin of attraction, it calculates the norm of the pertubation and compares it to the best pertubation found so far.
    If the norm is smaller, it updates the best pertubation found so far.
    It repeats this process Ntries times and returns the best pertubation found.

"""

function random_pertubation(mapper::AttractorMapper, X, search_area, dim,  Ntries=10000)
    best_dist = Inf
    best_shock = nothing
    Xsp = StateSpaceSet([X])
    for _ in 1:Ntries
        u0 = search_area[1] .+ rand(dim) .* (search_area[2] .- search_area[1])
        shock = X + u0
        if !(check_basin(mapper, shock, X))
            dist = LinearAlgebra.norm(u0)
            if dist < best_dist
                best_dist = dist
                best_shock = u0
            end
        end
    end
    return best_shock, best_dist
end


""" 
    mfs_brute_force(mapper, X, best_shock, best_dist, dim,  Ntries=10000)

Input:
    mapper: AttractorMapper
    X: initial point to test
    best_shock: vector of the best shock found so far
    best_dist: norm of the best shock found so far
    dim: dimension of the system
    Ntries: number of random pertubations to try, default = 10000

Output:
    best_shock: vector of the best shock found
    best_dist: norm of the best shock found

Description:
    This function works on the results obtained by random_pertubation(). It starts from the best shock found so far and tries to find a better one 
    by continuously reducing the radius of the sphere on the surface of which it generates random pertubations. If pertubation with the same basin of attraction is found,
    it updates the best shock found so far and reduces the radius of the sphere. 
    It repeats this process Ntries times and returns the best pertubation found.
"""


function mfs_brute_force(mapper, X, best_shock, best_dist, dim,  Ntries=10000)
    
    radius_dif = 100.

    for _ in 1:Ntries

        pertubation = random_point_on_sphere((best_dist-best_dist/radius_dif), dim)
        new_shock = X + pertubation
        
        if !(check_basin(mapper, X, new_shock))
            best_dist = set_distance(StateSpaceSet([X]), StateSpaceSet([new_shock]), Centroid())
            best_shock = pertubation
        end
    end

    return best_shock, best_dist
end


"""
    mfs(mapper, search_area, X, dim, n_iterations = 10000)

Input:
    mapper: AttractorMapper
    search_area: array of two points defining the search area
    X: initial point to test
    dim: dimension of the system
    n_iterations: number of iterations to run each step of the algorithm, default = 10000

Output:
    best_shock: vector of the best shock found
    best_dist: norm of the best shock found

Description:
    This function runs the minimal fatal shock algorithm on the initial point X.
    It does so by first generating a random pertubation of X and then trying to find a better one by using the brute force method.
"""


function mfs(mapper::AttractorMapper, search_area, X, dim, n_iterations = 10000)
    best_shock, best_dist = random_pertubation(mapper, X, search_area, dim,  n_iterations)
    best_shock, best_dist = mfs_brute_force(mapper, X, best_shock, best_dist, dim, n_iterations)

    return best_shock, best_dist
end

