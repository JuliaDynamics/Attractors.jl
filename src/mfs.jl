


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
    random_pertubation(mapper::AttractorMapper, X, search_area, dim,  total_iterations=10000)


Generate a random pertubation of the initial point "X" within specified "search_area" and checks if it is in the same basin of attraction.
It does so by generating a random vector of length dim and then adding it to X.
If the pertubation is not in the same basin of attraction, it calculates the norm of the pertubation and compares it to the best pertubation found so far.
If the norm is smaller, it updates the best pertubation found so far.
It repeats this process total_iterations times and returns the best pertubation found.

mapper: AttractorMapper
X: initial point to be tested
search_area: array of two points defining the search area
dim: dimension of the system
total_iterations: number of random pertubations to try, default = 10000

Output:
    best_shock: vector of the best shock found
    best_dist: norm of the best shock found
"""

function random_pertubation(mapper::AttractorMapper, X, search_area, dim, id_x;  total_iterations=10000)
    best_dist = Inf
    best_shock = nothing
    for _ in 1:total_iterations
        u0 = search_area[1] .+ rand(dim) .* (search_area[2] .- search_area[1])
        shock = X + u0
        if !(id_x == mapper(shock))
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
    mfs_brute_force(mapper, X, best_shock, best_dist, dim,  total_iterations=10000)

Input:
    mapper: AttractorMapper
    X: initial point to test
    best_shock: vector of the best shock found so far
    best_dist: norm of the best shock found so far
    dim: dimension of the system
    total_iterations: number of random pertubations to try, default = 10000

Output:
    best_shock: vector of the best shock found
    best_dist: norm of the best shock found

Description:
    This function works on the results obtained by random_pertubation(). It starts from the best shock found so far and tries to find a better one 
    by continuously reducing the radius of the sphere on the surface of which it generates random pertubations. If pertubation with the same basin of attraction is found,
    it updates the best shock found so far and reduces the radius of the sphere. 
    It repeats this process total_iterations times and returns the best pertubation found.
"""


function mfs_brute_force(mapper::AttractorMapper, X, best_shock, best_dist, dim, id_x,  total_iterations=10000)
    
    radius_dif = 100.

    for _ in 1:total_iterations

        pertubation = random_point_on_sphere((best_dist-best_dist/radius_dif), dim)
        new_shock = X + pertubation
        
        if !(id_x == mapper(new_shock))
            best_dist = LinearAlgebra.norm(pertubation)
            best_shock = pertubation
        end
    end

    return best_shock, best_dist
end




function mfs_objective(perturbation, X, idX, mapper::AttractorMapper, penalty=1000.0) 

    dist = LinearAlgebra.norm(perturbation)
    if dist == 0
        return penalty
    end
    if mapper(X + perturbation) == idX
        return dist + penalty
    else
        return dist
    end
end


"""
    minimal_fatal_shock(mapper, search_area, X, dim, total_iterations = 10000)

Input:
    mapper: AttractorMapper
    search_area: array of two points defining the search area
    X: initial point to test
    dim: dimension of the system
    total_iterations: number of iterations to run each step of the algorithm, default = 10000

Output:
    best_shock: vector of the best shock found
    best_dist: norm of the best shock found

Description:
    This function runs the minimal fatal shock algorithm on the initial point X.
    It does so by first generating a random pertubation of X and then trying to find a better one by using the brute force method.
"""


function minimal_fatal_shock(mapper::AttractorMapper, X, search_area; total_iterations = 10000, algorithm = :random)


    id_x = mapper(X)
    dim = mapper.ds

    function objective_function(perturbation)
        return mfs_objective(perturbation, X, id_x, mapper)
    end

    if algorithm == random 
        best_shock, best_dist = random_pertubation(mapper, X, search_area, dim, id_x;  total_iterations)
        best_shock, best_dist = mfs_brute_force(mapper, X, best_shock, best_dist, dim, id_x,  total_iterations)
    elseif algorithm == bboxopt
        result = bboptimize(objective_function; SearchRange = [search_area for _ in dim], NumDimensions = dim)
        best_shock = best_candidate(result)
        best_dist = LinearAlgebra.norm(best_shock)
    elseif algorithm == combined
        best_shock, best_dist = random_pertubation(mapper, X, search_area, dim, id_x;  total_iterations)
        best_shock, best_dist = mfs_brute_force(mapper, X, best_shock, best_dist, dim, id_x,  total_iterations)
        result = bboptimize(objective_function, best_shock; SearchRange = [search_area for _ in dim], NumDimensions = dim)
        best_shock = best_candidate(result)
        best_dist = LinearAlgebra.norm(best_shock)
    else
        error("Algorithm not recognized")
    end

    return best_shock, best_dist
end



