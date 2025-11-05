export BasinsOfAttraction, 
    ArrayBasinsOfAttraction, 
    SampledBasinsOfAttraction,
    extract_basins,
    extract_basins,
    extract_domain
    
#########################################################################################
# Basins of Attraction structure definition
#########################################################################################
"""
    BasinsOfAttraction{ID}

A subtype of `BasinsOfAttraction` is a convenient structure that stores a representation
of the `basins` of attraction, their associated `attractors`, and a representation of the `domain` 
over which the basin is defined. For example, this domain could be a `Grid` subtype matching
the size of the basins as an array or a set of points sampled from the state space. These fields 
can be accessed using the [`extract_basins`](@ref), [`extract_attractors`](@ref), 
and [`extract_domain`](@ref) functions respectively

The type parameter `ID` specifies the type of values stored in the basins. Typically,
this would be `Int`, corresponding to integer valued attractor labels.

Currently available subtypes:

* [`ArrayBasinsOfAttraction`](@ref)
* [`SampledBasinsOfAttraction`](@ref)

All `BasinsOfAttraction` subtypes can be used with [`basins_fractions`](@ref) provided
that the basins are represented as subtypes of `AbstractArray`. Additionally, all 
`BasinsOfAttraction` subtypes are iterable in the sense: `basins, attractors = BoA`, this 
was done to ensure backwards compatibility for functions whose original return format was
`basins, attractors` but has since been replaced with a `BasinsOfAttraction` type.

## For developers

`BasinsOfAttraction{ID}` defines an extendable interface. A new type needs to subtype
`BasinsOfAttraction{ID}` and implement [`extract_basins`](@ref), [`extract_attractors`](@ref),
and [`extract_domain`](@ref).
"""
abstract type BasinsOfAttraction{ID} end


"""
    ArrayBasinsOfAttraction{ID, D} <: BasinsOfAttraction{ID}
    ArrayBasinsOfAttraction(basins, attractors, grid::Tuple)
    ArrayBasinsOfAttraction(basins, attractors, grid::Grid)

A subtype of [`BasinsOfAttraction`](@ref) whose `basins` of attraction are represented by an array, 
or more specifically a subtype of `AbstractArray` whose values have type `ID` and has `D` 
number of dimensions. The `attractors` take the form of a dictionary mapping attractor labels
to `StateSpaceSet`'s with the points of each set being of length `D`. The domain here is `grid`
and has type `Grid`, each grid "cell" has a unique corresponding entry in the basins. 

The first constructor simply converts the tuple representation of the grid into type `Grid`, 
and then calls the second constructor. This `grid` tuple is a tuple of ranges defining the grid 
of initial conditions that partition the state space. The `grid` has to be of the same dimensionality
as the points representing the attractors and of the `basins`. 
"""
struct ArrayBasinsOfAttraction{ID, D, B <: AbstractArray{ID,D}, T, V, G <: Grid, AK} <: BasinsOfAttraction{ID}
    basins::B 
    attractors::Dict{AK, StateSpaceSet{D, T, V}}
    grid::G
end
# The definition of the first constructor can be found in `basins/basins_utilities.jl`.

"""
    SampledBasinsOfAttraction{ID, D} <: BasinsOfAttraction{ID}
    SampledBasinsOfAttraction(basins, attractors, sampled_points::Vector{U}) where {U <: AbstractVector}
    SampledBasinsOfAttraction(basins, attractors, sampled_points::S) where {S <: AbstractStateSpaceSet}

A subtype of [`BasinsOfAttraction`](@ref) whose `basins` of attraction are represented by a vector
whose values are of type `ID`. The `attractors` take the form of a dictionary mapping attractor 
labels to `StateSpaceSet`'s with the points of each set being of length `D`. The domain of this 
basin type is `sampled_points` which is also a `StateSpaceSet` with the same dimensionality `D`, 
element type, and vector type as those used to represent the attractors. Each sampled point has a 
unique corresponding entry in the `basins` vector. The length of `basins` and `sampled_points` must
be equal.

The first constructor converts a vector of vectors into a `StateSpaceSet`, then calls the second constructor.
That is, assuming each has the same dimensionality `D`, and the same element type as the `StateSpaceSet`'s 
representing the attractors.
"""
struct SampledBasinsOfAttraction{ID, D, T, V <: AbstractVector, AK} <: BasinsOfAttraction{ID}
    points_ids::Vector{ID}
    attractors::Dict{AK, StateSpaceSet{D, T, V}}
    sampled_points::StateSpaceSet{D, T, V}
end
# The definition of the first constructor can be found in `basins/basins_utilities.jl`.

#########################################################################################
# Basins of Attraction Convenience functions
######################################################################################### 
Base.iterate(BoA::BasinsOfAttraction, state=1) = state == 1 ? (extract_basins(BoA), 2) : state == 2 ? (extract_attractors(BoA), 3) : nothing

"""
    extract_basins(BoA::BasinsOfAttraction) → basins

Returns the basins component of a `BasinsOfAttraction`. For developing a new basin subtype
this function should be extended.
"""
extract_basins(BoA::ArrayBasinsOfAttraction) = BoA.basins
extract_basins(BoA::SampledBasinsOfAttraction) = BoA.points_ids

"""
    extract_attractors(BoA::BasinsOfAttraction) → attractors

Returns the attractors component of a `BasinsOfAttraction`. Which is a dictionary mapping
attractor labels to attractors represented as `StateSpaceSet`'s.

For developing a new basin subtype
this function should be extended.
"""
extract_attractors(BoA::ArrayBasinsOfAttraction) = BoA.attractors
extract_attractors(BoA::SampledBasinsOfAttraction) = BoA.attractors

"""
    extract_domain(BoA::BasinsOfAttraction) → domain 

Returns the domain component of a `BasinsOfAttraction`. For developing a new basin subtype
this function should be extended.
"""
extract_domain(BoA::SampledBasinsOfAttraction) = BoA.sampled_points
extract_domain(BoA::ArrayBasinsOfAttraction) = BoA.grid

#####################################################################################
# Pretty printing
#####################################################################################
function Base.show(io::IO, BoA::ArrayBasinsOfAttraction)
    ps = 14
    println(io, "$(nameof(typeof(BoA)))")
    println(io, rpad(" ID type: ", ps), typeof(BoA).parameters[1])
    println(io, rpad(" basin size: ", ps), size(BoA.basins))
    println(io, rpad(" grid: ", ps), extract_domain(BoA))
    attstrings = split(sprint(show, MIME"text/plain"(), extract_attractors(BoA)), '\n')
    println(io, rpad(" attractors: ", ps), attstrings[1])
    for j in 2:size(attstrings,1)
        println(io, rpad(" ", ps), attstrings[j])
    end
    return
end

# Add first few vectors to printing followed by ...
function Base.show(io::IO, BoA::SampledBasinsOfAttraction)
    ps = 14
    println(io, "$(nameof(typeof(BoA)))")
    println(io, rpad(" ID type: ", ps), typeof(BoA).parameters[1])
    println(io, rpad(" basin length: ", ps), length(BoA))
    attstrings = split(sprint(show, MIME"text/plain"(), extract_attractors(BoA)), '\n')
    println(io, rpad(" attractors: ", ps), attstrings[1])
    for j in 2:size(attstrings,1)
        println(io, rpad(" ", ps), attstrings[j])
    end
    return
end
