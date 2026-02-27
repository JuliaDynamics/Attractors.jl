using Neighborhood
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
    BasinsOfAttraction

A subtype of `BasinsOfAttraction` is a convenient structure that stores a representation
of the `basins` of attraction, their associated `attractors`, and a representation of the domain
over which the basin is defined. For example, this domain could be a `Grid` subtype matching
the size of the basins as an array or a set of points sampled from the state space. These fields
can be accessed using the [`extract_basins`](@ref), [`extract_attractors`](@ref),
and [`extract_domain`](@ref) functions respectively.

Currently available subtypes:

* [`ArrayBasinsOfAttraction`](@ref)
* [`SampledBasinsOfAttraction`](@ref)

All `BasinsOfAttraction` subtypes can be used with [`basins_fractions`](@ref) provided
that the basins are represented as subtypes of `AbstractArray`. Additionally, all
`BasinsOfAttraction` subtypes are iterable in the sense: `basins, attractors = BoA`, this
was done to ensure backwards compatibility for functions whose original return format was
`basins, attractors` but has since been replaced with a `BasinsOfAttraction` type.

The [`map_to_basin`](@ref) function provides simple interpolation of a point in state space
to determine which basin of attraction it is likely to belong to.
"""
abstract type BasinsOfAttraction{ID} end


"""
    ArrayBasinsOfAttraction(basins, attractors, grid)

A subtype of [`BasinsOfAttraction`](@ref) whose `basins` of attraction are represented by an
`array::AbstractArray`, that has `D` number of dimensions. The `attractors` take the form of a
dictionary mapping attractor labels to `StateSpaceSet`'s with the points of each set being of
length `D`. The `grid` represents the spatial domain, and can be anything given to [`AttractorsViaRecurrences`](@ref)
as a grid, i.e., a tuple of ranges or a `Grid` type.

"""
struct ArrayBasinsOfAttraction{ID, D, B <: AbstractArray{ID, D}, G <: Grid, K, S <: StateSpaceSet} <: BasinsOfAttraction{ID}
    basins::B
    attractors::Dict{K, S}
    grid::G

    function ArrayBasinsOfAttraction(basins::AbstractArray{ID, D}, attractors::Dict{K, S}, grid::G) where {D, ID, G <: Grid, K, S <: StateSpaceSet}
        # Dimensionality checks
        length(grid.grid) != ndims(basins) && error("The basins and the grid must have the same number of dimensions")
        # Attractor state space sets have the same type so same dimensions, can compare grid with any of them
        if !isempty(attractors) #
            length(grid.grid) != length(valtype(collect(values(attractors))[1])) && error("The attractor points and the grid must have the same number of dimensions")
        end
        B = typeof(basins)
        return new{ID, D, B, G, K, S}(basins, attractors, grid)
    end
end
# The definition of other constructors can be found in `basins/basins_utilities.jl`.

"""
    SampledBasinsOfAttraction(basins, attractors, sampled_points)

A subtype of [`BasinsOfAttraction`](@ref) whose `basins` of attraction are represented by a `vector::AbstractVector`.
The `attractors` take the form of a dictionary mapping attractor labels to `StateSpaceSet`'s with
the points of each set being of equal length. The spatial domain of this basin type is `sampled_points` which
can be a `StateSpaceSet` with the same dimensionality, element type, and vector type as those
used to represent the attractors or alternatively a vector of points with the aforementioned requirements.

Additional keyword arguments may be specified for use in the construction of a search structure which [`map_to_basin`](@ref)
uses to interpolate state space points to their nearest basin. These arguments are:

* `tree`: search tree constructor (e.g. `KDTree`, `BallTree`)
* `metric`: distance metric (e.g. `Euclidean()`, `Chebyshev()`)
* `searchstructure_kwargs...`: additional keyword arguments passed to `searchstructure`
"""
struct SampledBasinsOfAttraction{ID, D, T, V <: AbstractVector, AK, S <: StateSpaceSet{D, T, V}, ss <: Neighborhood.SearchType} <: BasinsOfAttraction{ID}
    points_ids::Vector{ID}
    attractors::Dict{AK, S}
    sampled_points::S
    search_struct::ss

    function SampledBasinsOfAttraction(basins::Vector{ID}, attractors::Dict{AK, S}, sampled_points::S; tree = KDTree, metric = Euclidean(), ss_kwargs...) where
        {ID, D, T, V <: AbstractVector, AK, S <: StateSpaceSet{D, T, V}}
        # Dimensionality checks
        length(basins) != sampled_points && error("The basins and the sampled points must have equal length")
        search_struct = searchstructure(tree, BoA.sampled_points, metric, ss_kwargs...)
        return new{ID, D, T, V, AK, S, typeof(search_struct)}(basins, attractors, sampled_points, search_struct)
    end
end
# The definition of other constructors can be found in `basins/basins_utilities.jl`.

#########################################################################################
# Basins of Attraction Convenience functions
#########################################################################################
Base.iterate(BoA::BasinsOfAttraction, state = 1) = state == 1 ? (extract_basins(BoA), 2) : state == 2 ? (extract_attractors(BoA), 3) : nothing

"""
    extract_basins(BoA::BasinsOfAttraction) → basins

Returns the basins component of a `BasinsOfAttraction` object.
"""
extract_basins(BoA::ArrayBasinsOfAttraction) = BoA.basins
extract_basins(BoA::SampledBasinsOfAttraction) = BoA.points_ids

"""
    extract_attractors(BoA::BasinsOfAttraction) → attractors

Returns the attractors component of a `BasinsOfAttraction` object. Which is a dictionary mapping
attractor labels to attractors represented as `StateSpaceSet`'s.
"""
extract_attractors(BoA::ArrayBasinsOfAttraction) = BoA.attractors
extract_attractors(BoA::SampledBasinsOfAttraction) = BoA.attractors

"""
    extract_domain(BoA::BasinsOfAttraction) → domain

Returns the domain component of a `BasinsOfAttraction` object.
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
    for j in 2:size(attstrings, 1)
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
    for j in 2:size(attstrings, 1)
        println(io, rpad(" ", ps), attstrings[j])
    end
    return
end
