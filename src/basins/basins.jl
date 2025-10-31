#########################################################################################
# Basin of Attraction structure definition
#########################################################################################


# All sub-types should contain the basins and attractors (maybe change to mapper??)
abstract type BasinsOfAttraction{ID} end

Base.iterate(BoA::BasinsOfAttraction, state=1) = state == 1 ? (BoA.basins, 2) : state == 2 ? (BoA.attractors, 3) : nothing

struct ArrayBasinsOfAttraction{ID, D, B <: AbstractArray{ID,D}, T, V, G <: Grid} <: BasinsOfAttraction{ID}
    basins::B 
    attractors::Dict{Int, StateSpaceSet{D, T, V}}
    grid::G
end

Base.size(BoA::ArrayBasinsOfAttraction, dim::Int) = size(BoA.basins, dim)

# Array basin with tuple grid
function ArrayBasinsOfAttraction(basins::B, attractors::Dict{Int, S}, grid_tup::Tuple) where {ID, B <: AbstractArray{ID}, S <: AbstractStateSpaceSet}
    if all(t -> t isa AbstractRange, grid_tup) && all(axis -> issorted(axis), grid_tup) # regular
        grid = RegularGrid(grid_tup)
    elseif any(t -> t isa AbstractVector, grid_tup) && all(axis -> issorted(axis), grid_tup) # irregular
        grid = IrregularGrid(grid_tup)
    else
        error("Incorrect grid specification!")
    end
    ArrayBasinsOfAttraction(basins, attractors, grid)
end

struct SampledBasinsOfAttraction{ID, D, T, V <: AbstractVector} <: BasinsOfAttraction{ID}
    basins::Vector{ID}
    attractors::Dict{Int, StateSpaceSet{D, T, V}}
    sampled_points::StateSpaceSet{D, T, V}
end

function SampledBasinsOfAttraction(basins::Vector{ID}, attractors::Dict{Int, StateSpaceSet{D, T, V}}, 
                                        sampled_points::Vector{U}) where {ID, D, T, V <: AbstractVector, U <: AbstractVector}
    S = StateSpaceSet(sampled_points)
    eltype(S) != V && error("The attractor points and sampled points must be represented" + 
    "by the same type of vector, that is they must have the same element type and length")
    SampledBasinsOfAttraction(basins, attractors, S)
end

Base.length(BoA::SampledBasinsOfAttraction) = length(BoA.basins)



#####################################################################################
# Pretty printing
#####################################################################################
function Base.show(io::IO, BoA::ArrayBasinsOfAttraction)
    ps = 14
    println(io, "$(nameof(typeof(BoA)))")
    println(io, rpad(" ID type: ", ps), typeof(BoA).parameters[1])
    println(io, rpad(" basin size: ", ps), size(BoA.basins))
    println(io, rpad(" grid: ", ps), BoA.grid)
    attstrings = split(sprint(show, MIME"text/plain"(), BoA.attractors), '\n')
    println(io, rpad(" attractors: ", ps), attstrings[1])
    for j in 2:size(attstrings,1)
        println(io, rpad(" ", ps), attstrings[j])
    end
    return
end

function Base.show(io::IO, BoA::SampledBasinsOfAttraction)
    ps = 14
    println(io, "$(nameof(typeof(BoA)))")
    println(io, rpad(" ID type: ", ps), typeof(BoA).parameters[1])
    println(io, rpad(" basin size: ", ps), length(BoA.basins))
    attstrings = split(sprint(show, MIME"text/plain"(), BoA.attractors), '\n')
    println(io, rpad(" attractors: ", ps), attstrings[1])
    for j in 2:size(attstrings,1)
        println(io, rpad(" ", ps), attstrings[j])
    end
    return
end

#########################################################################################
# Includes
#########################################################################################
include("fractality_of_basins.jl")
include("wada_test.jl")
