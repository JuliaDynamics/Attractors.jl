"""
    basins_fractions(basins::AbstractArray) → fs::Dict

Calculate the state space fraction of the basins of attraction encoded in `basins`.
The elements of `basins` are integers, enumerating the attractor that the entry of
`basins` converges to (i.e., like the output of [`basins_of_attraction`](@ref)).
Return a dictionary that maps attractor IDs to their relative fractions.

In [Menck2013](@cite) the authors use these fractions to quantify the stability of a basin of
attraction, and specifically how it changes when a parameter is changed.
For this, see [`continuation`](@ref).
"""
function basins_fractions(basins::AbstractArray, ids = unique(basins))
    fs = Dict{eltype(basins), Float64}()
    N = length(basins)
    for ξ in ids
        B = count(isequal(ξ), basins)
        fs[ξ] = B/N
    end
    return fs
end

"""
    convergence_and_basins_of_attraction(mapper::AttractorMapper, grid)

An extension of [`basins_of_attraction`](@ref).
Return `basins, attractors, iterations`, with `basins, attractors` as in
`basins_of_attraction`, and `iterations` being an array with same shape
as `basins`. It contains the number of iterations each initial condition took
to converge to its attractor. Multiply this with the time step timescale
`Δt` given to `mapper` to obtain normal time units.
`iterations` is useful to give to [`shaded_basins_heatmap`](@ref).

See also [`iterations_to_converge`](@ref).

# Keyword arguments

- `show_progress = true`: show progress bar.
"""
function convergence_and_basins_of_attraction(mapper::AttractorMapper, grid; show_progress = true)
    if length(grid) != dimension(mapper.ds)
        @error "The mapper and the grid must have the same dimension"
    end
    basins = zeros(length.(grid))
    iterations = zeros(Int, length.(grid))
    I = CartesianIndices(basins)
    progress = ProgressMeter.Progress(
        length(basins); desc = "Basins and convergence: ", dt = 1.0
    )

    for (k, ind) in enumerate(I)
        show_progress && ProgressMeter.update!(progress, k)
        u0 = Attractors.generate_ic_on_grid(grid, ind)
        basins[ind] = mapper(u0)
        iterations[ind] = iterations_to_converge(mapper)
    end
    attractors = extract_attractors(mapper)
    return basins, attractors, iterations
end

"""
    convergence_and_basins_of_attraction(mapper::AttractorsViaRecurrences, grid)

An extension of [`basins_fractions`](@ref).
TODO.
"""
