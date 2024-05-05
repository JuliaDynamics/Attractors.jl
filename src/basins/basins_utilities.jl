"""
    basins_fractions(basins::AbstractArray [,ids]) → fs::Dict

Calculate the state space fraction of the basins of attraction encoded in `basins`.
The elements of `basins` are integers, enumerating the attractor that the entry of
`basins` converges to (i.e., like the output of [`basins_of_attraction`](@ref)).
Return a dictionary that maps attractor IDs to their relative fractions.
Optionally you may give a vector of `ids` to calculate the fractions of only
the chosen ids (by default `ids = unique(basins)`).

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
Return `basins, attractors, convergence`, with `basins, attractors` as in
`basins_of_attraction`, and `convergence` being an array with same shape
as `basins`. It contains the time each initial condition took
to converge to its attractor.
It is useful to give to [`shaded_basins_heatmap`](@ref).

See also [`convergence_time`](@ref).

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
        iterations[ind] = convergence_time(mapper)
    end
    attractors = extract_attractors(mapper)
    return basins, attractors, iterations
end

"""
    convergence_and_basins_fractions(mapper::AttractorMapper, ics::StateSpaceSet)

An extension of [`basins_fractions`](@ref).
Return `fs, labels, convergence`. The first two are as in `basins_fractions`,
and `convergence` is a vector containing the time each initial condition took
to converge to its attractor.
Only usable with mappers that support `id = mapper(u0)`.

See also [`convergence_time`](@ref).

# Keyword arguments

- `show_progress = true`: show progress bar.
"""
function convergence_and_basins_fractions(mapper::AttractorMapper, ics::AbstractStateSpaceSet;
        show_progress = true,
    )
    N = size(ics, 1)
    progress = ProgressMeter.Progress(N;
        desc="Mapping initial conditions to attractors:", enabled = show_progress
    )
    fs = Dict{Int, Int}()
    labels = Vector{Int}(undef, N)
    iterations = Vector{typeof(current_time(mapper.ds))}(undef, N)

    for i ∈ 1:N
        ic = _get_ic(ics, i)
        label = mapper(ic; show_progress)
        fs[label] = get(fs, label, 0) + 1
        labels[i] = label
        iterations[i] = convergence_time(mapper)
        show_progress && ProgressMeter.next!(progress)
    end
    # Transform count into fraction
    ffs = Dict(k => v/N for (k, v) in fs)
    return ffs, labels, iterations
end
