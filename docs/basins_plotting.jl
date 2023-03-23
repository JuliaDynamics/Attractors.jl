# Functions dedicated into plotting basins related stuff
# Some default colors are needed...
COLORS = [
    "#7143E0",
    "#191E44",
    "#0A9A84",
    "#C0A12B",
    "#701B80",
    "#2E6137",
]


function animate_attractors_continuation(
        ds, attractors_info, prange, pidx;
        savename = "test.mp4", access = [1,2],
        limits = (-1,3,-2,2),
        framerate = 4, markersize = 10
    )
    ukeys = unique_keys(attractors_info)
    K = length(ukeys)
    fig = Figure()
    ax = Axis(fig[1,1]; limits)
    fracax = Axis(fig[1,2]; width = 50, limits = (0,1,0,1))
    hidedecorations!(fracax)

    colors = Dict(k => (to_color(COLORS[i]), 0.75) for (i, k) in enumerate(ukeys))
    att_obs = Dict(k => Observable(Point2f[]) for k in ukeys)
    for k in ukeys
        scatter!(ax, att_obs[k]; color = colors[k],
        label = "$k", markersize = markersize + rand(-4:4))
    end
    axislegend(ax)

    # setup fractions axis
    heights = Observable(fill(0.1, K))
    colors = [to_color(COLORS[i]) for i in ukeys]
    barplot!(fracax, fill(0.5, K), heights; width = 1, gap = 0, stack=1:K, color = colors)
    display(fig)

    record(fig, savename, eachindex(prange); framerate) do i
        p = prange[i]
        ax.title = "p = $p"
        attractors = attractors_info[i]
        fractions = fractions_curves[i]
        set_parameter!(ds, pidx, p)
        heights[] = [get(fractions, k, 0) for k in ukeys]

        for (k, att) in attractors
            tr = trajectory(ds, 1000, rand(vec(att)); Δt = 1)
            att_obs[k][] = vec(tr[:, access])
            notify(att_obs[k])
        end
        # also ensure that attractors that don't exist are cleared
        for k in setdiff(ukeys, collect(keys(attractors)))
            att_obs[k][] = Point2f[]; notify(att_obs[k])
        end
    end

end

function plot_attractors(attractors::Dict;  access = [1,2], markersize = 12)
    fig = Figure()
    ax = Axis(fig[1,1])
    ukeys = keys(attractors)
    colors = Dict(k => (to_color(COLORS[i]), 0.75) for (i, k) in enumerate(ukeys))
    for k in ukeys
        scatter!(ax, vec(attractors[k][:, access]); color = colors[k],
        label = "$k", markersize = markersize + rand(-2:4))
    end
    axislegend(ax)
    return fig
end

function basins_curves_plot(fractions_curves, prange; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    basins_curves_plot!(ax, fractions_curves, prange; kwargs...)
    return fig
end

"""
    basins_curves_plot!(ax::Axis, fractions_curves, prange; kwargs...)

Plot the fractions of basins of attraction versus a parameter range,
i.e., visualize the output of [`continuation`](@ref).
Keywords:
```julia
labels = Dict(ukeys .=> ukeys),
colors = colors_from_keys(ukeys),
separatorwidth = 1,
separatorcolor = "white",
add_legend = length(ukeys) < 8,
```
"""
function basins_curves_plot!(ax, fractions_curves, prange = 1:length(fractions_curves);
        ukeys = unique_keys(fractions_curves), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        separatorwidth = 1, separatorcolor = "white",
        add_legend = length(ukeys) < 7,
    )
    if !(prange isa AbstractVector{<:Real})
        error("!(prange <: AbstractVector{<:Real})")
    end
    bands = fractions_to_cumulative(fractions_curves, prange, ukeys)
    for (j, k) in enumerate(ukeys)
        if j == 1
            l, u = 0, bands[j]
            l = fill(0f0, length(u))
        else
            l, u = bands[j-1], bands[j]
        end
        band!(ax, prange, l, u; color = colors[k], label = "$(labels[k])", linewidth = 4)
        if separatorwidth > 0 && j < length(ukeys)
            lines!(ax, prange, u; color = separatorcolor, linewidth = separatorwidth)
        end
    end
    ylims!(ax, 0, 1); xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; position = :lt)
    return
end

using Random: shuffle!, Xoshiro
function colors_from_keys(ukeys)
    if length(ukeys) ≤ length(COLORS)
        colors = [COLORS[i] for i in eachindex(ukeys)]
    else # keep colorscheme, but add extra random colors
        n = length(ukeys) - length(COLORS)
        colors = shuffle!(Xoshiro(123), collect(cgrad(:darktest, n+1; categorical = true)))
        colors = append!(to_color.(COLORS), colors[1:(end-1)])
    end
    return Dict(k => colors[i] for (i, k) in enumerate(ukeys))
end

function fractions_to_cumulative(fractions_curves, prange, ukeys = unique_keys(fractions_curves))
    bands = [zeros(length(prange)) for _ in ukeys]
    for i in eachindex(fractions_curves)
        for (j, k) in enumerate(ukeys)
            bands[j][i] = get(fractions_curves[i], k, 0)
        end
    end
    # transform to cumulative sum
    for j in 2:length(bands)
        bands[j] .+= bands[j-1]
    end
    return bands
end

function attractors_curves_plot(attractors_info, attractor_to_real, prange = 1:length(attractors_info); kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    attractors_curves_plot!(ax, attractors_info, attractor_to_real, prange; kwargs...)
    return fig
end

function attractors_curves_plot!(ax, attractors_info, attractor_to_real, prange = 1:length(attractors_info);
        ukeys = unique_keys(attractors_info), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
    )
    for i in eachindex(attractors_info)
        attractors = attractors_info[i]
        for (k, A) in attractors
            val = attractor_to_real(A)
            scatter!(ax, prange[i], val; color = colors[k])
        end
    end
    xlims!(ax, minimum(prange), maximum(prange))
    if add_legend
        labels = map(k -> labels[k], ukeys)
        ele = map(k -> MarkerElement(; marker = :circle, color = colors[k]), ukeys)
        # From the source code of `axislegend`:
        Legend(ax.parent, ele, labels;
            bbox = ax.scene.px_area,
            Makie.legend_position_to_aligns(:lt)...,
            margin = (10, 10, 10, 10),
        )
    end
    return
end

function attractor_type(A)
    return "len=$(length(A))"
    # if length(A) == 1
    #     l =  "fixed p."
    # else
    #     # fractal dimension
    #     D = grassberger_proccacia_dim(A)
    #     @show D
    #     if D < 0.1
    #         l =  "fixed p."
    #     elseif D < 1
    #         l =  "limit c."
    #     else
    #         l =  "chaotic"
    #     end
    # end
    # return "$(l)"
end

# TODO: This function needs to be improved somehow
function create_attractor_names(ukeys, attractors_info, attractor_name)
    map(ukeys) do k
        # find first attractor with this key
        di = findlast(d -> haskey(d, k), attractors_info)
        A = attractors_info[di][k]
        label = attractor_name(A)
        return "$(k): $(label)"
    end
end


function basins_attractors_curves_plot(fractions_curves, attractors_info, attractor_to_real, prange = 1:length(attractors_info);
        ukeys = unique_keys(fractions_curves), # internal argument
        colors = colors_from_keys(ukeys),
    )
    fig = Figure()
    axb = Axis(fig[1,1])
    axa = Axis(fig[2,1])
    basins_curves_plot!(axb, fractions_curves, prange; ukeys, colors)
    # Okay here we need to delete key `-1` if it exists:
    if -1 ∈ ukeys
        popfirst!(ukeys)
        delete!(colors, -1)
    end
    attractors_curves_plot!(axa, attractors_info, attractor_to_real, prange;
        ukeys, colors, add_legend = false, # coz its true for fractions
    )
    hidexdecorations!(axb, ticks = false, grid = false)
    axa.xlabel = "parameter"
    axa.ylabel = "attractors"
    axb.ylabel = "basins %"
    return fig
end
