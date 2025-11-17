##########################################################################################
# Auto colors/markers
##########################################################################################
using Random: shuffle!, Xoshiro
function colors_from_keys(ukeys)
    # Unfortunately, until `to_color` works with `Cycled`,
    # we need to explicitly add here some default colors...
    COLORS = [
        "#7143E0",
        "#191E44",
        "#0A9A84",
        "#AF9327",
        "#5F166D",
        "#6C768C",
    ]
    if length(ukeys) ≤ length(COLORS)
        colors = [COLORS[i] for i in eachindex(ukeys)]
    else # keep colorscheme, but add extra random colors
        n = length(ukeys) - length(COLORS)
        colors = shuffle!(Xoshiro(123), collect(cgrad(COLORS, n+1; categorical = true)))
        colors = append!(to_color.(COLORS), colors[1:(end-1)])
    end
    return Dict(k => colors[i] for (i, k) in enumerate(ukeys))
end

function markers_from_keys(ukeys)
    MARKERS = [:circle, :dtriangle, :rect, :star5, :xcross, :diamond,
    :hexagon, :cross, :pentagon, :ltriangle, :rtriangle, :hline, :vline, :star4,]
    markers = Dict(k => MARKERS[mod1(i, length(MARKERS))] for (i, k) in enumerate(ukeys))
    return markers
end

##########################################################################################
# Attractors
##########################################################################################
function Attractors.plot_attractors(a; access = SVector(1, 2), kw...)
    fig = Figure()
    AX = length(access) == 2 ? Axis : Axis3
    ax = AX(fig[1,1])
    plot_attractors!(ax, a; access, kw...)
    return fig
end

function Attractors.plot_attractors!(ax, attractors;
        ukeys = sort(collect(keys(attractors))), # internal argument just for other keywords
        colors = colors_from_keys(ukeys),
        markers = markers_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
        axislegend_kwargs = NamedTuple(),
        access = SVector(1, 2),
        sckwargs = (strokewidth = 0.5, strokecolor = :black,)
    )
    for k in ukeys
        k ∉ keys(attractors) && continue
        A = attractors[k]
        scatter!(ax, A[:, access];
            color = (colors[k], 0.9), markersize = 20,
            marker = markers[k],
            label = "$(labels[k])",
            sckwargs...
        )
    end
    add_legend && axislegend(ax; axislegend_kwargs...)
    return
end


##########################################################################################
# Basins
##########################################################################################
function Attractors.heatmap_basins_attractors(grid, basins::AbstractArray, attractors; kwargs...)
    if length(size(basins)) != 2
        error("Heatmaps only work in two dimensional basins!")
    end
    fig = Figure()
    ax = Axis(fig[1,1])
    heatmap_basins_attractors!(ax, grid, basins, attractors; kwargs...)
    return fig
end

function Attractors.heatmap_basins_attractors(BoA::ArrayBasinsOfAttraction; kwargs...) 
    return Attractors.heatmap_basins_attractors(BoA.grid, BoA.basins, BoA.attractors; kwargs...)
end

function Attractors.heatmap_basins_attractors!(ax, grid, basins, attractors;
        ukeys = unique(basins), # internal argument just for other keywords
        colors = colors_from_keys(ukeys),
        markers = markers_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
        access = SVector(1, 2),
        sckwargs = (strokewidth = 1.5, strokecolor = :white,)
    )

    sort!(ukeys) # necessary because colormap is ordered
    # Set up the (categorical) color map and colormap values
    cmap = cgrad([colors[k] for k in ukeys], length(ukeys); categorical = true)
    # Heatmap with appropriate colormap values. We need to transform
    # the basin array to one with values the sequential integers,
    # because the colormap itself also has as values the sequential integers
    ids = 1:length(ukeys)
    replace_dict = Dict(k => i for (i, k) in enumerate(ukeys))
    basins_to_plot = replace(basins, replace_dict...)
    heatmap!(ax, grid..., basins_to_plot;
        colormap = cmap,
        colorrange = (ids[1]-0.5, ids[end]+0.5),
    )
    # Scatter attractors
    plot_attractors!(ax, attractors;
        ukeys, colors, access, markers,
        labels, add_legend, sckwargs
    )
    return ax
end

function Attractors.heatmap_basins_attractors!(ax, BoA::ArrayBasinsOfAttraction;
        ukeys = unique(BoA.basins), # internal argument just for other keywords
        colors = colors_from_keys(ukeys),
        markers = markers_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
        access = SVector(1, 2),
        sckwargs = (strokewidth = 1.5, strokecolor = :white,)
    )
    return Attractors.heatmap_basins_attractors!(ax, BoA.grid, BoA.basins, BoA.attractors;
        ukeys = ukeys, # internal argument just for other keywords
        colors = colors,
        markers = markers,
        labels = labels,
        add_legend = add_legend,
        access = access,
        sckwargs = sckwargs
    )
end
##########################################################################################
# Shaded basins
##########################################################################################
function Attractors.shaded_basins_heatmap(grid, basins::AbstractArray, attractors, iterations;
    show_attractors = true,
    maxit = maximum(iterations),
    kwargs...)
    if length(size(basins)) != 2
        error("Heatmaps only work in two dimensional basins!")
    end
    fig = Figure()
    ax = Axis(fig[1,1]; kwargs...)
    shaded_basins_heatmap!(ax, grid, basins, iterations, attractors; maxit, show_attractors)
    return fig
end

function Attractors.shaded_basins_heatmap(BoA::ArrayBasinsOfAttraction, iterations;
    show_attractors = true,
    maxit = maximum(iterations),
    kwargs...
    )
    return Attractors.shaded_basins_heatmap(BoA.grid, BoA.basins, BoA.attractors, iterations;
    show_attractors = show_attractors, maxit = maxit, kwargs...)
end

function Attractors.shaded_basins_heatmap!(ax, grid, basins, iterations, attractors;
    ukeys = unique(basins),
    show_attractors = true,
    maxit = maximum(iterations))

    sort!(ukeys) # necessary because colormap is ordered
    ids = 1:length(ukeys)
    replace_dict = Dict(k => i for (i, k) in enumerate(ukeys))
    basins_to_plot = replace(basins.*1., replace_dict...)
    access = SVector(1,2)

    cmap, colors = custom_colormap_shaded(ukeys)
    markers = markers_from_keys(ukeys)
    labels = Dict(ukeys .=> ukeys)
    add_legend = length(ukeys) < 7

    it = findall(iterations .> maxit)
    iterations[it] .= maxit
    for i in ids
        ind = findall(basins_to_plot .== i)
        mn = minimum(iterations[ind])
        mx = maximum(iterations[ind])
        basins_to_plot[ind] .= basins_to_plot[ind] .+ 0.99.*(iterations[ind].-mn)/mx
    end

    # The colormap is constructed in such a way that the first color maps
    # from id to id+0.99, id is an integer describing the current basin.
    # Each id has a specific color associated and the gradient goes from
    # light color (value id) to dark color (value id+0.99). It produces
    # a shading proportional to a value associated to a specific pixel.
    heatmap!(ax, grid..., basins_to_plot;
        colormap = cmap,
        colorrange = (ids[1], ids[end]+1),
    )
    # Scatter attractors
    if show_attractors
        for (i, k) ∈ enumerate(ukeys)
            k ≤ 0 && continue
            A = attractors[k]
            x, y = columns(A[:, access])
            scatter!(ax, x, y;
                color = colors[k], markersize = 20,
                marker = markers[k],
                strokewidth = 1.5, strokecolor = :white,
                label = "$(labels[k])"
            )
        end
        # Add legend using colors only
        add_legend && axislegend(ax)
    end
    return ax
end

function Attractors.shaded_basins_heatmap!(ax, BoA::ArrayBasinsOfAttraction, iterations;
    ukeys = unique(basins), show_attractors = true, maxit = maximum(iterations))
    return Attractors.shaded_basins_heatmap!(ax, BoA.grid, BoA.basins, iterations, BoA.attractors;
        ukeys = ukeys, show_attractors = show_attractors, maxit = maxit)
end

function custom_colormap_shaded(ukeys)
    # Light and corresponding dark colors for shading of basins of attraction
    colors = colors_from_keys(ukeys)
    n = length(colors)
    # Note that the factor to define light and dark color
    # is arbitrary.
    LIGHT_COLORS = [darken_color(colors[k],0.3) for k in ukeys]
    DARK_COLORS = [darken_color(colors[k],1.7) for k in ukeys]
    v_col = Array{typeof(LIGHT_COLORS[1]),1}(undef,2*n)
    vals = zeros(2*n)
    for k in eachindex(ukeys)
        v_col[2*k-1] = LIGHT_COLORS[k]
        v_col[2*k] = DARK_COLORS[k]
        vals[2*k-1] = k-1
        vals[2*k] = k-1+0.9999999999
    end
    return cgrad(v_col, vals/maximum(vals)), colors
end

"""
    darken_color(c, f = 1.2)
Darken given color `c` by a factor `f`.
If `f` is less than 1, the color is lightened instead.
"""
function darken_color(c, f = 1.2)
    c = to_color(c)
    return RGBAf(clamp.((c.r/f, c.g/f, c.b/f, c.alpha), 0, 1)...)
end

##########################################################################################
# Continuation
##########################################################################################
function Attractors.plot_basins_curves(fractions_cont, args...; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    ax.xlabel = "parameter"
    ax.ylabel = "basins %"
    plot_basins_curves!(ax, fractions_cont, args...; kwargs...)
    return fig
end

function Attractors.plot_basins_curves!(ax, fractions_cont, prange = 1:length(fractions_cont);
        ukeys = unique_keys(fractions_cont), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        separatorwidth = 1, separatorcolor = "white",
        add_legend = length(ukeys) < 7,
        axislegend_kwargs = (position = :lt,),
        series_kwargs = NamedTuple(),
        markers = markers_from_keys(ukeys),
        style = :band,
        filler = NaN,
    )
    if !(prange isa AbstractVector{<:Real})
        error("!(prange <: AbstractVector{<:Real})")
    end
    bands = continuation_series(fractions_cont, style == :band ? 0.0 : filler, ukeys)
    if style == :band
        # transform to cumulative sum
        for j in 2:length(bands)
            bands[ukeys[j]] .+= bands[ukeys[j-1]]
        end
        for (j, k) in enumerate(ukeys)
            if j == 1
                l, u = 0, bands[k]
                l = fill(0f0, length(u))
            else
                l, u = bands[ukeys[j-1]], bands[ukeys[j]]
            end
            band!(ax, prange, l, u;
                color = colors[k], label = "$(labels[k])", series_kwargs...
            )
            if separatorwidth > 0 && j < length(ukeys)
                lines!(ax, prange, u; color = separatorcolor, linewidth = separatorwidth, linestyle = :solid)
            end
        end
        ylims!(ax, 0, 1)
    elseif style == :lines
        for k in ukeys
            scatterlines!(ax, prange, bands[k];
                color = colors[k], label = "$(labels[k])", marker = markers[k],
                markersize = 10, linewidth = 3, series_kwargs...
            )
        end
    else
        error("Incorrect style specification for basins fractions curves")
    end

    xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; axislegend_kwargs...)
    return
end

function Attractors.plot_attractors_curves(attractors_cont, attractor_to_real, prange = 1:length(attractors_cont); kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    ax.xlabel = "parameter"
    ax.ylabel = "attractors"
    plot_attractors_curves!(ax, attractors_cont, attractor_to_real, prange; kwargs...)
    return fig
end

function Attractors.plot_attractors_curves!(ax, attractors_cont, attractor_to_real, prange = 1:length(attractors_cont);
        kwargs...
    )
    # make the continuation info values and just propagate to the main function
    continuation_info = map(attractors_cont) do dict
        Dict(k => attractor_to_real(A) for (k, A) in dict)
    end
    plot_continuation_curves!(ax, continuation_info, prange; kwargs...)
end

function Attractors.plot_continuation_curves!(ax, continuation_info, prange = 1:length(continuation_info);
        ukeys = unique_keys(continuation_info), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
        markers = markers_from_keys(ukeys),
        series_kwargs = NamedTuple(),
        axislegend_kwargs = (position = :lt,)
    )

    series = continuation_series(continuation_info, NaN, ukeys)
    for (k, v) in series
        v[isinf.(v)] .= NaN
    end

    for k in ukeys
        scatterlines!(ax, prange, series[k];
            color = colors[k], label = "$(labels[k])", marker = markers[k],
            markersize = 10, linewidth = 3, series_kwargs...
        )
    end
    xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; axislegend_kwargs..., unique = true)
    return
end

function Attractors.plot_continuation_curves(args...; kw...)
    fig = Figure()
    ax = Axis(fig[1,1])
    plot_continuation_curves!(ax, args...; kw...)
    return fig
end

# Mixed: basins and attractors
function Attractors.plot_basins_attractors_curves(fractions_cont, attractors_cont, a2r::Function, prange = 1:length(attractors_cont);
        kwargs...
    )
    return Attractors.plot_basins_attractors_curves(fractions_cont, attractors_cont, [a2r], prange; kwargs...)
end

# Special case with multiple attractor projections:
function Attractors.plot_basins_attractors_curves(
        fractions_cont, attractors_cont,
        a2rs::Vector, prange = 1:length(attractors_cont);
        ukeys = unique_keys(fractions_cont), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        markers = markers_from_keys(ukeys),
        style = :band,
        kwargs...
    )
    # generate figure and axes; add labels and stuff
    fig = Figure()
    axb = Axis(fig[1,1])
    A = length(a2rs)
    axs = [Axis(fig[1+i, 1]; ylabel = "attractors_$(i)") for i in 1:A]
    linkxaxes!(axb, axs...)
    if A == 1 # if we have only 1, make it pretier
        axs[1].ylabel = "attractors"
    end
    axs[end].xlabel = "parameter"
    axb.ylabel = "basins %"
    hidexdecorations!(axb; grid = false)
    for i in 1:A-1
        hidexdecorations!(axs[i]; grid = false)
    end
    # plot basins and attractors
    plot_basins_curves!(axb, fractions_cont, prange; ukeys, colors, labels, style, kwargs...)
    for (axa, a2r) in zip(axs, a2rs)
        plot_attractors_curves!(axa, attractors_cont, a2r, prange;
            ukeys, colors, markers, kwargs..., add_legend = false, # can be true for fractions
        )
    end
    return fig
end

# This function is kept for backwards compatibility only, really.
function Attractors.plot_basins_attractors_curves!(axb, axa, fractions_cont, attractors_cont,
        attractor_to_real, prange = 1:length(attractors_cont);
        ukeys = unique_keys(fractions_cont), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        kwargs...
    )

    if length(fractions_cont) ≠ length(attractors_cont)
        error("fractions and attractors don't have the same amount of entries")
    end

    plot_basins_curves!(axb, fractions_cont, prange; ukeys, colors, labels, kwargs...)

    plot_attractors_curves!(axa, attractors_cont, attractor_to_real, prange;
        ukeys, colors, add_legend = false, # coz its true for fractions
    )
    return
end
