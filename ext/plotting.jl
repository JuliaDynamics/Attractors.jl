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
        colors = shuffle!(Xoshiro(123), collect(cgrad(:darktest, n+1; categorical = true)))
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


function Attractors.heatmap_basins_attractors!(ax, grid, basins, attractors;
        ukeys = unique(basins), # internal argument just for other keywords
        colors = colors_from_keys(ukeys),
        markers = markers_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
        access = SVector(1, 2)
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
    for (i, k) ∈ enumerate(ukeys)
        k ≤ 0 && continue
        A = attractors[k]
        x, y = columns(A[:, access])
        scatter!(ax, x, y;
            color = colors[k], markersize = 20,
            marker = markers[k],
            strokewidth = 1.5, strokecolor = :white,
            label = "$(labels[k])",
        )
    end
    # Add legend using colors only
    add_legend && axislegend(ax)
    return ax
end

# function colors_from_keys_shaded(ukeys)
#     # Unfortunately, until `to_color` works with `Cycled`,
#     # we need to explicitly add here some default colors...
#     COLORS =[ 
#         :black,
#         :red,
#         :green, 
#         :blue,
#         :purple,
#         :yellow,]
#     n = length(COLORS)
#     v_col = []
#     vals = zeros(2*length(ukeys))
#     for k in eachindex(ukeys)
#         push!(v_col, COLORS[mod(k-1,n)+1])
#     end

#     return Dict(k => v_col[i] for (i, k) in enumerate(ukeys))
# end

"""
    darken_color(c, f = 1.2)
Darken given color `c` by a factor `f`.
If `f` is less than 1, the color is lightened instead.
"""
function darken_color(c, f = 1.2)
    c = to_color(c)
    return RGBAf(clamp.((c.r/f, c.g/f, c.b/f, c.alpha), 0, 1)...)
end

function custom_colormap_shaded(ukeys)
    # Light and corresponding dark colors for shading of basins of attraction
    colors = colors_from_keys(ukeys)
    n = length(colors)
    LIGHT_COLORS = [darken_color(colors[k],0.5) for k in 1:n]
    DARK_COLORS = [darken_color(colors[k],1.5) for k in 1:n]
    v_col = Array{typeof(LIGHT_COLORS[1]),1}(undef,2*n)
    vals = zeros(2*n)
    for k in eachindex(ukeys)
        v_col[2*k-1] = LIGHT_COLORS[k]
        v_col[2*k] = DARK_COLORS[k]
        vals[2*k-1] = k-1
        vals[2*k] = k-1+0.9999999999
    end
    return cgrad(v_col, vals/maximum(vals))
end

function Attractors.shaded_basins_heatmap(grid, basins::AbstractArray, iterations, attractors; 
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

function Attractors.shaded_basins_heatmap!(ax, grid, basins, iterations, attractors; 
    ukeys = unique(basins), 
    show_attractors = true, 
    maxit = maximum(iterations)) 

    sort!(ukeys) # necessary because colormap is ordered
    ids = 1:length(ukeys)
    replace_dict = Dict(k => i for (i, k) in enumerate(ukeys))
    basins_to_plot = replace(basins.*1., replace_dict...)
    access = SVector(1,2)
    
    cmap = custom_colormap_shaded(ukeys)
    colors = colors_from_keys(ukeys)
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


##########################################################################################
# Continuation
##########################################################################################
function Attractors.plot_basins_curves(fractions_curves, prange; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    ax.xlabel = "parameter"
    ax.ylabel = "basins %"
    plot_basins_curves!(ax, fractions_curves, prange; kwargs...)
    return fig
end

function Attractors.plot_basins_curves!(ax, fractions_curves, prange = 1:length(fractions_curves);
        ukeys = unique_keys(fractions_curves), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        separatorwidth = 1, separatorcolor = "white",
        add_legend = length(ukeys) < 7,
        axislegend_kwargs = (position = :lt,),
        series_kwargs = (markersize = 5, linewidth = 3),
        markers = markers_from_keys(ukeys),
        style = :band,
    )
    if !(prange isa AbstractVector{<:Real})
        error("!(prange <: AbstractVector{<:Real})")
    end
    bands = fractions_series(fractions_curves, ukeys)
    if style == :band
        # transform to cumulative sum
        for j in 2:length(bands)
            bands[j] .+= bands[j-1]
        end
        for (j, k) in enumerate(ukeys)
            if j == 1
                l, u = 0, bands[j]
                l = fill(0f0, length(u))
            else
                l, u = bands[j-1], bands[j]
            end
            band!(ax, prange, l, u;
                color = colors[k], label = "$(labels[k])", linewidth = 4,
                series_kwargs...
            )
            if separatorwidth > 0 && j < length(ukeys)
                lines!(ax, prange, u; color = separatorcolor, linewidth = separatorwidth)
            end
        end
        ylims!(ax, 0, 1)
    elseif style == :lines
        for (j, k) in enumerate(ukeys)
            scatterlines!(ax, prange, bands[j];
                color = colors[k], label = "$(labels[k])", marker = markers[k],
                series_kwargs...
            )
        end
    else
        error()
    end

    xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; axislegend_kwargs...)
    return
end

function fractions_series(fractions_curves, ukeys = unique_keys(fractions_curves))
    bands = [zeros(length(fractions_curves)) for _ in ukeys]
    for i in eachindex(fractions_curves)
        for (j, k) in enumerate(ukeys)
            bands[j][i] = get(fractions_curves[i], k, 0)
        end
    end
    return bands
end

function Attractors.plot_attractors_curves(attractors_info, attractor_to_real, prange = 1:length(attractors_info); kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    plot_attractors_curves!(ax, attractors_info, attractor_to_real, prange; kwargs...)
    return fig
end

function Attractors.plot_attractors_curves!(ax, attractors_info, attractor_to_real, prange = 1:length(attractors_info);
        ukeys = unique_keys(attractors_info), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        add_legend = length(ukeys) < 7,
        markers = markers_from_keys(ukeys),
        axislegend_kwargs = (position = :lt,)
    )
    for i in eachindex(attractors_info)
        attractors = attractors_info[i]
        for (k, A) in attractors
            val = attractor_to_real(A)
            scatter!(ax, prange[i], val;
                color = colors[k], markers = markers[k], label = string(labels[k]),
            )
        end
    end
    xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; axislegend_kwargs..., unique = true)
    return
end

function Attractors.plot_basins_attractors_curves(fractions_curves, attractors_info, attractor_to_real, prange = 1:length(attractors_info);
        kwargs...
    )
    fig = Figure()
    axb = Axis(fig[1,1])
    axa = Axis(fig[2,1])
    axa.xlabel = "parameter"
    axa.ylabel = "attractors"
    axb.ylabel = "basins %"
    hidexdecorations!(axb; grid = false)

    Attractors.plot_basins_attractors_curves!(axb, axa, fractions_curves, attractors_info,
        attractor_to_real, prange; kwargs...,
    )
    return fig
end

function Attractors.plot_basins_attractors_curves!(axb, axa, fractions_curves, attractors_info,
        attractor_to_real, prange = 1:length(attractors_info);
        ukeys = unique_keys(fractions_curves), # internal argument
        colors = colors_from_keys(ukeys),
        labels = Dict(ukeys .=> ukeys),
        kwargs...
    )

    if length(fractions_curves) ≠ length(attractors_info)
        error("fractions and attractors don't have the same amount of entries")
    end

    plot_basins_curves!(axb, fractions_curves, prange; ukeys, colors, labels, kwargs...)

    plot_attractors_curves!(axa, attractors_info, attractor_to_real, prange;
        ukeys, colors, add_legend = false, # coz its true for fractions
    )
    return
end

##########################################################################################
# Videos
##########################################################################################
function Attractors.animate_attractors_continuation(
        ds::DynamicalSystem, attractors_info, fractions_curves, prange, pidx;
        savename = "test.mp4", access = SVector(1, 2),
        limits = (0,1,0,1),
        framerate = 4, markersize = 10,
        ukeys = unique_keys(attractors_info),
        colors = colors_from_keys(ukeys),
        markers = markers_from_keys(ukeys),
        Δt = isdiscretetime(ds) ? 1 : 0.05,
        T = 100,
    )
    length(access) ≠ 2 && error("Need two indices to select two dimensions of `ds`.")
    K = length(ukeys)
    fig = Figure()
    ax = Axis(fig[1,1]; limits)
    fracax = Axis(fig[1,2]; width = 50, limits = (0,1,0,1), ylabel = "fractions", yaxisposition = :right)
    hidedecorations!(fracax)
    fracax.ylabelvisible = true

    # setup attractor axis (note we change colors to be transparent)
    att_obs = Dict(k => Observable(Point2f[]) for k in ukeys)
    plotf! = isdiscretetime(ds) ? scatter! : scatterlines!
    for k in ukeys
        plotf!(ax, att_obs[k]; color = (colors[k], 0.75), label = "$k", markersize, marker = markers[k])
    end
    axislegend(ax)

    # setup fractions axis
    heights = Observable(fill(0.1, K))
    barcolors = [colors[k] for k in ukeys]
    barplot!(fracax, fill(0.5, K), heights; width = 1, gap = 0, stack=1:K, color = barcolors)

    record(fig, savename, eachindex(prange); framerate) do i
        p = prange[i]
        ax.title = "p = $p"
        attractors = attractors_info[i]
        fractions = fractions_curves[i]
        set_parameter!(ds, pidx, p)
        heights[] = [get(fractions, k, 0) for k in ukeys]

        for (k, att) in attractors
            tr, tvec = trajectory(ds, T, first(vec(att)); Δt)
            att_obs[k][] = vec(tr[:, access])
            notify(att_obs[k])
        end
        # also ensure that attractors that don't exist are cleared
        for k in setdiff(ukeys, collect(keys(attractors)))
            att_obs[k][] = Point2f[]; notify(att_obs[k])
        end
    end
    return fig
end


