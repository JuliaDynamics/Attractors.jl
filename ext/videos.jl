function Attractors.animate_attractors_continuation(
    ds::DynamicalSystem, attractors_cont, fractions_cont, prange, pidx; kw...)
    pcurve = [[pidx => p] for p in prange]
    return animate_attractors_continuation(ds, attractors_cont, fractions_cont, pcurve; prange, kw...)
end

function Attractors.animate_attractors_continuation(
        ds::DynamicalSystem, attractors_cont, fractions_cont, pcurve;
        savename = "attracont.mp4", access = SVector(1, 2),
        limits = auto_attractor_lims(attractors_cont, access),
        framerate = 4, markersize = 10,
        ukeys = unique_keys(attractors_cont),
        colors = colors_from_keys(ukeys),
        markers = markers_from_keys(ukeys),
        Δt = isdiscretetime(ds) ? 1 : 0.05,
        T = 100,
        figure = NamedTuple(), axis = NamedTuple(), fracaxis = NamedTuple(),
        legend = NamedTuple(),
        add_legend = length(ukeys) ≤ 6,
        # for the extra plot
        a2rs = nothing,
        prange = 1:length(attractors_cont),
        series_kwargs = NamedTuple(),
        a2rs_ratio = 0.5,
        a2rs_ylabels = ["" for i in 1:(isnothing(a2rs) ? 1 : length(a2rs))],
        parameter_name = "parameter",
        vline_kwargs = (linestyle = :dash, linewidth = 3, color = "black"),
    )
    length(access) ≠ 2 && error("Need two indices to select two dimensions of `ds`.")
    K = length(ukeys)
    fig = Figure(; figure...)
    ax = Axis(fig[1,:][1,1]; limits, axis...)
    fracax = Axis(fig[1,:][1,2]; width = 50, limits = (0,1,0,1), ylabel = "fractions",
        yaxisposition = :right, fracaxis...
    )
    hidedecorations!(fracax; label = false)

    # setup attractor axis (note we change colors to be transparent)
    att_obs = Dict(k => Observable(Point2f[]) for k in ukeys)
    plotf! = isdiscretetime(ds) ? scatter! : scatterlines!
    for k in ukeys
        plotf!(ax, att_obs[k]; color = (colors[k], 0.75), label = "$k", markersize, marker = markers[k])
    end
    if add_legend
        axislegend(ax; legend...)
    end

    # setup fractions axis
    heights = Observable(fill(0.1, K))
    barcolors = [colors[k] for k in ukeys]
    barplot!(fracax, fill(0.5, K), heights; width = 1, gap = 0, stack=1:K, color = barcolors)

    # add additional attractor series plots if chosen
    if !isnothing(a2rs)
        a2rs isa Vector || error("`a2rs` must be a `Vector`.")
        vline_obs = Observable(prange[1])
        n = length(a2rs)
        extra_axs = [Axis(fig[2,:][i, :]) for i in 1:n]
        extra_axs[1].title = "global continuation: features"
        linkxaxes!(extra_axs...)
        extra_axs[end].xlabel = string(parameter_name)
        rowsize!(fig.layout, 2, Relative(a2rs_ratio))
        for (j, (axa, a2r)) in enumerate(zip(extra_axs, a2rs))
            plot_attractors_curves!(axa, attractors_cont, a2r, prange;
                ukeys, colors, markers, add_legend = false, series_kwargs,
            )
            axa.ylabel = a2rs_ylabels[j]
            j < length(a2rs) && hidexdecorations!(axa; grid = false)
            # vline observable
            vlines!(axa, vline_obs; vline_kwargs...)
        end
        # layouting improvements
        content(fig[1,:]).alignmode = Mixed(left = 0, right = 0)
        content(fig[2,:]).alignmode = Mixed(left = 0, right = 0)
        Makie.update_state_before_display!(fig)
        yspace = maximum(tight_yticklabel_spacing!, extra_axs)
        for ax in extra_axs
            ax.yticklabelspace = yspace + 2
        end
        xlims!(extra_axs[end], nothing, nothing)
    end

    # okay, recording loop now...
    record(fig, savename, eachindex(attractors_cont); framerate) do i
        p = pcurve[i]
        text = sprint(show, p; context=:compact => true)
        ax.title = "global continuation: attractors\np: $text"
        attractors = attractors_cont[i]
        fractions = fractions_cont[i]
        set_parameters!(ds, p)
        heights[] = [get(fractions, k, 0) for k in ukeys]
        vline_obs[] = prange[i]
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

function auto_attractor_lims(attractors_cont, access)
    xmin = ymin = Inf
    xmax = ymax = -Inf
    for atts in attractors_cont
        for (k, A) in atts
            P = A[:, access]
            mini, maxi = minmaxima(P)
            xmin > mini[1] && (xmin = mini[1])
            ymin > mini[2] && (ymin = mini[2])
            xmax < maxi[1] && (xmax = maxi[1])
            ymax < maxi[2] && (ymax = maxi[2])
        end
    end
    dx = xmax - xmin
    dy = ymax - ymin
    return (xmin - 0.1dx, xmax + 0.1dx, ymin - 0.1dy, ymax + 0.1dy)
end

