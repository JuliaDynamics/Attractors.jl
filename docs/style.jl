# This file is included in `index.md`

# %% Color theme definitions
mutable struct CyclicContainer
    c::Vector
    n::Int
end
CyclicContainer(c) = CyclicContainer(c, 0)

Base.length(c::CyclicContainer) = length(c.c)
Base.iterate(c::CyclicContainer, state=1) = Base.iterate(c.c, state)
Base.getindex(c::CyclicContainer, i) = c.c[mod1(i, length(c))]
Base.getindex(c::CyclicContainer, i::AbstractRange) = c.c[i]

function Base.getindex(c::CyclicContainer)
    c.n += 1
    c[c.n]
end

COLORSCHEME = [
    "#7143E0",
    "#191E44",
    "#0A9A84",
    "#C0A12B",
    "#701B80",
    "#2E6137",
]

COLORS = CyclicContainer(COLORSCHEME)
LINESTYLES = CyclicContainer(["-", ":", "--", "-."])

# %% Makie styling
# other styling elements for Makie
set_theme!(;
    palette = (color = COLORSCHEME,),
    fontsize = 22,
    figure_padding = 4,
    resolution = (1000, 500),
    linewidth = 3.0,
)

# %% Plotting default functions for basin stuff
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

function basins_fractions_plot!(ax, fractions_curves, prange;
        ukeys = unique_keys(fractions_curves), # internal argument
        labels = Dict(ukeys .=> ukeys),
        colors = colors_from_keys(ukeys),
        separatorwidth = 0, separatorcolor = "white",
        add_legend = length(ukeys) < 8,
    )
    bands = fractions_to_cumulative(fractions_curves, prange, ukeys)
    for (j, k) in enumerate(ukeys)
        if j == 1
            l, u = 0, bands[j]
        else
            l, u = bands[j-1], bands[j]
        end
        band!(ax, prange, l, u; color = colors[j], label = "$(labels[k])", linewidth = 4)
        if separatorwidth > 0 && j < length(ukeys)
            lines!(ax, prange, u; color = separatorcolor, linewidth = separatorwidth)
        end
    end
    ylims!(ax, 0, 1); xlims!(ax, minimum(prange), maximum(prange))
    add_legend && axislegend(ax; position = :lt)
    return
end

using Random: Xoshiro
function colors_from_keys(ukeys)
    if length(ukeys) ≤ length(COLORS)
        COLORS
    else
        shuffle!(Xoshiro(123), collect(cgrad(:darktest, length(ukeys); categorical = true)))
    end
end

function fractions_to_cumulative(fractions_curves, prange, ukeys = unique_keys(fractions_curves))
    @show ukeys
    bands = [zeros(length(prange)) for _ in ukeys]
    for i in eachindex(fractions_curves)
        for (j, k) in enumerate(ukeys)
            @show i, j, k
            bands[j][i] = get(fractions_curves[i], k, 0)
        end
    end
    # transform to cumulative sum
    for j in 2:length(bands)
        bands[j] .+= bands[j-1]
    end
    return bands
end


function basins_fractions_plot(fractions_curves, prange; kwargs...)
    fig = Figure()
    ax = Axis(fig[1,1])
    basins_fractions_plot!(ax, fractions_curves, prange; kwargs...)
    return fig
end