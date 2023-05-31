
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
