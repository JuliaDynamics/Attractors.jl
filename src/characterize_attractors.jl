
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