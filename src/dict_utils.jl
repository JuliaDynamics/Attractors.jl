export unique_keys, swap_dict_keys!, next_free_id

# Utility functions for managing dictionary keys that are useful
# in continuation and attractor matching business
# Thanks a lot to Valentin (@Seelengrab) for generous help in the key swapping code.
# Swapping keys and ensuring that everything works "as expected" is,
# surprising, one of the hardest things to code in Attractors.jl.

"""
    swap_dict_keys!(d::Dict, matching_map::Dict)

Swap the keys of a dictionary `d` given a `matching_map`
which maps old keys to new keys. Also ensure that a swap can happen at most once,
e.g., if input `d` has a key `4`, and `rmap = Dict(4 => 3, 3 => 2)`,
then the key `4` will be transformed to `3` and not further to `2`.
"""
function swap_dict_keys!(fs::Dict, _rmap::AbstractDict)
    isempty(_rmap) && return
    # Transform rmap so it is sorted in decreasing order,
    # so that double swapping can't happen
    rmap = sort!(collect(_rmap); by = x -> x[2])
    cache = Tuple{keytype(fs), valtype(fs)}[]
    for (oldkey, newkey) in rmap
        haskey(fs, oldkey) || continue
        oldkey == newkey && continue
        tmp = pop!(fs, oldkey)
        if !haskey(fs, newkey)
            fs[newkey] = tmp
        else
            push!(cache, (newkey, tmp))
        end
    end
    for (k, v) in cache
        fs[k] = v
    end
    return fs
end

"""
    overwrite_dict!(old::Dict, new::Dict)

In-place overwrite the `old` dictionary for the key-value pairs of the `new`.
"""
function overwrite_dict!(old::Dict, new::Dict)
    empty!(old)
    for (k, v) in new
        old[k] = v
    end
end

"""
    additive_dict_merge!(d1::Dict, d2::Dict)

Merge keys and values of `d2` into `d1` additively: the values of the same keys
are added together in `d1` and new keys are given to `d1` as-is.
"""
function additive_dict_merge!(d1::Dict, d2::Dict)
    z = zero(valtype(d1))
    for (k, v) in d2
        d1[k] = get(d1, k, z) + v
    end
    return d1
end


"""
    retract_keys_to_consecutive(v::Vector{<:Dict}) → rmap

Given a vector of dictionaries with various positive integer keys, retract all keys so that
consecutive integers are used. So if the dictionaries have overall keys 2, 3, 42,
then they will transformed to 1, 2, 3.

Return the replacement map used to replace keys in all dictionaries with
[`swap_dict_keys!`](@ref).

As this function is used in attractor matching in [`global_continuation`](@ref)
it skips the special key `-1`.
"""
function retract_keys_to_consecutive(v::Vector{<:Dict})
    ukeys = unique_keys(v)
    ukeys = setdiff(ukeys, [-1]) # skip key -1 if it exists
    rmap = Dict(k => i for (i, k) in enumerate(ukeys) if i != k)
    return rmap
end

"""
    unique_keys(v::Iterator{<:AbstractDict})
Given a vector of dictionaries, return a sorted vector of the unique keys
that are present across all dictionaries.
"""
function unique_keys(v)
    unique_keys = Set(keytype(first(v))[])
    for d in v
        for k in keys(d)
            push!(unique_keys, k)
        end
    end
    return sort!(collect(unique_keys))
end


"""
    next_free_id(new::Dict, old::Dict)

Return the minimum key of the "new" dictionary
that doesn't exist in the "old" dictionary.
If one of the two dictionaries are empty, return its maximum key + 1.
If both are empty, return 1.

The function assumes tha the dictionary keys are integers.
"""
next_free_id(a₊::AbstractDict, a₋::AbstractDict) = next_free_id(keys(a₊), keys(a₋))
function next_free_id(keys₊, keys₋)
    if length(keys₋) == 0 && length(keys(keys₊)) == 0
        return 1
    elseif length(keys₋) == 0
        return maximum(keys₊) + 1
    elseif length(keys₊) == 0
        return maximum(keys₋) + 1
    else
        s = setdiff(keys₊, keys₋)
        return isempty(s) ? maximum(keys₋) + 1 : minimum(s)
    end
end
