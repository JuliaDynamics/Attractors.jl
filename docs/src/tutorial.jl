# # Attractors.jl Tutorial

# [`Attractors`](@ref) is a module of the **DynamicalSystems.jl** library.
# This tutorial will walk you through its main functionality.
# That is, given a `DynamicalSystem` instance, find all its attractors and then
# continue these attractors, and their stability properties, across a parameter value.
# It also offers various functions that compute nonlocal stability properties for an
# attractor, any of which can be used in the continuation to quantify stability.

# Besides this, there are some other stuff, like for example [`edgestate`](@ref),
# but we won't cover anything else in this introductory tutorial. See the [examples](@ref examples)
# page instead.

# ## Outline




using Attractors

# ```@docs; canonical=false
# Attractors
# ```