@doc raw"""
  __integral_phi__ = Dict{Float64,Dict{Float64,uwreal}}();

  Cache variable for the integral inside `phi(g,bcoef::Vector)`.
  To allow independent computation with different loop approximation
  and/or different beta coefficient, the  structure of `__integral_phi__` 
  is `__integral_phi__[last beta coefficient][g]`
"""
__integral_phi__ = Dict{Float64,Dict{Float64,uwreal}}();

@doc raw"""
  reset_cache()

  Reset the cached variable. 
"""
function reset_cache(; key::Union{Nothing,Float64})
  if isnothing(key)
    global __integral_phi__ = nothing
    global __integral_phi__ = Dict{Float64,Dict{Float64,uwreal}}();
  else
    global __integral_phi__[key] = nothing
    global __integral_phi__[key] = Dict{Float64,Dict{Float64,uwreal}}();

  end
end