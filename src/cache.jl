@doc raw"""
  __integral_phi__ = Dict{Float64,Dict{Float64,Any}}();

  Cache variable for the integral inside `phi(g,bcoef::Vector)`.
  To allow independent computation with different loop approximation
  and/or different beta coefficient, the  structure of `__integral_phi__` 
  is `__integral_phi__[last beta coefficient][g]`
"""
__integral_phi__ = Dict{Float64,Dict{Float64,Any}}();

@doc raw"""
  reset_cache(; key::Union{Nothing,Float64} = nothing)

  Reset the cached variable. 
"""
function reset_cache(; key::Union{Nothing,Float64}=nothing)
  if isnothing(key)
    empty!(__integral_phi__) 
  else
    delete!(__integral_phi__,key)
  end
end
