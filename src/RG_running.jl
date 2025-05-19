module RG_running

using ADerrors, LaTeXStrings, ForwardDiff, Roots, QuadGK

import QuadGK: kronrod

"""An issue between QuadGK (#122) and ForwardDiff, probaly linked to ForwardDiff-issue #624 make this overloading necessary.
If and when the issue is solved, this line should be removed"""

QuadGK.kronrod(::Type{<:ForwardDiff.Dual{T,V,N}}, n::Integer) where {T,V,N} = QuadGK.kronrod(V, n)
#QuadGK.cachedrule(::Type{<:ForwardDiff.Dual{<:Any, T}}, n::Integer) where {T<:Number} =
#    QuadGK._cachedrule(typeof(float(real(one(T)))), Int(n))

include("cache.jl")
include("parameters.jl")
include("alpha_running.jl")

export reset_cache
export HyperPar,beta_function_coeff,MSbar,MSbar_float
export alpha


end #end module RG_running
