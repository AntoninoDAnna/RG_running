@doc raw"""  
    phi(g,bcoef::Vector)

It computes the function \phi(y=g_s^2). For numerical stability it was decided to compute it directly in terms of y = g_s^2. 
For rapidity, the results of the integrals are cached. When `phi(g_s^2,bcoef)` is called, it first check if the integral was already computed. 
If so it retrieve the cached results and return `phi`, otherwise it looks for the highest `g^2<g_s^2` for which the integral was computed, retrieve 
`cached_int[g^2]` and compute the integral from g^2 to g_s^2. The sum of these is cached such that each entry is the integral from `0` to the respective `g^2`.
  
WARNING: `phi` always return a Float64, and `bcoef` must not be an uwreal. g must be positive otherwise an `ArgumentError` is thrown
"""
function phi(g,bcoef::Vector)
  key = bcoef[end]
  if !(key in keys(__integral_phi__))
    __integral_phi__[key] = Dict{Float64,uwreal}(0.0=>0.0) ## make sure that there is at least one entry
  end
  if g<0
    throw(ArgumentError("Negative g in phi_y"))
  end
  gsqr= g^2
  #y = x^2
  cached_int = __integral_phi__[key]
  aux1(y,p) = sum([y^(i-1)*p[i] for i = eachindex(p)])  ## beta = x^3 * aux1
  aux2(y,p) = sum([y^(i-3)*p[i] for i = 3:lastindex(p)])
  integrand(y,p) =((p[1]-p[2]*y)*aux2(y,p) - p[2]^2)/(p[1]^2*aux1(y,p))

  gsqr_used = [k for (k,_) in cached_int]
  if !(gsqr in gsqr_used)
    aux = findall(x->x<gsqr,gsqr_used)
    lbound = maximum(gsqr_used[aux]) ## since gsqr is >0, worst case scenario is that lbound is 0.0
    cached_int[gsqr] = (cached_int[lbound] + int_error(integrand,lbound,gsqr ,uwreal.(bcoef)))
  end
  N = (bcoef[1] * gsqr)^(bcoef[2]/(2*bcoef[1]^2)) * exp(1 / (2*bcoef[1]*gsqr));

  return N*exp(value(0.5* cached_int[gsqr]))
end

@doc raw"""
    g_from_RG_eq(mu, bcoef::Vector{Float64};Lambda = MSbar_float().Lambda, nl=5,g0=1.0, c=0.5,verbose=false)

    g_from_RG_eq(mu, bcoef::Vector{uwreal};Lambda = MSbar_float().Lambda, nl=5,g0=1.0, c=0.5,verbose=false)  

It computes `g_s(\mu)` using the RG equation.`\mu` is assumed to be in MeV, `nl` is number of loops
at which one wants the beta function, `g0` is the intial guess for the
rooting routine (`root_error` from `ADerrors`), `c` shouldn't be changed.

The function only uses the first nl termns in bcoef.
`c` is an auxilary variable. If the rooting routine request `phi` computed at negative g
`g_from_RG_eq` catch the error and restart the rooting with initial guess `g0*c` and `c = 1.5*c`
In this case, if verbose=true, it prints a warning.

WARNING: At the moment `Lambda_MSBar = 341` without errors
"""
function g_from_RG_eq(mu, bcoef::Vector{Float64}; Lambda = MSbar_float().Lambda, nl=5,g0=1.0, c=0.5,verbose=false)  
  bcoef = bcoef[1:nl]
  mu_over_lambda = mu/Lambda
  gbar = 0.0;
  if mu_over_lambda isa uwreal
    aux_phi(g,p) = p[1] - phi(g,bcoef)
    try
      gbar =  root_error(aux_phi,g0,[mu_over_lambda])  
    catch e
      if e isa ArgumentError
        if verbose
          @warn """
          negative g requested.
          Automatic abortion of the rooting_routine
          and initial value changed to $(c*g0)
          """
        end
        gbar = g_from_RG_eq(mu,bcoef,Lambda=Lambda,nl=nl,g0 = c*g0, c = 0.5*c)
      else
        error(e)
      end
    end
  else
    f(g) =  mu_over_lambda - phi(g,bcoef)
    D(f) = x -> ForwardDiff.derivative(f,float(x))
    try
      gbar = Roots.find_zero((f,D(f)), g0, Roots.Newton())
    catch e 
      if e isa ArgumentError
        if verbose
          @warn """
          negative g requested.
          Automatic abortion of the rooting_routine
          and initial value changed to $(c*g0)
          """
        end
        gbar = g_from_RG_eq(mu,bcoef,Lambda = Lambda,nl=nl,g0 = c*g0, c = 0.5*c)
      elseif e isa DomainError
        throw(e)
      else
        @error e
        gbar = g_from_RG_eq(mu,bcoef,Lambda = Lambda,nl=nl,g0 = c*g0, c = 0.5*c)
      end
    end
  end
  return gbar
end

function g_from_RG_eq(mu, bcoef::Vector{uwreal};Lambda = MSbar.Lambda,nl=5,g0=1.0, c=0.5,verbose=false)
  if !(Lambda isa uwreal)
    @warn "Lambda is not an uwreal, are you sure you passing the correct arguments?"
  end
  if nl < length(bcoef) #check that ensure that the beta function is at the  correct loop approx
    bcoef = bcoef[1:nl];
  end

  mu_over_lambda = mu/Lambda
  
  if !(mu_over_lambda isa uwreal)
    mu_over_lambda = uwreal([mu_over_lambda,0.], "mu_over_lambda")
  end
  aux_phi(g,p) = p[1] - phi(g,p[2:end]);   
  gbar = 0.0
  try
    gbar =  root_error(aux_phi,g0,vcat(mu_over_lambda,bcoef))  
  catch e
    if e isa ArgumentError
      if verbose
        @warn """
        negative g requested.
        Automatic abortion of the rooting_routine
        and initial value changed to $(c*g0)
        """
      end
      gbar = g_from_RG_eq(mu,bcoef,nl=nl,g0 = c*g0, c = 1.5*c)
    else
      error(e)  
    end
  end
  return gbar
end

@doc raw"""
    alpha(mu,bcoef; nl=5, g0=1.0,c=0.5, verbose=false) 
  
    alpha(mu::AbstractArray,bcoef; nl=5,g0=1.0,c=0.5,verbose=false)

It returns `\alpha_s(\mu)` using the RG equation.
Assume `\mu` is in MeV.
See `g_from_RG_eq` for documentation.

    ```@example 
    bcoef = beta_function_coeff(3,3,nl=5)
    alpha_5l(10000,bcoef,nl=5)
    alpha_4l(10000,bcoef,nl=4)
    ```
"""
alpha(mu,bcoef;Lambda=MSbar.Lambda, nl=5, g0=1.0,c=0.5, verbose=false) = 
  g_from_RG_eq(mu,bcoef,Lambda =Lambda,nl=nl, g0=g0, c=c, verbose=verbose)^2/(4*pi)

function alpha(mu::AbstractArray,bcoef;Lambda=MSbar.Lambda, nl=5,g0=1.0,c=0.5,verbose=false)
  return [g_from_RG_eq(x,bcoef,Lambda =Lambda,nl=nl,g0=g0,c=c,verbose=verbose) for x in mu]
end


