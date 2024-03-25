@doc raw"""
  beta_function_coeff(nc,nf;nl = 5)

  Computes and return beta function coefficient. 
  `nc` is the number of color, `nf` is the number of active fermions
  `nl` is the number of loops at which one want the beta function. 
  
  The coefficient up to beta_4 are taken from 1701.01404, in particular they where 
  taken from the file attached to the paper in arXiv. 
  
  The function always computes all the coefficient but returns only the firt nl terms 
"""
function beta_function_coeff(nc, nf; nl=5)
  zeta3 = 1.20205690315959428539; zeta4 = pi^4/90;  zeta5 = 1.03692775514336992633;

  tf = 1/2; ca = nc; cf = (nc^2 - 1)/(2*nc); na = nc^2 - 1; nr = nc; 
  
  dada = nc^2/24*(nc^2 + 36)*na;
  dfda = (1/48)*(nc - 1)*(nc + 1)*(nc^2 + 6)*nr;
  dfdf = (nc - 1)*(nc + 1)*(nc^4 - 6*nc^2 + 18)/96/nc^3*nr; 

  bcoefs = [((11/3*ca - 4/3*tf*nf)/(4*pi)^2), 
  
  (1/(4*pi)^4)*(34/3*ca^2 - 4*cf*tf*nf - 20/3*ca*tf*nf), 
  
  (1/(4*pi)^6)*(2857/54*ca^3 + 2*cf^2*tf*nf - 205/9*cf*ca*tf*nf - 1415/27*ca^2*tf*nf + 44/9*cf*tf^2*nf^2 + 158/27*ca*tf^2*nf^2), 
  
  (1/(4*pi)^8)*(ca^4*(150653/486 - 44/9*zeta3) + ca^3*tf*nf*(-39143/81 + 136/3*zeta3) + ca^2*cf*tf*nf*(7073/243 - 656/9*zeta3) +ca*cf^2*tf*nf*(-4204/27 + 352/9*zeta3) + 46*cf^3*tf*nf +
  ca^2*tf^2*nf^2*(7930/81 + 224/9*zeta3) + cf^2*tf^2*nf^2*(1352/27 - 704/9*zeta3) + ca*cf*tf^2*nf^2*(17152/243 + 448/9*zeta3) + 424/243*ca*tf^3*nf^3 +
  1232/243*cf*tf^3*nf^3 + (dada/na)*(-80/9 + 704/3*zeta3) + nf*(dfda/na)*(512/9 - 1664/3*zeta3) +  nf^2*(dfdf/na)*(-704/9 + 512/3*zeta3)), 
  
  (1/(4*pi)^10)*(ca^5 * ( 8296235/3888 - 1630/81*zeta3 + 121/6*zeta4 - 1045/9*zeta5 ) + ca*(dada/na) * ( - 514/3 + 18716/3*zeta3 - 968*zeta4 - 15400/3*zeta5 ) 
  + nf*tf*ca^4 * ( - 5048959/972 + 10505/81*zeta3 - 583/3*zeta4 + 1230*zeta5 ) + nf*tf*cf*ca^3 * ( 8141995/1944 + 146*zeta3 + 902/3*zeta4 - 8720/3*zeta5 )
  + nf*tf*cf^2*ca^2 * ( - 548732/81 - 50581/27*zeta3 - 484/3*zeta4 + 12820/3*zeta5 )  + nf*tf*cf^3*ca * ( 3717 + 5696/3*zeta3 - 7480/3*zeta5 ) + nf*tf*cf^4 * ( - 4157/6 - 128*zeta3 )
  + nf*ca*(dfda/na) * ( 11312/9 - 127736/9*zeta3 + 2288*zeta4 + 67520/9*zeta5 )  + nf*cf*(dfda/na) * ( - 320 + 1280/3*zeta3 + 6400/3*zeta5 )
  + nf*tf*(dada/na) * ( 904/9 - 20752/9*zeta3 + 352*zeta4 + 4000/9*zeta5 )  + nf^2*tf^2*ca^3 * ( 843067/486 + 18446/27*zeta3 - 104/3*zeta4 - 2200/3*zeta5 )
  + nf^2*tf^2*cf*ca^2 * ( 5701/162 + 26452/27*zeta3 - 944/3*zeta4 + 1600/3*zeta5 )  + nf^2*tf^2*cf^2*ca * ( 31583/18 - 28628/27*zeta3 + 1144/3*zeta4 - 4400/3*zeta5 )
  + nf^2*tf^2*cf^3 * ( - 5018/9 - 2144/3*zeta3 + 4640/3*zeta5 )+ nf^2*ca*(dfdf/na) * ( - 7184/3 + 40336/9*zeta3 - 704*zeta4 + 2240/9*zeta5 )
  + nf^2*cf*(dfdf/na) * ( 4160/3 + 5120/3*zeta3 - 12800/3*zeta5 )  + nf^2*tf*(dfda/na) * ( - 3680/9 + 40160/9*zeta3 - 832*zeta4 - 1280/9*zeta5 )
  + nf^3*tf*(dfdf/na) * ( 3520/9 - 2624/3*zeta3 + 256*zeta4 + 1280/3*zeta5 )  + nf^3*tf^3*ca^2 * ( - 2077/27 - 9736/81*zeta3 + 112/3*zeta4 + 320/9*zeta5 )
  + nf^3*tf^3*cf*ca * ( - 736/81 - 5680/27*zeta3 + 224/3*zeta4 )  + nf^3*tf^3*cf^2 * ( - 9922/81 + 7616/27*zeta3 - 352/3*zeta4 )
  + nf^4*tf^4*ca * ( 916/243 - 640/81*zeta3 )  + nf^4*tf^4*cf * ( - 856/243 - 128/27*zeta3 ))];
  
  return bcoefs[1:nl]
end

#= function __phi(g,bcoef::Vector;)
  key = bcoef[end]
  if !(key in keys(__integral__))
    __integral__[key] = Dict{Float64,uwreal}()
  end
  cached_int = __integral__[key]
  aux1(x,p) = sum([x^(2*(i-1))*p[i] for i=eachindex(p)])  ## beta = x^3 * aux1
  aux2(x,p) = sum([x^(2*(i-1)-3)*p[i] for i = 3:lastindex(p)])
  integrand(x,p) =((p[1]-p[2]*x^2)*aux2(x,p) - p[2]^2*x)/(p[1]^2*aux1(x,p))

  g_used = [k for (k,_) in cached_int]
  if !(g in g_used)
    aux = findall(x->x<g,g_used)
    lbound = isempty(aux) ? 0 : maximum(g_used[aux])
    _bcoef = [uwreal([b,0.],"beta coeff") for b in bcoef];
    cached_int[g] = lbound == 0.0 ? int_error(integrand,0,g,_bcoef) : (cached_int[lbound] + int_error(integrand,lbound,g,_bcoef))
  end
  int_val = cached_int[g]
  N = (bcoef[1] * g^2)^(bcoef[2]/(2*bcoef[1]^2)) * exp(1 / (2*bcoef[1]*g^2));
  return N*exp(value(int_val))
end =#

@doc raw"""
  phi(g,bcoef::Vector)

  it computes the function $(L"$\phi(y=g_s^2)$"). For numerical stability it was decided to 
  compute it directly in terms of $(L"$y = g_s^2$"). For rapidity, the value of the integrals is
  stored into a cache variable `__integral_phi__`. When `phi(g_s,bcoef)` is called, it first 
  check if the integral was already computed. If so it retrieve the cached results and return `phi`,
  otherwise it looks for the highest `g<g_s` for which the integral was computed, retrieve 
  `cached_int[g]` and compute the integral from g to g_s. The sum of these is cached in 
  `__integral_phi__`, such that each entry is the integral from `0` to the respective `g`.
  
  WARNING: `phi` alway return a Float64, and `bcoef` must not be an uwreal. g must be positive
  otherwise an `ArgumentError` is thrown
"""
function phi(g,bcoef::Vector)
  key = bcoef[end]
  if !(key in keys(__integral_phi__))
    __integral_phi__[key] = Dict{Float64,uwreal}()
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
    lbound = isempty(aux) ? 0 : maximum(gsqr_used[aux])
    cached_int[gsqr] = lbound == 0.0 ? int_error(integrand,0,gsqr,uwreal.(bcoef)) : (cached_int[lbound] + int_error(integrand,lbound,gsqr ,uwreal.(bcoef)))
  end
  int_val = cached_int[gsqr]
  N = (bcoef[1] * gsqr)^(bcoef[2]/(2*bcoef[1]^2)) * exp(1 / (2*bcoef[1]*gsqr));

  return N*exp(value(0.5*int_val))
end

@doc raw"""
  g\_from\_RG\_eq(mu, bcoef::Vector{Float64}; nl=5,g0=1.0, c=0.5,verbose=false)
  
  g\_from\_RG\_eq(mu, bcoef::Vector{uwreal}; nl=5,g0=1.0, c=0.5,verbose=false)  

  It computes `g_s(\mu)` using the RG equation.`\mu` is assumed to be in MeV `nl` is number of loops
  at which one wants the beta function, `g0` is the intial guess for the
  rooting routine (`root_error` from `ADerrors`), `c` shouldn't be changed.
  
  The function is buildt such that it checks whether your are passing the correct
  number of beta coefficient. If more coefficient are given, it keeps only the first `nl`.
  In this way you can store only one vector of beta coefficient and modify the approximation
  now. 
  `c` is an auxilary variable. If the rooting routine request `phi` computed at negative g
  `g_from_RG_eq` catch the error and restart the rooting with initial guess `g0*c` and `c = 1.5*c`
  In this case, if verbose=true, it prints a warning.

  WARNING: At the moment `Lambda_MSBar = 341` without errors
"""
function g_from_RG_eq(mu, bcoef::Vector{Float64}; nl=5,g0=1.0, c=0.5,verbose=false)  
  if nl < length(bcoef) #check that ensure that the beta function is at the  correct loop approx
    bcoef = bcoef[1:nl];
  end
  mu_over_lambda = mu/uwreal([341.0,0.0],"Lambda_MSBar")
  aux_phi(g,p) = p[1] - phi(g,bcoef)
  gbar = 0.0;
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
      gbar = g_from_RG_eq(mu,bcoef,nl=nl,g0 = c*g0, c = 1.5*c)
    else
      error(e.msg)
    end
  end
  return gbar
end

function g_from_RG_eq(mu, bcoef::Vector{uwreal};nl=5,g0=1.0, c=0.5,verbose=false)
  if nl < length(bcoef) #check that ensure that the beta function is at the  correct loop approx
    bcoef = bcoef[1:nl];
  end
  mu_over_lambda = mu/uwreal([341.0,0.0],"Lambda_MSBar")
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
      error(e.msg)  
    end
  end
  return gbar
end

@doc raw"""
  alpha(mu,bcoef; nl=5, g0=1.0,c=0.5, verbose=false) 
  
It returns `\alpha_s(\mu)` using the RG equation.
Assume `\mu` is in MeV.
See `g_from_RG_eq` for documentation.

```@example 
bcoef = beta_function_coeff(3,3,nl=5)
alpha_5l(10000,bcoef,nl=5)
alpha_4l(10000,bcoef,nl=4)
```
"""
alpha(mu,bcoef; nl=5, g0=1.0,c=0.5, verbose=false) = 
  g_from_RG_eq(mu,bcoef,nl=nl, g0=g0, c=c, verbose=verbose)^2/(4*pi)