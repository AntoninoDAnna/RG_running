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

struct HyperPar
  F2::Union{Float64,uwreal}
  F3::Union{Float64,uwreal}
  Lambda::Union{Float64,uwreal}
end

function msbar()
  F2=uwreal([1.7505,0.0089], "F2")
  F3_val=0.5226 
  Lambda_val=341.0
  Lambda_F3_cov= [[150.7862014579557, -0.01716315895017156] [ -0.01716315895017156, 0.00001830314004469306]]
  vals = cobs([Lambda_val, F3_val], Lambda_F3_cov, "Lambda_MSBar & F3")
  return HyperPar(F2, vals[2], vals[1])
end

MSbar = msbar()
MSbar_Float = HyperPar(1.7505,0.5226,341.0)