using RG_running, ADerrors

mu   = 1000 #MeV
beta = beta_function_coeff(3,3)
a = alpha(mu,beta)

mu_AD = uwreal([1000.,10.],"mu");
a_AD = alpha(mu_AD,beta)
uwerr(a_AD)


print("""
    alpha_s(1GeV)       = $a;
    alpha_s(1(0.01)GeV) = $(a_AD);
""")
