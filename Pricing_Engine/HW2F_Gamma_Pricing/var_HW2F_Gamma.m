function sigma = var_HW2F_Gamma(Params,Gamma,dt)

sigma = ((Gamma*Params(2)).^2)*(1-exp(-2*Params(1)*(dt)))/(2*Params(1));

end