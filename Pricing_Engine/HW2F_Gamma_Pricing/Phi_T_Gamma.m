function [r0,phi] = Phi_T_Gamma(Params,fwd_rate,Gamma,T)

a     = Params(1);
sigma = Params(2);
b     = Params(3);
eta   = Params(4);
rho   = Params(5);

sigma_t = Gamma*sigma;
eta_t   = Gamma*eta;

phi    = fwd_rate + (B_HW2F(a,0,T).^2).*((sigma_t.^2)/2) + (B_HW2F(b,0,T).^2).*((eta_t.^2)/2) + B_HW2F(a,0,T).*B_HW2F(b,0,T).*rho.*sigma_t.*eta_t;
r0     = phi(1);
end