function [Mx,My]=MT_xy_HW2_Gamma(Pars,T,s,t)
 
% The function compute MT_xy with fwd measure T
 
% input Pars=[a sigma b eta rho] or [T1 Gamma1_T T1 Gamma2_T 0]
%   T: scalar, fwd measure
%   s: first time, vector ns
%   t: second time, vector nt
%   output: matrix nt x ns


a    = Pars(1,1);
sigma= Pars(1,2);
b    = Pars(1,3);
eta  = Pars(1,4);
rho  = Pars(1,5);
Ti   = Pars(2:end,1);
Gamma_i = Pars(2:end,2);
nGamma = size(Ti,1);


Mx=MT_HW2F_Gamma(a,sigma,b,eta,rho,Ti,Gamma_i,nGamma,T,s,t);
My=MT_HW2F_Gamma(b,eta,a,sigma,rho,Ti,Gamma_i,nGamma,T,s,t);
