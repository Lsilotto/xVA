function [Mx,My]=MT_xy_HW2(Pars,T,s,t)

% The function compute MT_xy with fwd measure T
 
% input Pars=[a sigma b eta rho]
%   T: scalar, fwd measure
%   s: first time, vector ns
%   t: second time, vector nt
%   output: matrix nt x ns

a=Pars(1,1);
sigma=Pars(1,2);
b=Pars(1,3);
eta=Pars(1,4);
rho=Pars(1,5);

Mx=MT_HW2F(a,sigma,T,s,t,b,eta,rho);
My=MT_HW2F(b,eta,T,s,t,a,sigma,rho);
