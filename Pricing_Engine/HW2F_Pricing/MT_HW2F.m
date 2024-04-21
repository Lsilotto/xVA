function [M]=MT_HW2F(a,sigma,T,s,t,b,eta,rho)
 
% The function computes M for G2++; params a-sigma-b-eta-rho,
% T forward measure
% s first time, t second time
 
BtT=BtT_HW2F(a,s,t);
ns=size(s,1);
nt=size(t,1);
M=zeros(nt,ns);
 
t1=repmat(t,1,ns);
s1=repmat(s',nt,1);
 
X1=BtT*(sigma^2/a+rho*sigma*eta/b);
X2=0.5*sigma^2/a^2.*(exp(-a*(T-t1))-exp(-a*(T+t1-2*s1)));
X3=rho*sigma*eta/(b*(a+b)).*(exp(-b*(T-t1))-exp(-b*T-a*t1+(a+b)*s1));
M=X1-X2-X3;