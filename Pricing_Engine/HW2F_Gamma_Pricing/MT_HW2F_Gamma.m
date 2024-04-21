function [M]=MT_HW2F_Gamma(a,sigma,b,eta,rho,Ti,Gamma_i,nGamma,T,s,t)
 
% The function computes M for G2++; params a-sigma-b-eta-rho,
% T forward measure
% s first time, t second time

% G2++ parameters:
% a=Pars(1,1);
% sigma=Pars(1,2);
% b=Pars(1,3);
% eta=Pars(1,4);
% rho=Pars(1,5);
% Ti=Pars(2:end,1);
% Gamma_i=Pars(2:end,2);
% nGamma=size(Ti,1);
 
%BtT=BtT_HW2F(a,s,t);
ns=size(s,1);
nt=size(t,1);
M=zeros(nt,ns);
 
for j=1:ns
% first time
 T_imin=s(j);
 index_min=find(Ti<s(j));
  
 if isempty(index_min)==1
 % s<=Ti(1)
 Gamma_imin=Gamma_i(1);
 index_min(1)=0;
 else
      
 if size(index_min,1)==nGamma
 Gamma_imin=Gamma_i(nGamma);
 else
 Gamma_imin=Gamma_i(index_min(end)+1);
 end
  
 end
  
    for k=1:nt  
         
    T_imax=t(k);
    index_max=find(Ti> t(k));
    if isempty(index_max)==1
    % t(k)>=Ti(end)
    Gamma_imax=Gamma_i(end);    
    index_max=nGamma;
    else
    Gamma_imax=Gamma_i(index_max(1));
    end 
    Gamma_T=Gamma_i(index_min(end)+1:index_max(1)-1,1);
    T_T=Ti(index_min(end)+1:index_max(1)-1,1);
    Gamma_T=[Gamma_T;Gamma_imax];
    T_T=[T_imin;T_T;T_imax];
    
    x=1/a*exp(-a*t(k))*(exp(a*T_T(2:end))-exp(a*T_T(1:end-1)))...
        - exp(-a*(t(k)+T))*(exp(2*a*T_T(2:end))-exp(2*a*T_T(1:end-1)))/(2*a);
    
    M(k,j)=sigma^2/a*sum(Gamma_T.^2.*x);
    
    x=1/a*exp(-a*t(k))*(exp(a*T_T(2:end))-exp(a*T_T(1:end-1)))... 
         - exp(-(a*t(k)+b*T))*(exp((a+b)*T_T(2:end))-exp((a+b)*T_T(1:end-1)))/(a+b);
     
    M(k,j)= M(k,j) + rho*sigma*eta/b*sum(Gamma_T.^2.*x);
    end
end