function VarX=VarY_HW2F_Gamma(Pars,t,T)
 
% function for the variance of x_t+y_t
% input     Pars=[[a sigma] ; [Ti Gamm_Ti]]
% Gamma stepwise, value Gamma_Ti(1) from 0 to Ti.
%   t: first time
%   T: second time
% output: matrice nT x nt
 

b=Pars(1,1);
eta=Pars(1,2);

Ti=Pars(2:end,1);
Gamma_i=Pars(2:end,2);
nGamma=size(Ti,1);
 
nT=size(T,1);
nt=size(t,1);
VarX=zeros(nT,nt);
 
% t(j) fixed
for j=1:nt
 T_imin=t(j);
 index_min=find(Ti<t(j));
 if isempty(index_min)==1
 % t<=Ti(1)
 Gamma_imin=Gamma_i(1);
 index_min(1)=0;
 else
 if size(index_min,1)==nGamma
 Gamma_imin=Gamma_i(nGamma);
 else
 Gamma_imin=Gamma_i(index_min(end)+1);
 end
 end
    for k=1:nT  
    T_imax=T(k);
    index_max=find(Ti>T(k));
    if isempty(index_max)==1
    % T(k)>=Ti(end)
    Gamma_imax=Gamma_i(end);    
    index_max=nGamma;
    else
    Gamma_imax=Gamma_i(index_max(1));
    end 
    Gamma_T=Gamma_i(index_min(end)+1:index_max(1)-1,1);
    T_T=Ti(index_min(end)+1:index_max(1)-1,1);
    Gamma_T=[Gamma_T;Gamma_imax];
    T_T=[T_imin;T_T;T_imax];
     
    x=exp(2*b*T_T(2:end))-exp(2*b*T_T(1:end-1));
      VarX(k,j)=0.5*eta^2/b*sum(Gamma_T.^2.*x)*exp(-2*b*T(k));
     
    end     
end