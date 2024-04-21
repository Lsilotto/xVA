function V=V_HW2F_Gamma(Pars,t,T)
 
% function for the variance for the primitive of x_t+y_t
% input Pars=[a sigma b eta rho] or [T1 Gamma1_T T1 Gamma2_T 0]
% Gamma_1 = Gamma_2 - avoid overparametrization
%   t: first time
%   T sec time
% output: array nTxnt
 
a=Pars(1,1);
sigma=Pars(1,2);
b=Pars(1,3);
eta=Pars(1,4);
rho=Pars(1,5);
Ti=Pars(2:end,1);
Gamma_i=Pars(2:end,2);
nGamma=size(Ti,1);
 
 
nT=size(T,1);
nt=size(t,1);
V=zeros(nT,nt);
 
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
    x=diff(T_T) -2/a*exp(-a*T(k))*(exp(a*T_T(2:end))-exp(a*T_T(1:end-1)))+exp(-2*a*T(k))*(exp(2*a*T_T(2:end))-exp(2*a*T_T(1:end-1)))/(2*a);
    % 1st factor vol
    V(k,j)=sigma^2/a^2*sum(Gamma_T.^2.*x);
    x=diff(T_T) -2/b*exp(-b*T(k))*(exp(b*T_T(2:end))-exp(b*T_T(1:end-1)))+exp(-2*b*T(k))*(exp(2*b*T_T(2:end))-exp(2*b*T_T(1:end-1)))/(2*b);
    % 2nd
    V(k,j)=eta^2/b^2*sum(Gamma_T.^2.*x)+V(k,j);
    
% %     
    % correlations
    x=diff(T_T) -1/a*exp(-a*T(k))*(exp(a*T_T(2:end))-exp(a*T_T(1:end-1)))-1/b*exp(-b*T(k))*(exp(b*T_T(2:end))-exp(b*T_T(1:end-1)));
    x=x+exp(-(a+b)*T(k))*(exp((a+b)*T_T(2:end))-exp((a+b)*T_T(1:end-1)))/(a+b);
    V(k,j)=2*eta*sigma*rho/(a*b)*sum(Gamma_T.^2.*x)+V(k,j);
    end
         
end
