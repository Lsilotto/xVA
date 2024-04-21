function V=V_HW2F(Pars,t,T)

% function for the variance for the primitive of x_t+y_t
% input Pars=[a sigma b eta rho]
% Gamma_1 = Gamma_2 - avoid overparametrization
%   t: first time
%   T sec time
% output: array nTxnt

% pag 145 Brigo-Mercurio
a=Pars(1,1);
sigma=Pars(1,2);
b=Pars(1,3);
eta=Pars(1,4);
rho=Pars(1,5);

nT=size(T,1);
nt=size(t,1);
V=zeros(nT,nt);

T1=repmat(T,1,nt);
t1=repmat(t',nT,1);
DT=T1-t1;
eadt=exp(-a*DT);
ebdt=exp(-b*DT);
s1=DT+2*eadt/a-0.5*eadt.*eadt/a-1.5/a;
s2=DT+2*ebdt/b-0.5*ebdt.*ebdt/b-1.5/b;
s3=DT+(eadt-1)/a+(ebdt-1)/b-(eadt.*ebdt-1)/(a+b);

V = sigma^2*s1/a^2+eta^2*s2/b^2+2*rho*eta*sigma*s3/(a*b);



