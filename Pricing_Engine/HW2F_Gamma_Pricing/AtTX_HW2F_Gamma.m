function A=AtTX_HW2F_Gamma(Pars,t,T)
% The function computes A = exp(1/2*(V(t,T)-V(0,T)+V(0,t)))
% ZCB(t,T)= DF_M(T)/DF_M(t)*A*exp(-B(a,t,T)*x_t)
% input Pars=[a sigma b eta rho] or [T1 Gamma1_T T1 Gamma2_T zeros(size(Gamma))]
% 	t: first time
% 	T: second time
% output: matrix nTxnt

a=Pars(1,1);
sigma=Pars(1,2);
b=Pars(1,3);
eta=Pars(1,4);
rho=Pars(1,5);

nT=size(T,1);
nt=size(t,1);
A=zeros(nT,nt);
V=V_HW2F_Gamma(Pars,t,T);
V1=V_HW2F_Gamma(Pars,0,T);
V2=V_HW2F_Gamma(Pars,0,t);
% first matrix is nTxnt, second nTx1, third ntX1;

B=0.5*(V-repmat(V1,1,nt)+repmat(V2',nT,1));
A=exp(B);
