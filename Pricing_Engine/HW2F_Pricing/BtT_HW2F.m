function BtT=BtT_HW2F(x,t,T)
% The function computes B(a,t,T)
% input:x G2++ mean reversion parameter
% 		t first time
%		T second time
% output	matrix nT x nt

% T>t

% pag. 148 Brigo-Mercurio
nT=size(T,1);
nt=size(t,1);
BtT=zeros(nT,nt);
BtT=(1-exp(-x.*(repmat(T,1,nt)-repmat(t',nT,1))))./x;
