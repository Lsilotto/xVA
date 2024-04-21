function [ZCB_tT] = ZCB_HW2F_Gamma(t,T,param,DF_mkt_0T,DF_mkt_0t,x_t,y_t)
% The function compute ZCB(t,T) for G2++. 
% Formula 4.14 pag. 146 Brigo-Mercurio        

a           = param(1,1);
b           = param(1,3);

V_tT        = V_HW2F_Gamma(param,t,T);
V_0T        = V_HW2F_Gamma(param,0,T);
V_0t        = V_HW2F_Gamma(param,0,t);


A_tT        = 0.5.*(V_tT-V_0T+V_0t)-((1-exp(-a.*(T-t)))/a).*x_t-((1-exp(-b.*(T-t)))/b).*y_t;
ZCB_tT      = DF_mkt_0T./DF_mkt_0t.*exp(A_tT);