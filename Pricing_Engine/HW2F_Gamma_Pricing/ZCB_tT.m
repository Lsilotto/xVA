function[ZCB_tT] = ZCB_tT(r_mkt,t,T,DF_mkt_tau,interp_method,param,x_t,y_t)

rM_0t  = interp1(DF_mkt_tau,r_mkt,t,interp_method,'extrap');
PM_0t  = exp(-rM_0t*t);

rM_0T  = interp1(DF_mkt_tau,r_mkt,T,interp_method,'extrap');
PM_0T  = exp(-rM_0T.*T);

for i = 1:size(PM_0T,2)
    ZCB_tT(:,i) = ZCB_HW2F_Gamma(t,T(i),param,PM_0T(i),PM_0t,x_t,y_t);
end
