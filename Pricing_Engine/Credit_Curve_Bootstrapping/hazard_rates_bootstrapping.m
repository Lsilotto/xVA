function hazard_rates = hazard_rates_bootstrapping(evaluation_date,effective_date,first_coupon_date,cds_curve,period, basis,recovery_rate, yields)

options = optimset;
options.Display = 'off';
options.TolX = 10^-10;

hazard_rates=[];
siz=size(cds_curve,1);

for i=1:siz   
    maturity=datestr(datenum(evaluation_date)+cds_curve(i,1));
    coupon=cds_curve(i,2);
    last_hazard_rate=fzero(@(last_hazard_rate) cds_price(evaluation_date,effective_date,maturity,first_coupon_date,coupon,period, basis,recovery_rate, yields, hazard_rates, last_hazard_rate),0.001);
    hazard_rates=[hazard_rates;cds_curve(i,1),last_hazard_rate ]; 
end    
