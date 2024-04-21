function MtM = cds_price(evaluation_date,effective_date,maturity,first_coupon_date,coupon,period, basis,recovery_rate, yields, hazard_rates, last_hazard_rate)


evaluation_date=datenum(evaluation_date);
effective_date=datenum(effective_date);
maturity=datenum(maturity);
first_coupon_date=datenum(first_coupon_date);
coupon_dates_pl= busdate(cfdates(effective_date, maturity, period, basis,NaN,NaN,first_coupon_date));
coupon_dates_pl(end)=min(coupon_dates_pl(end),maturity);    
coupon_dates_dl= busdate(cfdates(effective_date, maturity, period, basis));
coupon_dates_dl(end)=min(coupon_dates_dl(end),maturity);
coupon_terms_pl=days365(evaluation_date,[effective_date,coupon_dates_pl]);
coupon_terms_dl=days365(evaluation_date,[effective_date,coupon_dates_dl]);
discount_factors_pl=rate2disc(-1,interp1(yields(:,1),yields(:,2),coupon_terms_pl,'linear')',yearfrac(evaluation_date,[effective_date,coupon_dates_pl],3)');
discount_factors_dl=rate2disc(-1,interp1(yields(:,1),yields(:,2),coupon_terms_dl,'linear')',yearfrac(evaluation_date,[effective_date,coupon_dates_dl],3)');
full_hazard=[hazard_rates;days365(evaluation_date,maturity),last_hazard_rate];
surv_prob_pl=survival_probability(coupon_terms_pl', full_hazard);
surv_prob_dl=survival_probability(coupon_terms_dl', full_hazard);
PL =coupon/10000*sum(diff([0,coupon_terms_pl]/360).*discount_factors_pl'.*surv_prob_pl);
DL = (1 - recovery_rate/100)*sum((-diff([1,surv_prob_dl])).*discount_factors_dl');
MtM =PL-DL;







