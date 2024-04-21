function CDSCurve = CDSCurve(CDSCurveIN,recoveryRate,curveDate,firstCouponDate,maxDate,...
period, basis,yieldCurve)

% Build output curves
hazardRates = hazard_rates_bootstrapping(curveDate,curveDate,firstCouponDate,CDSCurveIN, period, basis,recoveryRate, yieldCurve);
tauD = [1:1:maxDate-curveDate]';      % output date grid (in days from curveDate)
survivalProbabilities=survival_probability(tauD,hazardRates);   
   
CDSCurve.cdsCurve                    = CDSCurveIN;
CDSCurve.curveDate                   = curveDate;
CDSCurve.tau                         = tauD;
CDSCurve.hazardRates                 = hazardRates;
CDSCurve.survivalProbabilities       = survivalProbabilities';
