function[DF_curve_CF_dates,r_out] = DF_CF_dates(DF_curve,DF_curve_pillars,CF_dates_grid,grid_point,interp_method,r_in)

tau_DF_curve                   = DF_curve_pillars - grid_point;
tau_CF                         = CF_dates_grid - grid_point;

if nargin < 6
    r_out                      = zcb2r(DF_curve,tau_DF_curve);
    
    if verLessThan ('matlab','9.4') == 0
        DF_curve_CF_dates          = exp(-interp1(tau_DF_curve,r_out',tau_CF,interp_method,'extrap')'.*tau_CF);
    else
        DF_curve_CF_dates          = exp(bsxfun(@times,-interp1(tau_DF_curve,r_out',tau_CF,interp_method,'extrap')',tau_CF));
    end
    
else
    if verLessThan ('matlab','9.4') == 0
        DF_curve_CF_dates          = exp(-interp1(tau_DF_curve,r_in',tau_CF,interp_method,'extrap')'.*tau_CF);
    else
        DF_curve_CF_dates          = exp(bsxfun(@times,-interp1(tau_DF_curve,r_in',tau_CF,interp_method,'extrap')',tau_CF));
    end
   
end


