function[SIMM_disc_sensitivity,SIMM_fwd_sensitivity] = SIMM_delta(shock_DF_disc,shock_DF_fwd,DF_disc,DF_fwd,V,id_instrument,ptf_data,schedule,evaluation_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise)

n_pillar_SIMM         = size(shock_DF_disc,1);
n_path                = size(shock_DF_disc{1},1);
V_disc_shock          = zeros(n_path,n_pillar_SIMM);
V_fwd_shock           = zeros(n_path,n_pillar_SIMM);
SIMM_disc_sensitivity = zeros(n_path,length(n_pillar_SIMM));
SIMM_fwd_sensitivity  = zeros(n_path,length(n_pillar_SIMM));


for i = 1:n_pillar_SIMM
    
    % compute delta for discounting curve 
    V_disc_shock(:,i)            = PortfolioExposure(id_instrument,ptf_data,schedule,DF_fwd,shock_DF_disc{i},evaluation_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise);
    SIMM_disc_sensitivity(:,i)   = V_disc_shock(:,i) - V;
    % compute delta for forwarding curve 
    V_fwd_shock(:,i)             = PortfolioExposure(id_instrument,ptf_data,schedule,shock_DF_fwd{i},DF_disc,evaluation_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise);
    SIMM_fwd_sensitivity(:,i)    = V_fwd_shock(:,i) - V;    
        
end



