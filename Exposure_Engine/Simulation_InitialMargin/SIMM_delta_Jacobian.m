function[SIMM_disc_sensitivity,SIMM_fwd_sensitivity] = SIMM_delta_Jacobian(shock_DF_disc,shock_DF_fwd,DF_disc,DF_fwd,V,id_instrument,ptf_data,schedule,evaluation_date,pricing_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise,Jacobian_disc,Jacobian_fwd_disc,Jacobian_fwd,fwd_curve_pillars,disc_curve_pillars)


n_pillar_disc    = size(shock_DF_disc,1);
n_pillar_fwd     = size(shock_DF_fwd,1);
n_path           = size(DF_fwd,1);
V_disc_shock     = zeros(n_path,n_pillar_disc);
V_fwd_shock      = zeros(n_path,n_pillar_fwd);
delta_disc       = zeros(n_path,n_pillar_disc);
delta_fwd        = zeros(n_path,n_pillar_fwd);

if disc_curve_pillars(1) == 0
    disc_curve_pillars   = disc_curve_pillars(2:end);
end
if fwd_curve_pillars(1)  == 0
    fwd_curve_pillars    = fwd_curve_pillars(2:end);
end


% Compute Sensitivities
% Forwarding curve
for i = 1:n_pillar_fwd
    V_fwd_shock(:,i)   = PortfolioExposure(id_instrument,ptf_data,schedule,shock_DF_fwd{i},DF_disc,evaluation_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise);
    delta_fwd(:,i)     = V_fwd_shock(:,i) - V;
end
fwd_curve_sensitivity  = delta_fwd * Jacobian_fwd;

% Discounting curve
for i = 1:n_pillar_disc
    V_disc_shock(:,i)  = PortfolioExposure(id_instrument,ptf_data,schedule,DF_fwd,shock_DF_disc{i},evaluation_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise);
    delta_disc(:,i)    = V_disc_shock(:,i) - V;
end
disc_curve_sensitivity = delta_disc * Jacobian_disc + delta_fwd * Jacobian_fwd_disc;


% Allocation on ISDA-SIMM pillars
SIMM_pillars          = (evaluation_date - pricing_date + [14 30 91 182 365 365*2 365*3 365*5 3650 365*15 365*20 365*30])./365;
SIMM_disc_sensitivity = zeros(n_path,length(SIMM_pillars));
SIMM_fwd_sensitivity  = zeros(n_path,length(SIMM_pillars));

% Forwarding curve
for i = 1:n_pillar_fwd
    index_min = max(find(SIMM_pillars <= fwd_curve_pillars(i)));
    index_max = min(find(SIMM_pillars >= fwd_curve_pillars(i)));
    
    if index_min == index_max % = j-th pillar
        SIMM_fwd_sensitivity(:,index_max)      = SIMM_fwd_sensitivity(:,index_max) + fwd_curve_sensitivity(:,i);
        
    else % allocation
        
        if isempty(index_min) == true       % allocated in the first pillar 
            SIMM_fwd_sensitivity(:,index_max)  = SIMM_fwd_sensitivity(:,index_max) + fwd_curve_sensitivity(:,i);
            
        elseif isempty(index_max) == true   % allocated in the last pillar
            SIMM_fwd_sensitivity(:,index_min)  = SIMM_fwd_sensitivity(:,index_min) + fwd_curve_sensitivity(:,i);
            
        else % proportional allocation
            alloc_min = (SIMM_pillars(index_max)-fwd_curve_pillars(i))/(SIMM_pillars(index_max)-SIMM_pillars(index_min));
            alloc_max = (fwd_curve_pillars(i)-SIMM_pillars(index_min))/(SIMM_pillars(index_max)-SIMM_pillars(index_min));
            
            SIMM_fwd_sensitivity(:,index_min)  = SIMM_fwd_sensitivity(:,index_min) + fwd_curve_sensitivity(:,i)*alloc_min;
            SIMM_fwd_sensitivity(:,index_max)  = SIMM_fwd_sensitivity(:,index_max) + fwd_curve_sensitivity(:,i)*alloc_max;

        end   
    end
end
            
% Discounting curve
for i = 1:n_pillar_disc
    index_min = max(find(SIMM_pillars <= disc_curve_pillars(i)));
    index_max = min(find(SIMM_pillars >= disc_curve_pillars(i)));
    
    if index_min == index_max % = j-th pillar
        SIMM_disc_sensitivity(:,index_max)    = SIMM_disc_sensitivity(:,index_max) + disc_curve_sensitivity(:,i);
        
    else % allocation
        
        if isempty(index_min) == true       % allocated in the first pillar 
            SIMM_disc_sensitivity(:,index_max) = SIMM_disc_sensitivity(:,index_max) + disc_curve_sensitivity(:,i);
            
        elseif isempty(index_max) == true   % allocated in the last pillar
            SIMM_disc_sensitivity(:,index_min) = SIMM_disc_sensitivity(:,index_min) + disc_curve_sensitivity(:,i);
            
        else % proportional allocation
            alloc_min = (SIMM_pillars(index_max)-disc_curve_pillars(i))/(SIMM_pillars(index_max)-SIMM_pillars(index_min));
            alloc_max = (disc_curve_pillars(i)-SIMM_pillars(index_min))/(SIMM_pillars(index_max)-SIMM_pillars(index_min));
            
            SIMM_disc_sensitivity(:,index_min) = SIMM_disc_sensitivity(:,index_min) + disc_curve_sensitivity(:,i)*alloc_min;
            SIMM_disc_sensitivity(:,index_max) = SIMM_disc_sensitivity(:,index_max) + disc_curve_sensitivity(:,i)*alloc_max;

        end   
    end
end
            
end



