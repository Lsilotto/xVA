function[vega_margin,curvature_margin,vega] = SIMM_vega_HW2F(parameters,V,swap_rates,swap_annuity,displacement,evaluation_date,id_instruments,ptf_data,schedules,DF_fwd_sim,DF_disc_sim,libor_fixing_in,date,swpt_exercise,FX_spot,check_ctp_IM)


n_sim                 = size(V,1); 
n_assets              = size(V,2);
swapt_expiry          = (ptf_data(:,2) - evaluation_date)/365;
strike                = ptf_data(:,4);
SIMM_tenors           = [14 30 91 182 365 365*2 365*3 365*5 3650 365*15 365*20 365*30]./365;
delta_V               = zeros(n_sim,n_assets);
vega                  = cell(1,length(id_instruments));
VR                    = cell(length(SIMM_tenors),n_assets);
CVR                   = cell(length(SIMM_tenors),n_assets);
tmp_net_VR            = cell(1,length(SIMM_tenors));
tmp_net_CVR           = cell(1,length(SIMM_tenors));


% Define ISDA-SIMM parameters
%+++++++++++ V 2.1
correlation_matrix    = [1 0.63 0.59 0.47 0.31 0.22 0.18 0.14 0.09 0.06 0.04 0.05;...
                         0.63 1 0.79 0.67 0.52 0.42 0.37 0.30 0.23 0.18 0.15 0.13;...
                         0.59 0.79 1 0.84 0.68 0.56 0.50 0.42 0.32 0.26 0.24 0.21;...
                         0.47 0.67 0.84 1 0.86 0.76 0.69 0.60 0.48 0.42 0.38 0.33;...
                         0.31 0.52 0.68 0.86 1 0.94 0.89 0.80 0.67 0.60 0.57 0.53;...
                         0.22 0.42 0.56 0.76 0.94 1 0.98 0.91 0.79 0.73 0.70 0.66;...
                         0.18 0.37 0.50 0.69 0.89 0.98 1 0.96 0.87 0.81 0.78 0.74;...
                         0.14 0.30 0.42 0.60 0.80 0.91 0.96 1 0.95 0.91 0.88 0.84;...
                         0.09 0.23 0.32 0.48 0.67 0.79 0.87 0.95 1 0.98 0.97 0.94;...
                         0.06 0.18 0.26 0.42 0.60 0.73 0.81 0.91 0.98 1 0.99 0.97;...
                         0.04 0.15 0.24 0.38 0.57 0.70 0.78 0.88 0.97 0.99 1 0.99;...
                         0.05 0.13 0.21 0.33 0.53 0.66 0.74 0.84 0.94 0.97 0.99 1];
VRW                   = 0.16;
VTb                   = (2200*10^6)/FX_spot;
SF                    = 0.5.*min(1,14./(SIMM_tenors*365));
HVRir                 = 0.62;

% Bump G2++ sigma and eta parameters
shock_parameters      = parameters;
shock_parameters(1,2) = 0.0507036160404055;
shock_parameters(1,4) = 0.00874074464664423;


for j = 1:n_assets
    
    switch id_instruments(j)
        
        case 2 % swaption
            
            if evaluation_date < schedules{j,1}(1,1)
                
                % Compute bumped prices
                [V_shock,swap_rate_shock,swap_annuity_shock] = PortfolioExposure(id_instruments(j),ptf_data(j,:),schedules(j,:),DF_fwd_sim,DF_disc_sim,evaluation_date,libor_fixing_in,date,shock_parameters,swpt_exercise);
                                
                
                % Compute implied volatilities
                bs_price                             = V(:,j)./swap_annuity(:,j);
                bs_price_shock                       = V_shock./swap_annuity_shock;
                
                if ptf_data(j,6) == 1 % payer swaption
                    
                    impl_vol{j}                          = blsimpv(swap_rates(:,j)+displacement,strike(j)+displacement,0,swapt_expiry(j),bs_price,10, 0, [], true); %10: limit (default), 0: dividend yield(default), true -> call
                    impl_vol_shock{j}                    = blsimpv(swap_rate_shock+displacement,strike(j)+displacement,0,swapt_expiry(j),bs_price_shock,10, 0, [], true);
                    
                elseif ptf_data(j,6) == -1 % receiver swaption
                    
                    impl_vol{j}                          = blsimpv(swap_rates(:,j)+displacement,strike(j)+displacement,0,swapt_expiry(j),bs_price,10, 0, [], false); %10: limit (default), 0: dividend yield(default), false -> put
                    impl_vol_shock{j}                    = blsimpv(swap_rate_shock+displacement,strike(j)+displacement,0,swapt_expiry(j),bs_price_shock,10, 0, [], false);
                    
                end
                
                delta_impl_vol{j}                    = impl_vol_shock{j} - impl_vol{j};
                
                % Compute Vega
                delta_V(:,j)                         = V_shock - V(:,j);
                
                vega{j}                              = delta_V(:,j)./delta_impl_vol{j};
                
                vega{j}  = vega{j}.*check_ctp_IM{j};
                
                % Manage deltaV = delta_impl_vol = 0
                tmp                                  = vega{j};
                tmp(isnan(tmp))                      = 0;
                vega{j}                              = tmp;
                
                %%%
                
                % formula par 10.(c) ISDA SIMM: vega risk (before allocating on pillar SIMM)
                tmp_VR(:,j) = impl_vol{j}.*vega{j};
                
                % formula par 11.(a) ISDA SIMM: curvature risk (before allocating on pillar SIMM)
                tmp_CVR(:,j) = impl_vol{j}.*vega{j};
                
                % Allocating on pillar SIMM
                index = find(swapt_expiry(j) == SIMM_tenors);
                if isempty(index) == false % expiry = pillar SIMM
                    VR{index,j}  = tmp_VR;
                    CVR{index,j} = tmp_CVR.*SF(index);
                else % expiry != pillar SIMM --> proportional allocation
                    index_min = max(find(SIMM_tenors < swapt_expiry(j)));
                    index_max = min(find(SIMM_tenors > swapt_expiry(j)));
                    
                    if isempty(index_min) == true
                        VR{index_max,j}   = tmp_VR(:,j);
                        CVR{index_max,j}  = tmp_CVR(:,j).*SF(index_max);
                    elseif isempty(index_max) == true
                        VR{index_min,j}   = tmp_VR(:,j);
                        CVR{index_min,j}  = tmp_CVR(:,j).*SF(index_min);
                    else
                        
                        VR{index_min,j} = tmp_VR(:,j).* (SIMM_tenors(index_max) - swapt_expiry(j))/(SIMM_tenors(index_max) - SIMM_tenors(index_min));
                        VR{index_max,j} = tmp_VR(:,j).* (swapt_expiry(j) - SIMM_tenors(index_min))/(SIMM_tenors(index_max) - SIMM_tenors(index_min));
                        
                        CVR{index_min,j} = (tmp_CVR(:,j).*SF(index_min))* (SIMM_tenors(index_max) - swapt_expiry(j))/(SIMM_tenors(index_max) - SIMM_tenors(index_min));
                        CVR{index_max,j} = (tmp_CVR(:,j).*SF(index_max))* (swapt_expiry(j) - SIMM_tenors(index_min))/(SIMM_tenors(index_max) - SIMM_tenors(index_min));
                        
                    end
                    
                end
                
                
            else
                
                delta_V(:,j)        = zeros(n_sim,1);
                delta_impl_vol{j}   = zeros(n_sim,1);
                vega{j}             = zeros(n_sim,1);
                impl_vol{j}         = zeros(n_sim,1);
                impl_vol_shock{j}   = zeros(n_sim,1);
                
            end
            
        case 1 % irs
            
            delta_V(:,j)        = zeros(n_sim,1);
            delta_impl_vol{j}   = zeros(n_sim,1);
            vega{j}             = zeros(n_sim,1);
            impl_vol{j}         = zeros(n_sim,1);
            impl_vol_shock{j}   = zeros(n_sim,1);
            
        otherwise
            
            disp('not already implemented')
            
    end
    
end



% formula par 10.(d) ISDA SIMM, aggregation by instruments, get net sensitivity
tmp_VCRb         = zeros(n_sim,1);
tmp_sum_CVRb     = zeros(n_sim,1);
tmp_abs_sum_CVRb = zeros(n_sim,1);

for b = 1:length(SIMM_tenors)
    
    % vega
    tmp_net_VR{b}  = sum(cell2mat(VR(b,1:n_assets)),2);
    if isempty(tmp_net_VR{b}) == true
        tmp_net_VR{b} = zeros(n_sim,1);
    end
    tmp_VCRb  = tmp_VCRb + tmp_net_VR{b};
    
    % curvature
    tmp_net_CVR{b} = sum(cell2mat(CVR(b,1:n_assets)),2);
    if isempty(tmp_net_CVR{b}) == true
        tmp_net_CVR{b} = zeros(n_sim,1);
    end
    tmp_sum_CVRb     = tmp_sum_CVRb     + tmp_net_CVR{b};
    tmp_abs_sum_CVRb = tmp_abs_sum_CVRb + abs(tmp_net_CVR{b});
    
end

% vega
VCRb   = max(1, sqrt(abs(sum(tmp_VCRb,2))/VTb));
if verLessThan ('matlab','9.4')     == 0
    net_VR = VRW * cell2mat(tmp_net_VR) .* VCRb;
else
    net_VR = VRW * bsxfun(@times,cell2mat(tmp_net_VR),VCRb);
end

    % curvature
net_CVR = cell2mat(tmp_net_CVR);
theta   = min(sum(tmp_sum_CVRb,2)./sum(tmp_abs_sum_CVRb,2),0);
lambda  = (norminv(0.995)^2 - 1).*(1+theta)-theta;

% formula par 10.(e) ISDA SIMM
for i = 1:n_sim
    Kb(i) = sqrt(net_VR(i,:)*correlation_matrix*net_VR(i,:)');
    vega_margin(i) = Kb(i);

    % curvature
    Kb_curv(i) = sqrt(net_CVR(i,:)*correlation_matrix.^2*net_CVR(i,:)');

end
vega_margin      = vega_margin';
curvature_margin = max(sum(net_CVR,2)+lambda.*Kb_curv',0).*HVRir^-2;  
end

