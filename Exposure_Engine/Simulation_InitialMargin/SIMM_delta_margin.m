function[delta_margin,DF_disc_sensitivity,DF_fwd_sensitivity] = SIMM_delta_margin(shock_DF_disc,shock_DF_fwd,DF_disc,DF_fwd,V,id_instruments,ptf_data,schedule,evaluation_date,pricing_date,simulated_libor_fixings,ptf_dates,parameters,swpt_exercise,FX_spot,Jacobian_disc,Jacobian_fwd_disc,Jacobian_fwd,fwd_curve_pillars,disc_curve_pillars)

n_assets               = length(id_instruments);
DF_disc_sensitivity    = cell(1,length(id_instruments));
DF_fwd_sensitivity     = cell(1,length(id_instruments));
n_path                 = size(DF_fwd,1);


% Compute sensitivity for instruments in the portfolio
if nargin > 15 % use Jacobian

    for j = 1:n_assets
        [DF_disc_sensitivity{j}, DF_fwd_sensitivity{j}] = SIMM_delta_Jacobian(shock_DF_disc,shock_DF_fwd,DF_disc,DF_fwd,V{j},id_instruments(j),ptf_data(j,:),schedule(j,:),evaluation_date,pricing_date,simulated_libor_fixings(j,:),ptf_dates,parameters,swpt_exercise(j),Jacobian_disc,Jacobian_fwd_disc,Jacobian_fwd,fwd_curve_pillars,disc_curve_pillars);
    end
    
else

    for j = 1:n_assets
        [DF_disc_sensitivity{j}, DF_fwd_sensitivity{j}] = SIMM_delta(shock_DF_disc,shock_DF_fwd,DF_disc,DF_fwd,V{j},id_instruments(j),ptf_data(j,:),schedule(j,:),evaluation_date,simulated_libor_fixings(j,:),ptf_dates,parameters,swpt_exercise(j));
    end
    
end


% Aggregate sensitivity for instruments in the portfolio (get s_k,i)               
% output -> matrx n_path by n_pillarSIMM containing net sensitivities

SIMM_pillars           = size(DF_disc_sensitivity{1},2);
% Discounting curve
tmp_s_disc_curve       = zeros(n_path,SIMM_pillars);
for j = 1:n_assets
    tmp_s_disc_curve   = tmp_s_disc_curve + DF_disc_sensitivity{j};
end
s_disc_curve           = tmp_s_disc_curve;

% Forwarding curve
tmp_s_fwd_curve        = zeros(n_path,SIMM_pillars);
for j = 1:n_assets
    tmp_s_fwd_curve    = tmp_s_fwd_curve + DF_fwd_sensitivity{j};
end
s_fwd_curve            = tmp_s_fwd_curve;


% Define ISDA-SIMM parameters
%+++++++++++ V 2.1 
RW                              = [114 115 102 71 61 52 50 51 51 51 54 62];  % Risk Weights
Tb                              = (210*10^6)/FX_spot;                        % Concentration threshold
correlation_matrix              = [1 0.63 0.59 0.47 0.31 0.22 0.18 0.14 0.09 0.06 0.04 0.05;...
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
sub_curves_correlation          = 0.98;


% Compute weighted sensitivity
% output -> matrix n_path by n_pillarSIMM containing weighted sensitivities
CR                                  = max(1,sqrt(abs(sum(s_disc_curve,2)+sum(s_fwd_curve,2))/Tb)); % Concentration Risk Factor
if verLessThan ('matlab','9.4')     == 0
    WS_fwd_curve                    = RW.*s_fwd_curve.*CR;
    WS_disc_curve                   = RW.*s_disc_curve.*CR;
else
    WS_curve                        = bsxfun(@times,RW,CR);
    WS_fwd_curve                    = bsxfun(@times,WS_curve,s_fwd_curve);
    WS_disc_curve                   = bsxfun(@times,WS_curve,s_disc_curve);
end


% Aggragate weighted sensitivities
for i = 1:n_path
    
    a = WS_fwd_curve(i,:)*correlation_matrix*WS_fwd_curve(i,:)';
    b = WS_disc_curve(i,:)*correlation_matrix*WS_disc_curve(i,:)';
    c = 2*(WS_fwd_curve(i,:)*correlation_matrix*WS_disc_curve(i,:)'*sub_curves_correlation);
    
    delta_margin(i,:) = sqrt(a+b+c);
    
end

end