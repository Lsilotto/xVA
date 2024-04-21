function [swpt_exercise,mtm_expiry] = find_swaption_exercise(id_instruments,ptf_data,x,y,xy_time_grid,pricing_date,r_disc,DF_disc_tau,r_fwd,DF_fwd_tau,interp_method,schedules,libor_fixings,parameters)
%{
% The function finds if the swaption is exercised for the various simulated MC scenarios by
% simulating the MtM of the underlying Swap at Swaption expiry date. If MtM > 0 -> Swaption is exercised
%
%%% input:
%
% id_instruments   = vector n_assets by 1 containing the id of the instrument: 
%                    1-swap
%                    2-swaption
%
% ptf_data         = matrix n_assets by 13 containing on the columns the following 
%                     [start_date, expiry_date (if swaption), tenor, strike, notional, omega
%                     (1->payer/-1->receiver), tenor_fixed_leg(mm),
%                     tenor_floating_leg(mm), dc_basis_fix, dc_basis_float,
%                     fixing_libor, next_cash_flow_date, tau_next_CF_date]
%
% x                = matrix n_path by length(xy_time_grid) -> simulated G2++ x
%                    process
%
% y                = matrix n_path by length(xy_time_grid) -> simulated G2++ y
%                    process
%
% xy_time_grid     = x and y processes daily simulation time grid
%
% pricing_date     = today (Matlb date format)
%
% r_disc           = vector n_pillar_DF_curve_input by 1 containing the input zero
%                    rates (discounting curve)
%
% DF_disc_tau      = vector n_pillar_DF_curve_input by 1 containing the year fraction
%                    for the input curve (discounting)
%
% r_fwd            = vector n_pillar_DF_curve_input by 1 containing the input zero
%                    rates (forwarding curve)
%
% DF_fwd_tau       = vector n_pillar_DF_curve_input by 1 containing the year fraction
%                    for the input curve (forwarding)
%
% interp_method    = chosen interpolation method (see Matlab documentation)
%
% schedules        = cell matrices n_assets by 3 for swap or by 6 for swaption
%                    containing:
%                     - IRS
%                         {1} : matrix n_dates by 3. col1 : start_date, CF
%                               dates. col2: flag 1 of CF fixed leg, 0 otherwise.
%                               col3: flag 1 of CF foating leg, 0 otherwise
%                         {2} : matrix n_dates by 3 containig the schedule
%                               of the fixed leg: [fixed_leg_dates, zeros, fixed_leg_tau];
%                         {3} : matrix n_dates by 4 containing th schedule
%                               of the floating leg: [floating_leg_dates zeros, zeros, floating_leg_tau];
%
%                     - Swaption
%                         {1} : matrix n_dates by 3. col1 : expiry_date, start_date, CF
%                               dates. col2: flag 1 if CF fixed leg, 0 otherwise.
%                               col3: flag 1 if CF foating leg, 0 otherwise
%                         {2} : vector n_dates by 1 containing the CF dates
%                               of the fixed leg: start_date,start_date,CF_dates
%                         {3} : vector n_dates by 1 containing the schedule
%                               of the floating leg: expiry_date, start_date, CF_dates
%                         From {3} to {6} schedule of the underlying swap as described above
%
% libor_fixings    = cell matrix n_assets by n_fixing(max) containing matrices of
%                    dimension n_sim by 3 as follows:
%                    col1 : LIBOR fixing 
%                    col2 : next LIBOR fixing date
%                    col3 : tau(current fixing date; next fixing date)
%
% parameters       = matrix 9 by 5 containing:
%                    - row1 calibrated parameters
%                      col1 from 2:end expiries of market swaptions used for calibration
%                    - col2 from 2:end Gamma calibrated values 
%
%%% output:
%
% swpt_exercise     = cell n_assets by 1 containing vectors n_paths by 1:
%                    'true'  -> swaption is exercised 
%                    'false' -> swaption not exercised
%
% mtm_expiry        = cell n_assets by 1 containing vectors n_paths by 1:
%                     MtM of the underlying Swap at Swaption expiry date
%}

n_assets                   = length(id_instruments);
n_paths                    = size(x,1);
swaption_exercise_date     = ptf_data(:,2); 
swaption_exercise_date_yf  = (swaption_exercise_date - pricing_date)/365;
swpt_exercise              = cell(n_assets,1);
mtm_expiry                 = cell(n_assets,1);
tmp_swpt_exercise          = true(n_paths,1);

for j = 1:n_assets
    
    if id_instruments(j) == 2
        % find x and y at Expiry date
        if verLessThan ('matlab','9.4') == 0
            index         = find(round(swaption_exercise_date_yf(j),4) == round(xy_time_grid,4));
        else
            index         = find(round(swaption_exercise_date_yf(j)*10000)/10000 == round(xy_time_grid*10000)/10000); % workaround per versioni matlab < 2018a 
        end
        
        % forwarding and discounting curves construction
        instrument_CF_dates = schedules{j,1}(:,1);
        CF_dates            = instrument_CF_dates(instrument_CF_dates >= swaption_exercise_date(j)); % date CF formato data matlab
        curve_pillars       = (CF_dates - pricing_date)/365;
        
        CF_disc_curve = ZCB_tT(r_disc,swaption_exercise_date_yf(j),curve_pillars',DF_disc_tau,interp_method,parameters,x(:,index),y(:,index));
        CF_fwd_curve  = ZCB_tT(r_fwd, swaption_exercise_date_yf(j),curve_pillars',DF_fwd_tau ,interp_method,parameters,x(:,index),y(:,index));
        
        % calculation of MtM of the underlying Swap at Swaption expiry date
        mtm_expiry{j}    = Exposure(id_instruments(j),ptf_data(j,:),schedules(j,:),CF_fwd_curve,CF_disc_curve,swaption_exercise_date(j),libor_fixings(j,:),CF_dates,parameters,tmp_swpt_exercise);
        
        mtm_expiry{j} = mtm_expiry{j} * ptf_data(j,5);
        
        for n = 1:n_paths
            if mtm_expiry{j}(n) > 0
                swpt_exercise{j}(n) = true;
            else
                swpt_exercise{j}(n) = false;
            end
        end
        swpt_exercise{j} = swpt_exercise{j}';
    end
end
        
            
            
        
        
    
        
        

