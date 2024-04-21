function[exposure,swap_rate,annuity,floating_leg,fixed_leg] = Exposure(id_instrument,inst_data,schedule,DF_fwd_sim,DF_disc_sim,pricing_date,libor_fixing_in,date,parameters,swpt_exercise)
%
% The function produces a vector n_sim by 1 containig the MtM of a Swap or Swaption
% at a given pricing date
%
%%% input: 
%
% id_instrument     = id of the instrument: 
%                     1-swap
%                     2-swaption
%
% inst_data         = vector 1 by 13 with the following structure:
%                     [start_date, expiry_date (if swaption), tenor, strike, notional, omega
%                     (1->payer/-1->receiver), tenor_fixed_leg(mm),
%                     tenor_floating_leg(mm), dc_basis_fix, dc_basis_float,
%                     fixing_libor, next_cash_flow_date, tau_next_CF_date]
%
% schedule          = cell 1 by 3 for swap or 1 by 6 for swaption
%                     containing:
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
% DF_fwd_sim        = vector n_sim by n_date containig the discount factors for the forwarding curve  
%                     P(t,T_i) simulated for each path, where t is the point of the time grid  
%                     and e T_i are the relevant dates of the portfolio.
%              
% DF_disc_sim       = vector n_sim by n_date containig the discount factors for the discounting curve  
%                     P(t,T_i) simulated for each path, where t is the point of the time grid  
%                     and e T_i are the relevant dates of the portfolio.
%
% dv                = vector containig holidays 
%
% pricing_date      = time grid point corresponding to evaluation date 
%          
% libor_fixing_in   = cell vector 1 by n_fixings conaining matrices n_sim by 3 built as follows:
%                     col1 : LIBOR fixing 
%                     col2 : next LIBOR fixing date
%                     col3 : tau(current fixing date; next fixing date)
%
% date              = vector containing the relevant dates in correspondence of which discount factors have been simulated
%
% parameters        = matrix 9 by 5 containing:
%                     - row1 calibrated parameters
%                       col1 from 2:end expiries of market swaptions used for calibration
%                     - col2 from 2:end Gamma calibrated values  
%
% swpt_exercise     = vector n_paths by 1:
%                     'true'  -> swaption is exercised 
%                     'false' -> swaption not exercised
%
%%% output:
%
% exposure          = vector n_sim x 1 n_sim by 1 containig the MtM
%                     at time t related to each instrument for all
%                     simulated MC scenarios
%

n_sim                     = size(DF_fwd_sim,1);
exposure                  = zeros(n_sim,1);
libor_fixing_in_tmp       = cell2mat(libor_fixing_in);                     
swaption_expiry           = schedule{1}(1,1);

for i = 1:n_sim
    switch id_instrument
        case 1                                                             % Swap pricing
            [path_i_exp,swap_rate(i,:),annuity(i,:),floating_leg(i,:),fixed_leg(i,:)]          = Pricer(id_instrument,pricing_date,inst_data,schedule,DF_disc_sim(i,:),DF_fwd_sim(i,:),libor_fixing_in_tmp(i,:),date,parameters);
            exposure(i)                                       = path_i_exp;
            
        case 2                                                             % Swaption pricing
            if pricing_date < swaption_expiry
                floating_leg(i,:) = NaN;
                fixed_leg(i,:)    = NaN;
                [exposure(i,:),swap_rate(i,:),annuity(i,:)]   = Pricer(id_instrument,pricing_date,inst_data,schedule,DF_disc_sim(i,:),DF_fwd_sim(i,:),libor_fixing_in_tmp(i,:),date,parameters);
            else
                if swpt_exercise(i) == true                                % Physical swaption exercised = Swap
                    id_instrument_tmp                         = 1;
                    schedule_tmp                              = schedule(4:end);
                    [path_i_exp,swap_rate(i,:),annuity(i,:),floating_leg(i,:),fixed_leg(i,:)]  = Pricer(id_instrument_tmp,pricing_date,inst_data,schedule_tmp,DF_disc_sim(i,:),DF_fwd_sim(i,:),libor_fixing_in_tmp(i,:),date,parameters);
                    exposure(i)                               = path_i_exp;
                    clear id_instrument_tmp schedule_tmp
                else
                    swap_rate(i,:) = NaN;
                    annuity(i,:)   = NaN;
                end
            end
            
    end
end


end