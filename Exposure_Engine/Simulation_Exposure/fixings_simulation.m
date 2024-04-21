function [simulated_libor_fixings] = fixings_simulation(ptf_data,schedules,pricing_date,time_grid,parameters,x,y,DF_fwd_mkt,id_instruments,interp_method)
%{ 
% The function computes LIBOR fixings from G2++ x and y processes simulated on a daily grid for Swaps and Swaptions
%
%%% input:
%
% ptf_data         = matrix n_assets by 13 containing on the columns the following 
%                     [start_date, expiry_date (if swaption), tenor, strike, notional, omega
%                     (1->payer/-1->receiver), tenor_fixed_leg(mm),
%                     tenor_floating_leg(mm), dc_basis_fix, dc_basis_float,
%                     fixing_libor, next_cash_flow_date, tau_next_CF_date]
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
% pricing_date
%
% time_grid        = time grid used for simulating the exposure  (in year fraction) 
%
% parameters       = matrix 9 by 5 containing:
%                    - row1 calibrated parameters
%                      col1 from 2:end expiries of market swaptions used for calibration
%                    - col2 from 2:end Gamma calibrated values
%
% x                = matrix n_path by length(xy_time_grid) -> simulated G2++ x
%                    process
%
% y                = matrix n_path by length(xy_time_grid) -> simulated G2++ y
%                    process
%
% DF_fwd_mkt       = market term strucuture of discount factors for forwarding curve
%
% id_instruments   = vector n_assets by 1 containing the id of the instrument: 
%                    1-swap
%                    2-swaption
%
%%% output:
%
% simulated_libor_fixings = cell matrix n_assets by n_fixing(max) containing matrices of
%                           dimension n_sim by 3 as follows:
%                           col1 : LIBOR fixing 
%                           col2 : next LIBOR fixing date
%                           col3 : tau(current fixing date; next fixing date)
%}                           
                          
n_sim                                       = size(x,1);
n_assets                                    = size(ptf_data,1);
tau     = DF_fwd_mkt(2:end,1)/365;
tmp_Rf  = - log(DF_fwd_mkt(2:end,2))./tau;
Rf      = [0;tmp_Rf];
interp_DF_fwd_mkt                           = exp(-interp1(DF_fwd_mkt(:,1)/365,Rf,time_grid,interp_method,'extrap').*time_grid);

tmp_time_grid                               = pricing_date+round(time_grid*365);

for j = 1:n_assets
    if id_instruments(j)     == 1
        start_date                          = schedules{j,1}(1,1);
        tmp_schedule                        = [start_date 0 0 0; schedules{j,3}];
    elseif id_instruments(j) == 2
        start_date                          = schedules{j,1}(2,1);
        tmp_schedule                        = [start_date 0 0 0; schedules{j,6}];
    end
    for i = 1:size(tmp_schedule,1)
        index(i)                            = find(tmp_schedule(i,1) == tmp_time_grid);
    end 
    % Compute simulated discount factors
    fixing_grid                             = time_grid(index);
    for t = 2:length(fixing_grid)
        DF_fwd_mkt_0T                       = interp_DF_fwd_mkt(index(t));
        DF_fwd_mkt_0t                       = interp_DF_fwd_mkt(index(t-1));
        DF_fwd_simulated(:,t-1)             = ZCB_HW2F_Gamma(fixing_grid(t-1),fixing_grid(t),parameters,DF_fwd_mkt_0T,...
                                              DF_fwd_mkt_0t,x(:,index(t-1)),y(:,index(t-1)));
    end
    % Compute simulated LIBOR fixings
    for i = 1:size(DF_fwd_simulated,2)
        simulated_libor_fixings{j,i}(:,1)   = (1./DF_fwd_simulated(:,i)-1)./tmp_schedule(i+1,end);
        
        simulated_libor_fixings{j,i}(:,2)   = tmp_schedule(i+1,1);                                % next fixing date
        simulated_libor_fixings{j,i}(:,3)   = tmp_schedule(i+1,end);                              % tau(T(i-1);T(i))
    end
    % If the firs LIBOR fixing is already known, read from input file 
    if ptf_data(j,11) ~= 0
        simulated_libor_fixings{j,1}        = repmat(ptf_data(i,11:end),n_sim,1);
    end
    clear index tmp_libor_fixing DF_fwd_simulated
end
