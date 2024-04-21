function [ptfExposure,swap_rate,annuity,floating_leg,fixed_leg] = PortfolioExposure(id_instruments,ptf_data,schedules,DF_fwd_sim,DF_disc_sim,pricing_date,libor_fixing_in,date,parameters,swpt_exercise)
%{
% The function produces the exposure for a portfolio of Swaps and Swaptions 
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
%
% DF_fwd_sim       = vector n_sim by n_date containig the discount factors for the forwarding curve  
%                    P(t,T_i) simulated for each path, where t is the point of the time grid  
%                    and e T_i are the relevant dates of the portfolio.
%              
% DF_disc_sim      = vector n_sim by n_date containig the discount factors for the discounting curve  
%                    P(t,T_i) simulated for each path, where t is the point of the time grid  
%                    and e T_i are the relevant dates of the portfolio.
%
% pricing_date     = time grid point corresponding to evaluation date 
%                    la exposure
%          
% libor_fixing_in  = cell vector 1 by n_fixings conaining matrices n_sim by 3 built as follows:
%                    col1 : LIBOR fixing 
%                    col2 : next LIBOR fixing date
%                    col3 : tau(current fixing date; next fixing date)
%
% date             = vector containing the relevant dates in correspondence of which discount factors have been simulated
%
% parameters       = matrix 9 by 5 containing:
%                    - row1 calibrated parameters
%                      col1 from 2:end expiries of market swaptions used for calibration
%                    - col2 from 2:end Gamma calibrated values  
%
% swpt_exercise    = cell vector n_assets by 1 containing vectors n_paths by 1:
%                    'true'  -> swaption is exercised 
%                    'false' -> swaption not exercised
% 
%%% output:
%
% ptfExposure      = vector n_sim by n_assets containing exposure values at
%                    time t for each instrument in all simulated MC scenarios
%
%%% assumptions:   - floating cash flow dates correspond to LIBOR fixing dates
%                  - if no LIBOR fixing values in input and t corresponds to start date 
%                    (or first evaluation date in case of swap forward
%                    starting) --> they are computed in the function
%}

n_sim                        = size(DF_disc_sim,1);
if nargin < 10
    swpt_exercise = zeros(length(n_sim),1);
end
n_assets                     = length(id_instruments);
ptfExposure                  = zeros(n_sim,n_assets);
swap_rate                    = zeros(n_sim,n_assets);
annuity                      = zeros(n_sim,n_assets);

for i = 1:n_assets
    [ptfExposure(:,i), swap_rate(:,i), annuity(:,i), floating_leg(:,i), fixed_leg(:,i)]    = Exposure(id_instruments(i),ptf_data(i,:),schedules(i,:),DF_fwd_sim,DF_disc_sim,pricing_date,libor_fixing_in(i,:),date,parameters,swpt_exercise{i});
end



