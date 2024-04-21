function[dates,schedules] = find_dates(id_instruments, ptf_data,dv)
%{
% The function finds the relevant date for pricing the instruments in the portfolio:
% start date, expiry date (for swaption) and cash flow dates. 
% Furthermore schedules of the instruments in the portfolio are defined.
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
% dv               = vector containig holidays
%
%%% output:
%
% dates            = vector n_relevant_dates by 1 containing dates t_i 
%                    used for computing simulated ZCB(s,t_i) (where s is the pricing_date)
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
%}

n_assets                      = length(id_instruments);
start_date                    = ptf_data(:,1);
swaption_expiry               = ptf_data(:,2);
tenor_fix                     = ptf_data(:,7);
tenor_float                   = ptf_data(:,8);
dc_basis_fix                  = ptf_data(:,9);
dc_basis_float                = ptf_data(:,10);
schedules                     = cell(n_assets,3);
dates                         = [];

for i = 1:length(id_instruments)
    
    if id_instruments(i)      == 1
        maturity              = (ptf_data(i,3)-start_date(i))/365;
        sched                 = Swap_Date_Generator(start_date(i),maturity*12,tenor_fix(i),tenor_float(i),dv);
        date                  = sched(:,1);
        schedules{i,1}        = sched;
        initial_date          = sched(1,1); 
        fixed_leg_dates       = sched(sched(:,2)==1);
        fixed_leg_tau         = yearfrac([initial_date;fixed_leg_dates(1:end-1)],fixed_leg_dates,dc_basis_fix(i));
        schedules{i,2}        = [fixed_leg_dates zeros(length(fixed_leg_tau),1) fixed_leg_tau];
        floating_leg_dates    = sched(sched(:,3)==1);
        floating_leg_tau      = yearfrac([initial_date;floating_leg_dates(1:end-1)],floating_leg_dates,dc_basis_float(i));
        schedules{i,3}        = [floating_leg_dates zeros(length(floating_leg_tau),1) zeros(length(floating_leg_tau),1) floating_leg_tau];    
    
    elseif id_instruments(i) == 2
        start_date_swap       = busdate(busdate(swaption_expiry(i),1,dv),1,dv);    
        maturity              = (ptf_data(i,3)-swaption_expiry(i))/365;
        sched                 = Swap_Date_Generator(start_date_swap,maturity*12,tenor_fix(i),tenor_float(i),dv);
        schedules{i,4}        = sched;
        sched                 = [[swaption_expiry(i) 0 0];sched];
        date                  = sched(:,1); 
        schedules{i,1}        = sched;
        fixed_leg_dates       = sched(sched(:,2)==1);
        tmp_start_date        = sched(2,1);
        schedules{i,2}        = [tmp_start_date;fixed_leg_dates]; % input pricer swaption
        fixed_leg_tau         = yearfrac([tmp_start_date;fixed_leg_dates(1:end-1)],fixed_leg_dates,dc_basis_fix(i));
        schedules{i,5}        = [fixed_leg_dates zeros(length(fixed_leg_tau),1) fixed_leg_tau];
        floating_leg_dates    = sched(sched(:,3)==1);
        schedules{i,3}        = [sched(1,1);tmp_start_date; floating_leg_dates]; % input pricer swaption
        floating_leg_tau      = yearfrac([tmp_start_date;floating_leg_dates(1:end-1)],floating_leg_dates,dc_basis_float(i));
        schedules{i,6}        = [floating_leg_dates zeros(length(floating_leg_tau),1) zeros(length(floating_leg_tau),1) floating_leg_tau];
    
    end
    dates                     = [dates;date];
end
dates                         = unique(dates);
