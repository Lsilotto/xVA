function [price,swap_rate,annuity,floating_leg,fixed_leg] = Pricer(id_instrument,pricing_date,inst_data,schedule,DF_disc,DF_fwd,libor_fixing_in,date,parameters)
%{
% The function calculates Swap and Swaption MtM
%
%%% input: 
%
% id_instrument     = id of the instrument: 
%                     1-swap
%                     2-swaption
%
% pricing_date 
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
% DF_disc           = vector 1 by dates containing the Discount Factors (for discounting curve)
%                     corresponfing to the relevant dates of the instruments in the portfolio
%
% DF_fwd            = vector 1 by dates containing the Discount Factors (for forwarding curve)
%                     corresponfing to the relevant dates of the instruments in the portfolio
%
% libor_fixing_in   = vector 1 by (3*n_fixings) containing the simulated LIBOR fixings                    
%
% date              = relevant dates of the instruments in the portfolio
%                     (start_date, expiry(if swaption), CF_dates) in Matlab format
%
% parameters        = matrix 9 by 5 containing:
%                     - row1 calibrated parameters
%                       col1 from 2:end expiries of market swaptions used for calibration
%                     - col2 from 2:end Gamma calibrated values                    
%
%%% output:
%
% price             = MtM of the instrument
%
% libor_fixing_out  = (for IRS) if pricing date = LIBOR fixing date, the
%                     function returns LIBOR value, otherwise the function returns the input vector
%}

if id_instrument                 == 1
    swap_schedule                = schedule{1};
    tmp                          = [date DF_disc'];                        % find DFs corresponding to instrument dates
    tmp1                         = [date DF_fwd'];
    k                            = 1;
    for i = 1:length(swap_schedule(:,1))
        index                    = find(swap_schedule(i,1) == tmp(:,1));
        if isempty(index)        == 0
            DF_disc_in(k)        = tmp(index,2);
            DF_fwd_in(k)         = tmp1(index,2);
            k                    = k+1; 
        end
    end
    if k == 1
            DF_disc_in(1)        = 0;
            DF_fwd_in(1)         = 0;
    end
    clear tmp tmp1
    k = 0;
    % Data reordering
    for i = 1:(length(libor_fixing_in)/3)
        libor_fixing_in_tmp(i,:) = libor_fixing_in(k+1:k+3);
        k = k+3;
    end
    [price,swap_rate,annuity,floating_leg,fixed_leg]    = PricerIRS(pricing_date,inst_data,schedule,DF_disc_in,DF_fwd_in,libor_fixing_in_tmp);

elseif id_instrument             == 2 
    float_schedule               = schedule{3};
    PayRec                       = 0; 
    Method                       = 1;
    dcbasis_Fix                  = inst_data(9);
    Strike                       = inst_data(4);
    tmp                          = [date DF_disc'];                        % find DFs corresponding to instrument dates
    tmp1                         = [date DF_fwd'];
    k                            = 1;
    for i = 1:length(float_schedule(:,1))
        index                    = find(float_schedule(i,1) == tmp(:,1));
        if isempty(index)    == 0
            DF_disc_in(k)        = tmp(index,2);
            DF_fwd_in(k)         = tmp1(index,2);
            k                    = k+1;
        end
    end
    clear tmp tmp1
    float_schedule               = float_schedule(2:end);
    [Price_P,Price_R]            = Price_Swaption_HW2F_Gamma_MULTI(pricing_date,inst_data(2),parameters,schedule{2},float_schedule,DF_disc_in',DF_fwd_in',Strike,dcbasis_Fix,PayRec,Method);

    switch inst_data(6)
        case 1 
            price                = Price_P*inst_data(5);                   % payer swaption
        case -1 
            price                = Price_R*inst_data(5);                   % receiver swaption
    end
    
    clear DF_disc_in DF_fwd_in
    
    tmp_schedule                 = schedule(4:end);
    DF_disc                      = DF_disc(2:end);
    DF_fwd                       = DF_fwd(2:end);
    swap_schedule                = tmp_schedule{1};
    tmp                          = [date(2:end) DF_disc'];                 % find DFs corresponding to instrument dates
    tmp1                         = [date(2:end) DF_fwd'];
    k                            = 1;
    for i = 1:length(swap_schedule(:,1))
        index                    = find(swap_schedule(i,1) == tmp(:,1));
        if isempty(index)        == 0
            DF_disc_in(k)        = tmp(index,2);
            DF_fwd_in(k)         = tmp1(index,2);
            k                    = k+1; 
        end
    end
    if k == 1
            DF_disc_in(1)        = 0;
            DF_fwd_in(1)         = 0;
    end
    clear tmp tmp1
    k = 0;
    % Data reordering
    for i = 1:(length(libor_fixing_in)/3)
        libor_fixing_in_tmp(i,:) = libor_fixing_in(k+1:k+3);
        k = k+3;
    end
    [~,swap_rate,annuity]        = PricerIRS(pricing_date,inst_data,tmp_schedule,DF_disc_in,DF_fwd_in,libor_fixing_in_tmp);
    
    
end


end