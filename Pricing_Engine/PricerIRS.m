function [irs,swap_rate,annuity,floating_leg,fixed_leg,fwd_rate] = PricerIRS(pricing_date,inst_data,swap_schedule,DF_disc,DF_fwd,libor_fixing_in)
%{
% The function calculates Swap MtM
%
%%% input: 
%
% pricing_date      
%
% inst_data         = vector 1 by 13 with the following structure:
%                     [start_date, expiry_date (if swaption), tenor, strike, notional, omega
%                     (1->payer/-1->receiver), tenor_fixed_leg(mm),
%                     tenor_floating_leg(mm), dc_basis_fix, dc_basis_float,
%                     fixing_libor, next_cash_flow_date, tau_next_CF_date]
%
% swap_schedule     = cell 1 by 3 containing:
%                     {1} : matrix n_dates by 3. col1 : start_date, CF
%                           dates. col2: flag 1 of CF fixed leg, 0 otherwise.
%                           col3: flag 1 of CF foating leg, 0 otherwise
%                     {2} : matrix n_dates by 3 containig the schedule
%                           of the fixed leg: [fixed_leg_dates, zeros, fixed_leg_tau];
%                     {3} : matrix n_dates by 4 containing th schedule
%                           of the floating leg: [floating_leg_dates zeros, zeros, floating_leg_tau];
%
% DF_disc           = vector 1 by dates containing the Discount Factors (for discounting curve)
%                     corresponfing to the relevant dates of the instruments in the portfolio
%
% DF_fwd            = vector 1 by dates containing the Discount Factors (for forwarding curve)
%                     corresponfing to the relevant dates of the instruments in the portfolio
%
% libor_fixing_in   = matrix n_fixings by 3 containing data related to LIBOR fixing
%                     [fixing_libor, next_CF_date, next_CF_tau]. (if not available =[0,0,0])
%
%%% output:
%
% irs              = MtM of the Swap
%
% libor_fixing_out = if pricing date = LIBOR fixing date, the
%                    function returns LIBOR value, otherwise the function returns the input vector
%
% swap_rate        = Swap par-rate
%}

initial_date                        = swap_schedule{1}(1,1);
end_date                            = swap_schedule{1}(end,1);

if pricing_date < end_date
    K                               = inst_data(4);
    notional                        = inst_data(5);
    omega                           = inst_data(6);

    % tmp is a variabile used for building DFs vector with respect to the
    % residual cash flow dates (0 if CF date is before the pricing date)
    tmp                             = size(swap_schedule{1},1) - length(DF_disc);
    tmp1                            = zeros(1,tmp);
    tmp1                            = [tmp1 DF_disc];
    tmp2                            = zeros(1,tmp);
    tmp2                            = [tmp2 DF_fwd];
    DF_disc                         = [swap_schedule{1} tmp1'];
    floating_leg_DF_fwd             = [swap_schedule{1} tmp2'];

    % schedule fixed leg:
    % [CF_dates, DF_discCurve, tau(T_i-1;T_i)]
    fixed_leg_DF_disc               = DF_disc(DF_disc(:,2)==1,4);
    fixed_leg_schedule              = swap_schedule{2};
    fixed_leg_schedule(:,2)         = fixed_leg_DF_disc;

    % schedule floating leg:
    % [CF_dates, DF_discCurve, DF_fwdCurve, tau(T_i-1;T_i)]
    floating_leg_DF_disc            = DF_disc(DF_disc(:,3)==1,4);
    floating_leg_DF_fwd             = floating_leg_DF_fwd(floating_leg_DF_fwd(:,3)==1,4);
    floating_leg_schedule           = swap_schedule{3};
    floating_leg_schedule(:,2)      = floating_leg_DF_disc;
    floating_leg_schedule(:,3)      = floating_leg_DF_fwd;

    % Remove the cash flows before the pricing date
    actual_schedule                 = floating_leg_schedule(floating_leg_schedule(:,1)>=pricing_date);

    if pricing_date == actual_schedule(1)
        % if pricing_date = fixing date consider the corresponding elements of the schedule 
        fixed_leg_schedule          = fixed_leg_schedule(fixed_leg_schedule(:,1)>=pricing_date,:);
        floating_leg_schedule       = floating_leg_schedule(floating_leg_schedule(:,1)>=pricing_date,:);
    else
        % else consider the next elements
        fixed_leg_schedule          = fixed_leg_schedule(fixed_leg_schedule(:,1)>pricing_date,:);
        floating_leg_schedule       = floating_leg_schedule(floating_leg_schedule(:,1)>pricing_date,:);
    end

    %%% Pricing Swap

    % fixed leg
    fixed_leg_DF                    = fixed_leg_schedule(:,2);
    fixed_leg_tau                   = fixed_leg_schedule(:,3);
    fixed_leg                       = sum(K*fixed_leg_DF.*fixed_leg_tau)*notional;

    % floating leg
    floating_leg_dates              = floating_leg_schedule(:,1);
    floating_leg_DF_disc            = floating_leg_schedule(:,2);
    floating_leg_DF_fwd             = floating_leg_schedule(:,3);
    floating_leg_tau                = floating_leg_schedule(:,4);

    index                           = min(find(libor_fixing_in(:,2)>=pricing_date)); % float payment


    if pricing_date < initial_date  % forward start
        floating_leg_DF_fwd         = [DF_fwd(1);floating_leg_DF_fwd];                                 % add DF_fwdCurve to the start date to get LIBOR
        fwd_libor                   = (floating_leg_DF_fwd(1:end-1)./floating_leg_DF_fwd(2:end)-1)./floating_leg_tau;
        floating_leg                = notional*sum(fwd_libor.*floating_leg_tau.*floating_leg_DF_disc); % calculate floating leg

        fwd_rate = fwd_libor;

    else

        first_cf                    = libor_fixing_in(index,1)*libor_fixing_in(index,3)*floating_leg_DF_disc(1);

        if length(DF_fwd) > 1
            floating_leg_tau(1)     = [];
            fwd_libor               = (floating_leg_DF_fwd(1:end-1)./floating_leg_DF_fwd(2:end)-1)./floating_leg_tau;
            floating_leg            = sum(fwd_libor.*floating_leg_tau.*floating_leg_DF_disc(2:end));
            floating_leg            = notional*(floating_leg+first_cf); % calculate floating leg

            fwd_rate = [libor_fixing_in(index,1); fwd_libor];

        else
            floating_leg            = notional*(first_cf); % calculate floating leg 

            fwd_rate = libor_fixing_in(index,1);

        end

    end
    irs                                 = omega*(floating_leg-fixed_leg);
    annuity                             = sum(fixed_leg_DF.*fixed_leg_tau)*notional;
    swap_rate                           = floating_leg/annuity;
else
    irs                                 = 0;
    swap_rate                           = NaN;
    annuity                             = NaN;
    floating_leg                        = NaN;
    fixed_leg                           = NaN;
    fwd_rate                            = NaN;
end
end