% This code simulates the exposure of Swaps and Swaptions (phisically settled) and calculates CVA, DVA, FCA, FBA and MVA 
% considering both Variation Margin and Initial Margin (according to ISDA-SIMM).
% This code works also for a portfolio of Swaps/Swaptions

%% Import market data from Excel file
[d1,d2,d3]                              = textread('TargetDates.txt','%s %s %s');     % holiday calendar
dv                                      = strcat(d1,d2,d3);
dv                                      = datenum(dv);
Date                                    = '28_12_2018';                               % market data from excel file
filename                                = ['MarketData_Input_' Date '.xls'];
DF_disc_input                           = xlsread(filename,'Discounting');            % discounting curve
DF_disc_input(:,1)                      = x2mdate(DF_disc_input(:,1));
DF_fwd_input                            = xlsread(filename,'Forwarding');             % forwarding curve
DF_fwd_input(:,1)                       = x2mdate(DF_fwd_input(:,1));

pricing_date                            = DF_disc_input(1,1);                                    % set pricing date
DF_disc                                 = [DF_disc_input(:,1)-pricing_date DF_disc_input(:,2)];  % discounting curve Discount Factors 
r_disc                                  = [0;-log(DF_disc(2:end,2))./(DF_disc(2:end,1)/365)]';   % discounting curve Zero Rates
DF_fwd                                  = [DF_fwd_input(:,1)-pricing_date DF_fwd_input(:,2)];    % forwarding curve Discount Factors for f
r_fwd                                   = [0;-log(DF_fwd(2:end,2))./(DF_fwd(2:end,1)/365)]';     % forwarding curve Zero Rates

clear('d1','d2','d3','Date','filename')

%% Set simulation parameters
flag_IM            = 1;         % set the flag equal to 1 for simulating Initial Margin, 0 otherwise
flag_VM            = 1;         % set the flag equal to 1 for simulating Variation Margin, 0 otherwise
n_sim              = 2;      % set the number of Monte Carlo scenarios (notice that the MC simulation considers antithetic variates, set this parameter to even number)
flag_seed          = 1;         % set the flag equal to 1 for fixing the seed of MC simulation, 0 otherwise
grid_basis         = 'month';   % set the time grid basis ('day','month')
grid_interval      = 12;         % set the interval between the points of the time grid
grid_type          = 'joint';   % set the parameter to 'joint' in order to add points on the time grid corresponding to the cash flows of the instrument, 'std' for an equally spaced time grid  
interp_method      = 'linear';  % set the interpolation method (see Matlab documentation)

%% Set CSA parameters
if (flag_VM == 1 || flag_IM == 1)
    MPOR_dd        = 2;   % Margin Period of Risk in days
end
if flag_VM == 1
    csa_threshold  = 0;   % CSA Threshold
    csa_mta        = 0;   % CSA Minimum Transfer Amount
end

%% Import instrument data from Excel file
ptf_data                                = xlsread('portfolio_28_12_18','Load_instrument');
id_instruments                          = ptf_data(:,1);                                           % 1 for Swap, 2 for Swaption
ptf_data(:,2)                           = x2mdate(ptf_data(:,2));                                  % Instrument start date
ptf_data(ptf_data(:,1)==2,3)            = busdate(x2mdate(ptf_data(ptf_data(:,1)==2,3))-1,1,dv);   % Instrument expiry date (for swpation only)
ptf_data(:,4)                           = x2mdate(ptf_data(:,4));                                  % Instrument end date
ptf_data(:,13)                          = x2mdate(ptf_data(:,4));                                  % Instrument next cash flow date
ptf_data                                = ptf_data(:,2:end);                                        
[instruments_dates, schedules]          = find_dates(id_instruments,ptf_data,dv);                  % Get the schedule of the instrument in order to calculate the corresponding Zero Coupon Bond
T_yy                                    = (instruments_dates(end) - pricing_date)/365;             % Instrument maturity in years

clear('dv')

%% G2++ parameters calibration
% load G2++ parameters calibrated within the scripts:
%   - 'main_calibration' for parameters calibrated on ATM volatilities;
%   - 'cube_calibration' for parameters calibrated on volatility cube;

load('28_12_2018_Pars_HW2F_Gamma_new_ATM_fullExpiries.mat');
parameters = Pars_HW2F_Gamma_new;

%% Simulation of G2++ x and y processes and LIBOR fixings
n_step                                 = round(365*T_yy);
dt                                     = T_yy/n_step;
time_grid                              = (0:dt:T_yy);
[~,~,x,y]                              = HW2F_xy_simulation(parameters,DF_disc(2:end,:),DF_fwd(2:end,:),T_yy,n_step,n_sim,time_grid,flag_seed);
simulated_libor_fixings                = fixings_simulation(ptf_data,schedules,pricing_date,time_grid,parameters,x,y,DF_fwd,id_instruments,interp_method);

clear('n_step','dt')

%% Exposure and Collateral time grids construction
% Standard grid: equally spaced time grid 
% Joint grid: add instrument cash flow dates + 1 day to the Standard grid

% Build Exposure time grid
if strcmp(grid_type,'joint') == 1
    tmp_instruments_CF_dates                               = instruments_dates(2:end-1)';
    [exp_evaluation_dates,exposure_grid_yf,x_exp,y_exp]    = build_exposure_grid(pricing_date,instruments_dates(end),grid_basis,grid_interval,x,y,time_grid,tmp_instruments_CF_dates);
else
    [exp_evaluation_dates,exposure_grid_yf,x_exp,y_exp]    = build_exposure_grid(pricing_date,instruments_dates(end),grid_basis,grid_interval,x,y,time_grid);
end
% Interpolating on market ZCB to get discount factors for the simulated Mark-to-Futures
DF_disc_curve_mkt_exp_grid                                 = exp(-interp1(DF_disc(:,1)/365,r_disc,exposure_grid_yf).*exposure_grid_yf);

% Bulding collateral grid
if flag_VM == 1 || flag_IM == 1
    [col_evaluation_dates,collateral_grid_yf,x_coll,y_coll] = build_collateral_grid(exp_evaluation_dates,MPOR_dd,x,y,time_grid);
    DF_disc_curve_mkt_coll_grid                             = exp(-interp1(DF_disc(:,1)/365,r_disc,collateral_grid_yf).*collateral_grid_yf);
end

%% Find MC scenarios in which the Swaption is excercised 
if any(id_instruments == 2) == false
    swpt_exercise = cell(length(id_instruments),1);
else
    [swpt_exercise,swap_mtm_expiry] = find_swaption_exercise(id_instruments,ptf_data,x,y,time_grid,pricing_date,r_disc,DF_disc(:,1)/365,r_fwd,DF_fwd(:,1)/365,interp_method,schedules,simulated_libor_fixings,parameters);
end

clear('time_grid','x','y')

%% Simulation of Mark-to-Futures, Variation Margin and Initial Margin (ISDA-SIMM)
n_assets             = length(id_instruments);
instruments_dates_yf = (instruments_dates-pricing_date)/365;
SIMM_pillars         = [14 30 91 182 365 365*2 365*3 365*5 3650 365*15 365*20 365*30]./365; % Definition of ISDA-SIMM tenors

% Initialize variables
ptf_uncollateralized_exposure = cell(n_assets,length(exposure_grid_yf));

if flag_VM == 1
    ptf_variation_margin      = cell(n_assets,length(collateral_grid_yf));
    ptf_exposure_tMPOR        = cell(n_assets,length(collateral_grid_yf));
end

if flag_IM ==1
    SIMM_mtf                  = cell(n_assets,length(collateral_grid_yf));
    delta_margin              = zeros(n_sim,length(collateral_grid_yf));
    vega_margin               = zeros(n_sim,length(collateral_grid_yf));
    curvature_margin          = zeros(n_sim,length(collateral_grid_yf));
    initial_margin            = zeros(n_sim,length(collateral_grid_yf));
    displacements             = 0.06;
    FX_spot                   = 1.145;
end


% Loop on time grid in order to get MtF(t), VM(t) and IM(t)
tic
for t = 1:length(exposure_grid_yf)
    disp('Exposure Simulation'); disp(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculation of MtF at time t for each MC scenario %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Forwarding and discounting curves construction
    residual_CF_dates_yf    = instruments_dates_yf(instruments_dates_yf >= exposure_grid_yf(t));
    CF_disc_curve_exp_grid  = ZCB_tT(r_disc,exposure_grid_yf(t),residual_CF_dates_yf',DF_disc(:,1)/365,interp_method,parameters,x_exp(:,t),y_exp(:,t));
    CF_fwd_curve_exp_grid   = ZCB_tT(r_fwd, exposure_grid_yf(t),residual_CF_dates_yf',DF_fwd(:,1)/365, interp_method,parameters,x_exp(:,t),y_exp(:,t));
    
    % Get MtF at time t
    residual_CF_dates       = instruments_dates(instruments_dates >= exp_evaluation_dates(t));
    [und_ptf_exposure, swap_rate(:,t)]...
                            = PortfolioExposure(id_instruments,ptf_data,schedules,CF_fwd_curve_exp_grid,CF_disc_curve_exp_grid,exp_evaluation_dates(t),simulated_libor_fixings,residual_CF_dates,parameters,swpt_exercise);  % tmp_ptf_exposure: matrice n_sim x n_assets
    tmp_ptf_exposure        = DF_disc_curve_mkt_exp_grid(t) * und_ptf_exposure;
    
    % Results reordered in a cell variable (n_assets by number of time grid points)
    % containing for each time step a vector (n_MC_senarios by 1) containig the discounted MtF for each simulated path
    for i = 1:n_assets
        ptf_uncollateralized_exposure{i,t} = tmp_ptf_exposure(:,i);
    end
    clear('und_ptf_exposure','tmp_ptf_exposure','residual_CF_dates_yf','CF_disc_curve_exp_grid','CF_fwd_curve_exp_grid','residual_CF_dates')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculation of VM at time t for each MC scenario %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag_VM == 1
        
        if collateral_grid_yf(t) < 0 % VM(0) = 0
            
            for i = 1:n_assets
                ptf_variation_margin{i,t} = zeros(n_sim,1);
            end
            
        else
            
            % Forwarding and discounting curves construction
            residual_CF_dates_yf_coll     = instruments_dates_yf(instruments_dates_yf >= collateral_grid_yf(t));
            CF_disc_curve_coll_grid       = ZCB_tT(r_disc,collateral_grid_yf(t),residual_CF_dates_yf_coll',DF_disc(:,1)/365,interp_method,parameters,x_coll(:,t),y_coll(:,t));
            CF_fwd_curve_coll_grid        = ZCB_tT(r_fwd, collateral_grid_yf(t),residual_CF_dates_yf_coll',DF_fwd(:,1)/365, interp_method,parameters,x_coll(:,t),y_coll(:,t));
            
            % Get MtF at time t-tMPOR
            residual_CF_dates_coll        = instruments_dates(instruments_dates >= col_evaluation_dates(t));
            tmp_ptf_exposure_tMPOR        = PortfolioExposure(id_instruments,ptf_data,schedules,CF_fwd_curve_coll_grid,CF_disc_curve_coll_grid,col_evaluation_dates(t),simulated_libor_fixings,residual_CF_dates_coll,parameters,swpt_exercise);
            tmp_ptf_exposure_tMPOR        = DF_disc_curve_mkt_coll_grid(t) * tmp_ptf_exposure_tMPOR;
            
            
            if MPOR_dd == 0 && t == 1
                
                for i = 1:n_assets
                    tmp_ptf_variation_margin{i,1} = zeros(n_sim,1);
                end
                % Capitalization factor for VM(t-1;t)
                coll_cptl_disc_curve          = ones(n_sim,1);
                
                for i = 1:n_assets
                    % Results reordered
                    ptf_exposure_tMPOR{i,t}   = tmp_ptf_exposure_tMPOR(:,i);
                    % VM calculation
                    ptf_variation_margin{i,t} = collateral_simulation(csa_threshold,csa_mta,ptf_exposure_tMPOR{i,t},tmp_ptf_variation_margin{:,1},coll_cptl_disc_curve);
                end
                
                clear('tmp_ptf_variation_margin')
                
            else
                
                % Capitalization factor for VM(t-1;t)
                coll_cptl_disc_curve          = ZCB_tT(r_disc,exposure_grid_yf(t-1),exposure_grid_yf(t),DF_disc(:,1)/365,interp_method,parameters,x_exp(:,t-1),y_exp(:,t-1));
                
                for i = 1:n_assets
                    % Results reordered
                    ptf_exposure_tMPOR{i,t}   = tmp_ptf_exposure_tMPOR(:,i);
                    % VM calculation
                    ptf_variation_margin{i,t} = collateral_simulation(csa_threshold,csa_mta,ptf_exposure_tMPOR{i,t},ptf_variation_margin{i,t-1},coll_cptl_disc_curve);
                end
                
            end
            
            clear('CF_disc_curve_coll_grid','CF_fwd_curve_coll_grid','tmp_ptf_exposure_tMPOR','coll_cptl_disc_curve')
            
        end
        
    else
        for i = 1:n_assets
            ptf_variation_margin{i,t} = zeros(n_sim,1);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculation of IM at time t for each MC scenario %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag_IM == 1
        
        if collateral_grid_yf(t) >= 0  % IM(0) = 0
            
            % Forwarding and discounting curves corresponding to ISDA-SIMM tenors
            SIMM_disc_curve_coll_grid                     = ZCB_tT(r_disc,collateral_grid_yf(t),collateral_grid_yf(t)+SIMM_pillars,DF_disc(:,1)/365,interp_method,parameters,x_coll(:,t),y_coll(:,t));
            SIMM_fwd_curve_coll_grid                      = ZCB_tT(r_fwd, collateral_grid_yf(t),collateral_grid_yf(t)+SIMM_pillars,DF_fwd(:,1)/365, interp_method,parameters,x_coll(:,t),y_coll(:,t));
            
            % Rates interpolated on cash flow dates
            residual_CF_dates_yf_coll                     = instruments_dates_yf(instruments_dates_yf >= collateral_grid_yf(t));
            [CF_SIMM_disc_curve_coll_grid, tmp_r_disc]    = DF_CF_dates(SIMM_disc_curve_coll_grid,collateral_grid_yf(t)+SIMM_pillars,residual_CF_dates_yf_coll',collateral_grid_yf(t),interp_method);
            [CF_SIMM_fwd_curve_coll_grid,  tmp_r_fwd]     = DF_CF_dates(SIMM_fwd_curve_coll_grid, collateral_grid_yf(t)+SIMM_pillars,residual_CF_dates_yf_coll',collateral_grid_yf(t),interp_method);
            
            % Get MtF at time t-tMPOR using curves interpolated on ISDA-SIMM tenors
            residual_CF_dates_coll                        = instruments_dates(instruments_dates >= col_evaluation_dates(t));
            [tmp_SIMM_und_mtf,...
                SIMM_swap_rate, SIMM_swap_annuity]        = PortfolioExposure(id_instruments,ptf_data,schedules,CF_SIMM_fwd_curve_coll_grid,CF_SIMM_disc_curve_coll_grid,col_evaluation_dates(t),simulated_libor_fixings,residual_CF_dates_coll,parameters,swpt_exercise);
            
            % Results reordered
            for i = 1:n_assets
                SIMM_mtf{i,t} = tmp_SIMM_und_mtf(:,i);        
            end
            
            %%%%%%%%%%%%%%%%%%%%
            %%% Delta Margin %%%
            %%%%%%%%%%%%%%%%%%%%

            % Construction of bumped curves
            for j = 1:size(SIMM_pillars,2) % loop on curve tenors from position 2 since position 1 is P(t,t) = 1
                r_disc_shock                              = [tmp_r_disc(:,1:j-1) tmp_r_disc(:,j)+0.0001 tmp_r_disc(:,j+1:end)];
                
                CF_SIMM_disc_curve_coll_grid_shock{j,1}   = DF_CF_dates(CF_SIMM_disc_curve_coll_grid,collateral_grid_yf(t)+SIMM_pillars,residual_CF_dates_yf_coll',collateral_grid_yf(t),interp_method,r_disc_shock);
                
                r_fwd_shock                               = [tmp_r_fwd(:,1:j-1) tmp_r_fwd(:,j)+0.0001 tmp_r_fwd(:,j+1:end)];
                CF_SIMM_fwd_curve_coll_grid_shock{j,1}    = DF_CF_dates(CF_SIMM_fwd_curve_coll_grid, collateral_grid_yf(t)+SIMM_pillars,residual_CF_dates_yf_coll',collateral_grid_yf(t),interp_method,r_fwd_shock);
            end
            clear j
            
            
            delta_margin(:,t)                             = SIMM_delta_margin(CF_SIMM_disc_curve_coll_grid_shock,CF_SIMM_fwd_curve_coll_grid_shock,CF_SIMM_disc_curve_coll_grid,CF_SIMM_fwd_curve_coll_grid,SIMM_mtf(:,t),id_instruments,ptf_data,schedules,...
                                                                              col_evaluation_dates(t),pricing_date,simulated_libor_fixings,residual_CF_dates_coll,parameters,swpt_exercise,FX_spot);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Vega and Curvature Margin %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Construction of bumped curves
            shock_parameters = parameters;
            SIMM_vega_disc_curve_coll_grid                = ZCB_tT(r_disc,collateral_grid_yf(t),collateral_grid_yf(t)+SIMM_pillars,DF_disc(:,1)/365,interp_method,shock_parameters,x_coll(:,t),y_coll(:,t));
            SIMM_vega_fwd_curve_coll_grid                 = ZCB_tT(r_fwd, collateral_grid_yf(t),collateral_grid_yf(t)+SIMM_pillars,DF_fwd(:,1)/365, interp_method,shock_parameters,x_coll(:,t),y_coll(:,t));
            
            % Rates interpolated on cash flow dates
            CF_SIMM_vega_disc_curve_coll_grid             = DF_CF_dates(SIMM_vega_disc_curve_coll_grid,collateral_grid_yf(t)+SIMM_pillars,residual_CF_dates_yf_coll',collateral_grid_yf(t),interp_method);
            CF_SIMM_vega_fwd_curve_coll_grid              = DF_CF_dates(SIMM_vega_fwd_curve_coll_grid, collateral_grid_yf(t)+SIMM_pillars,residual_CF_dates_yf_coll',collateral_grid_yf(t),interp_method);
            
            % Check on counterparty posting IM: -1 for the Counterperty, 1 for the Investor
            if flag_VM == 1
                for i = 1:n_assets
                    check_ctp_IM{i,t} = sign(ptf_uncollateralized_exposure{i,t} - ptf_variation_margin{i,t})*-1;
                end
            else
                for i = 1:n_assets
                    check_ctp_IM{i,t} = sign(ptf_uncollateralized_exposure{i,t})*-1;
                end
            end
          
            [vega_margin(:,t),curvature_margin(:,t),vega(:,t)] = SIMM_vega_HW2F(parameters,SIMM_mtf{:,t},SIMM_swap_rate,SIMM_swap_annuity,displacements,col_evaluation_dates(t),id_instruments,ptf_data,schedules,CF_SIMM_vega_fwd_curve_coll_grid,CF_SIMM_vega_disc_curve_coll_grid,...
                                                                                simulated_libor_fixings,residual_CF_dates_coll,swpt_exercise,FX_spot,check_ctp_IM(:,t));
            
            % Get discounted IM
            initial_margin(:,t) = DF_disc_curve_mkt_coll_grid(t)*(delta_margin(:,t) + vega_margin(:,t) + curvature_margin(:,t));
        end
        clear('SIMM_disc_curve_coll_grid','SIMM_fwd_curve_coll_grid','CF_SIMM_disc_curve_coll_grid','CF_SIMM_fwd_curve_coll_grid','tmp_r_disc','tmp_r_fwd','tmp_SIMM_und_mtf','SIMM_swap_rate','SIMM_swap_annuity','r_disc_shock',...
            'r_fwd_shock','CF_SIMM_disc_curve_coll_grid_shock','CF_SIMM_fwd_curve_coll_grid_shock','SIMM_vega_disc_curve_coll_grid','SIMM_vega_fwd_curve_coll_grid','residual_CF_dates_yf_coll','residual_CF_dates_coll')
    else
        initial_margin(:,t) = zeros(n_sim,1);
        
    end
    
end
elapsed_time = toc;

% Avoid NaN values for IM
if flag_IM == 1
    vega_margin(isnan(vega_margin))           = 0;
    curvature_margin(isnan(curvature_margin)) = 0;
    initial_margin(isnan(initial_margin))     = 0;
end

%% xVA calculation
% Get uncollateralized, collateralized with VM and collateralized with VM and IM exposures
exposure                 = find_exposure(ptf_uncollateralized_exposure,ptf_variation_margin,initial_margin);

% Get Expected Positive Exposure (EPE), Expected Negative Exposure (ENE),
% Expected Exposure (EE) and Expected Initial Margin (EIM) along with 3 sigma upper and lower bounds
n_sigma = 3;
[EPE,EPE_lb,EPE_ub,ENE,ENE_lb,ENE_ub,EE]   = find_metrics(exposure,n_sigma);
if flag_IM == 1
    [EIM,EIM_ub,EIM_lb]                        = find_metrics_IM(initial_margin,n_sigma);
end

% Get default probabilities from CDS curves of the Bank and of the Investor
filename_CDS             = 'CDS_spreads_28_12_18.xlsx';
CDSCurves                = xlsread(filename_CDS);           % Load CDS spreads from Excel file
CDSCurves_I              = CDSCurves(:,1:2);                % CDS curve of the Investor
CDSCurves_C              = [CDSCurves(:,1),CDSCurves(:,3)]; % CDS curve of the Counterperty
tau                      = (DF_disc_input(:,1)-pricing_date)/365;
ZeroData                 = [DF_disc_input(:,1) [1;-log(DF_disc_input(2:end,2))./tau(2:end)]];
MarketData_I             = [CDSCurves(:,1)+pricing_date CDSCurves(:,2)];
MarketData_C             = [CDSCurves(:,1)+pricing_date CDSCurves(:,3)];
RR_inv                   = 0.4;                             % define recovery rate for the Investor
RR_ctp                   = 0.4;                             % define recovery rate for the Counterperty
% Default probabilities
try 
    DP_inv_grid              = cdsbootstrap(ZeroData,MarketData_I,pricing_date,'RecoveryRate',RR_inv,'Period',4,'Basis',0,'BusDayConvention','actual','ProbDates',exp_evaluation_dates');
    DP_ctp_grid              = cdsbootstrap(ZeroData,MarketData_C,pricing_date,'RecoveryRate',RR_ctp,'Period',4,'Basis',0,'BusDayConvention','actual','ProbDates',exp_evaluation_dates');
catch ('Finance Toolbox not available, Default probabilities loaded from input file')
    load('DP_inv_grid_daily30y.mat');
    load('DP_ctp_grid_daily30y.mat')
    for i=1:length(exp_evaluation_dates)
        index(i)             = find(exp_evaluation_dates(i) == DP_inv_grid);
    end
    DP_inv_grid = DP_inv_grid(index,:);
    DP_ctp_grid = DP_ctp_grid(index,:);
end

DP_inv_grid              = DP_inv_grid(:,2);
SP_inv_grid              = 1-DP_inv_grid;

DP_ctp_grid              = DP_ctp_grid(:,2);
SP_ctp_grid              = 1-DP_ctp_grid;

% Construction of Funding Curve as CDS curve of the Investor + 100 bps 
MarketData_I_fund      = MarketData_I;
MarketData_I_fund(:,1) = MarketData_I_fund(:,1) - pricing_date;
MarketData_I_fund(:,2) = (MarketData_I_fund(:,2)+100)/10^4;

r_fund_exp_grid   = interp1(MarketData_I_fund(:,1)/365,MarketData_I_fund(:,2),exposure_grid_yf,interp_method,'extrap');
r_fund_exp_grid   = [0 r_fund_exp_grid(2:end)];
r_disc_exp_grid   = [0 -log(DF_disc_curve_mkt_exp_grid(2:end))./exposure_grid_yf(2:end)];

FundSpread_I_exp_grid = (r_fund_exp_grid(2:end)-r_disc_exp_grid(2:end)).*diff(exposure_grid_yf);
FundSpread_I_exp_grid = [FundSpread_I_exp_grid(1) FundSpread_I_exp_grid];

%%%%%%%%%%%%%%%%%%%
%%% CVA and DVA %%%
%%%%%%%%%%%%%%%%%%%
% Get CVA, DVA and related 3 sigma lower and upper bounds
[CVA,DVA,CVA_grid,DVA_grid] = find_cva_dva(EPE,ENE,RR_ctp,DP_ctp_grid,SP_ctp_grid,RR_inv,DP_inv_grid,SP_inv_grid);
[CVA_ub,DVA_ub]             = find_cva_dva(EPE_lb,ENE_lb,RR_ctp,DP_ctp_grid,SP_ctp_grid,RR_inv,DP_inv_grid,SP_inv_grid);
[CVA_lb,DVA_lb]             = find_cva_dva(EPE_ub,ENE_ub,RR_ctp,DP_ctp_grid,SP_ctp_grid,RR_inv,DP_inv_grid,SP_inv_grid);

%%%%%%%%%%%%%%%%%%%
%%% FVA and MVA %%%
%%%%%%%%%%%%%%%%%%%
% Get FCA, FBA and related 3 sigma lower and upper bounds
[FCA,FCA_grid]       = find_fmva(EPE   ,FundSpread_I_exp_grid,SP_ctp_grid,SP_inv_grid);
[FCA_lb,FCA_lb_grid] = find_fmva(EPE_ub,FundSpread_I_exp_grid,SP_ctp_grid,SP_inv_grid);
[FCA_ub,FCA_ub_grid] = find_fmva(EPE_lb,FundSpread_I_exp_grid,SP_ctp_grid,SP_inv_grid);

[FBA,FBA_grid]       = find_fmva(ENE   ,FundSpread_I_exp_grid,SP_ctp_grid,SP_inv_grid);
[FBA_lb,FBA_lb_grid] = find_fmva(ENE_ub,FundSpread_I_exp_grid,SP_ctp_grid,SP_inv_grid);
[FBA_ub,FBA_ub_grid] = find_fmva(ENE_lb,FundSpread_I_exp_grid,SP_ctp_grid,SP_inv_grid);

% Get MVA and related 3 sigma lower and upper bounds
if flag_IM == 1
    [MVA,MVA_grid]       = find_fmva(EIM   ,FundSpread_I_exp_grid.*SP_ctp_grid'.*SP_inv_grid');
    [MVA_lb,MVA_lb_grid] = find_fmva(EIM_ub,FundSpread_I_exp_grid.*SP_ctp_grid'.*SP_inv_grid');
    [MVA_ub,MVA_ub_grid] = find_fmva(EIM_lb,FundSpread_I_exp_grid.*SP_ctp_grid'.*SP_inv_grid');
end