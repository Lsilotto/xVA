% This code calibrates G2++ parameters on ATM swaption prices 

%% Import market data from Excel file
[d1,d2,d3] = textread('TargetDates.txt','%s %s %s'); % holiday calendar
dv = strcat(d1,d2,d3);
dv = datenum(dv);
Date = '28_12_2018';
filename = ['MarketData_Input_' Date '.xls'];        % market data from excel file

DF_input      = xlsread(filename,'Discounting');     % discounting curve
DF_input(:,1) = x2mdate(DF_input(:,1));
FWD_input      = xlsread(filename,'Forwarding');     % forwarding curve
FWD_input(:,1) = x2mdate(FWD_input(:,1));

Pricing_Date  = DF_input(1,1);           % set pricing date
[yy,mm,dd]    = datevec(Pricing_Date);
DF       = DF_input;                     % discounting curve Discount Factors
DF(:,1)  = DF(:,1) - Pricing_Date;
FWD      = FWD_input;                    % forwarding curve Discount Factors
FWD(:,1) = FWD(:,1) - Pricing_Date;

Volatilities  = xlsread(filename,'Vola_Physical'); % ATM swaption volatilities
Expiry        = Volatilities(2:end,1);
Tenor         = Volatilities(1,2:end)';
nexp = size(Expiry,1);
nten = size(Tenor,1);

Prices_Physical = xlsread(filename,'Price_Physical'); % ATM swaption prices
Prices_Physical = Prices_Physical(2:end,2:end);

Displacements = xlsread(filename,'Displacements'); % Volatility log-normal shifts
Displacements = Displacements(2:end,2:end)/100;

%% Exclusion of extreme points of volatility cube (illiquid points)
% Tenor 1Y
Tenor = Tenor(2:end);
nten = size(Tenor,1);
Volatilities = Volatilities(:,2:end);
Prices_Physical = Prices_Physical(:,2:end);
Displacements = Displacements(:,2:end);

% Expiry < 2Y
Expiry = Expiry(8:end);
nexp = size(Expiry,1);
Volatilities = Volatilities(8:end,:);
Prices_Physical = Prices_Physical(8:end,:);
Displacements = Displacements(8:end,:);

% Expiry > 15Y
Expiry = Expiry(1:11);
nexp = size(Expiry,1);
Volatilities = Volatilities(1:11,:);
Prices_Physical = Prices_Physical(1:11,:);
Displacements = Displacements(1:11,:);

% Zero rates
tau=(DF(2:end,1))/365;
Rf=-log(DF(2:end,2))./tau;
Rf=[0;Rf];

dcbasis_Fix=6;   % 30/360 European
dcbasis_Float=2; % actual/360

%% Calibration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% G2++ calibration constant volatility %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% guessing points
x = [0.8248 ;   0.0273   ; 0.0298 ;   0.0086 ;  -0.9762]';
% calibration upper & lower bound
LB = [0.1; 0.001; 0.0001; 0.001; -1]';
UB = [1.5; 1; 0.25; 0.25; -0.5]';

[Pars_HW2F_new, Swapt_ATM2F, residual_HW2F_new] = Calib_HW2F([x' LB' UB'], DF, FWD,Pricing_Date, Expiry, Tenor, 0,0,0, Swap_Rate_ATM, 0,Prices_Physical/10^4, 0, 0,0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% G2++ calibration time dep. volatility %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calibration upper & lower bound
LB_HW2F = [0.1; 0.001; 0.0001; 0.001; -1]';
UB_HW2F = [1.5; 1; 0.25; 0.25; -0.5]';
% guessing point
Pars_HW2F_Gamma = [[Pars_HW2F_new zeros(5,1)];[T_ATM_Swap(:,1) ones(nexp,1)]];

LB_Gamma = [LB_HW2F' ; zeros(nexp,1)];
UB_Gamma = [UB_HW2F' ; 10*ones(nexp,1)];

[Pars_HW2F_Gamma_new, Swapt_ATM2F_Gamma, Swapt_Cube_2F_Gamma, resnormHW2F_Gamma_new] = Calib_HW2F_Gamma([Pars_HW2F_Gamma LB_Gamma UB_Gamma] ,DF,FWD,...
                                                                                                         Pricing_Date,Expiry,Tenor,0,0,0,Swap_Rate_ATM,0,Prices_Physical/10^4,0,0,0);