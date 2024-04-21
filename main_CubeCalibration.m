% This code calibrates G2++ parameters on swaption prices cube 

%% Import market data from Excel file
[d1,d2,d3] = textread('TargetDates.txt','%s %s %s'); % holiday calendar
dv = strcat(d1,d2,d3);
dv = datenum(dv);
Date     = '28_12_2018';
filename = ['MarketData_Input_' Date '.xls'];        % market data from excel file

DF_input       = xlsread(filename,'Discounting');    % discounting curve
DF_input(:,1)  = x2mdate(DF_input(:,1));
FWD_input      = xlsread(filename,'Forwarding');     % forwarding curve
FWD_input(:,1) = x2mdate(FWD_input(:,1));

Pricing_Date  = DF_input(1,1);           % set pricing date
[yy,mm,dd]    = datevec(Pricing_Date);
DF       = DF_input;                     % discounting curve Discount Factors
DF(:,1)  = DF(:,1) - Pricing_Date;
FWD      = FWD_input;                    % forwarding curve Discount Factors
FWD(:,1) = FWD(:,1) - Pricing_Date;

filename = 'Input_28_12_EOM_CMS_calibrator.xlsx'; % swaption cube from excel file
% ATM 
ATM_Volatilities  = xlsread(filename,'Vola_Physical_ATM_input'); 
ATM_Expiry        = ATM_Volatilities(:,1);
ATM_Tenor         = ATM_Volatilities(:,2);
ATM_Volatilities  = ATM_Volatilities(:,3);
% Cube 
Prices_Physical = xlsread(filename,'Cube_Physical_Prices_input');
Cube_Expiry     = Prices_Physical(:,1);
Cube_Tenor      = Prices_Physical(:,2);
Cube_Omega      = Prices_Physical(:,3);
Cube_shiftK     = Prices_Physical(:,4)/10000;
Cube_Price_Physical = Prices_Physical(:,5);
ATM_Price_Physical  = Prices_Physical(6:11:end,5);
ATM_Omega           = Prices_Physical(6:11:end,3);
% Volatility log-normal shifts
ATM_Displacements = xlsread(filename,'Displacements_input');
ATM_Displacements = ATM_Displacements(:,3)/100;

%% Exclusion of extreme points of volatility cube (illiquid points)
% Expiry < 2Y
ATM_Displacements  = ATM_Displacements(26:end);
ATM_Expiry         = ATM_Expiry(26:end);
ATM_Omega          = ATM_Omega(26:end);
ATM_Price_Physical = ATM_Price_Physical(26:end);
ATM_Tenor          = ATM_Tenor(26:end);
ATM_Volatilities   = ATM_Volatilities(26:end);

ATM_nexp = size(ATM_Expiry,1);
ATM_nten = size(ATM_Tenor,1);

Cube_Expiry         = Cube_Expiry(276:end);
Cube_Omega          = Cube_Omega(276:end);
Cube_Price_Physical = Cube_Price_Physical(276:end);
Cube_Tenor          = Cube_Tenor(276:end);
Cube_shiftK         = Cube_shiftK(276:end);

% Expiry > 15Y
ATM_Displacements  = ATM_Displacements(1:20);
ATM_Expiry         = ATM_Expiry(1:20);
ATM_Omega          = ATM_Omega(1:20);
ATM_Price_Physical = ATM_Price_Physical(1:20);
ATM_Tenor          = ATM_Tenor(1:20);
ATM_Volatilities   = ATM_Volatilities(1:20);

ATM_nexp = size(ATM_Expiry,1);
ATM_nten = size(ATM_Tenor,1);

Cube_Expiry         = Cube_Expiry(1:220);
Cube_Omega          = Cube_Omega(1:220);
Cube_Price_Physical = Cube_Price_Physical(1:220);
Cube_Tenor          = Cube_Tenor(1:220);
Cube_shiftK         = Cube_shiftK(1:220);

% Zero rates
tau=(DF(2:end,1))/365;
Rf=-log(DF(2:end,2))./tau;
Rf=[0;Rf];

dcbasis_Fix=6;   % 30/360 European
dcbasis_Float=2; % actual/360

%% Calibration
% Get Swap rate from price cube
RePrices_Physical = ATM_Price_Physical*0;
ATM_Swap_Rate     = RePrices_Physical;
ATM_Swap_T        = zeros(ATM_nexp,1);

for i=1:ATM_nexp
    
    Expiry_Swaption = busdate(datenum([yy mm+ATM_Expiry(i) dd])-1,1,dv);
    Start_DateX     = busdate(busdate(Expiry_Swaption,1,dv),1,dv);
    
    Coupon_Index_S  = Swap_Date_Generator(Start_DateX, ATM_Tenor(i),12,6,dv);
    Notional_S      = Coupon_Index_S;
    [Swap_Rate, Swap_Rate1, FWDX, Annuity, Annuity_Cash] = IRS(DF,FWD,Pricing_Date,Notional_S,Coupon_Index_S,dcbasis_Fix,dcbasis_Float);
    Tmat = (Expiry_Swaption-Pricing_Date)/365;
    
    [call(i,1)]          = blsprice(Swap_Rate+ATM_Displacements(i),Swap_Rate+ATM_Displacements(i),0,Tmat,ATM_Volatilities(i));
    discount_expiry      = exp(-interp1(Pricing_Date+DF(:,1),Rf,Expiry_Swaption,'linear').*(Expiry_Swaption-Pricing_Date)/365);
    RePrices_Physical(i) = call(i,1)*Annuity/discount_expiry;
    
    ATM_Swap_Rate(i)= Swap_Rate;
    ATM_Swap_T(i)   = Tmat;
    Ann(i,1)        = Annuity;
    
end

tmp_Cube_Strike = [];
for i = 1:length(ATM_Swap_Rate)
    tmp_Cube_Strike = [tmp_Cube_Strike; ATM_Swap_Rate(i).*ones(size(unique(Cube_shiftK),1),1)];
end
Cube_Strike = tmp_Cube_Strike + Cube_shiftK;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% G2++ calibration constant volatility %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% guessing points
x = [0.8248 ;   0.0273   ; 0.0298 ;   0.0086 ;  -0.9762]';
% calibration upper & lower bound
LB = [0.1; 0.001; 0.0001; 0.001; -1]';
UB = [1.5; 1; 0.25; 0.25; -0.5]';

[Pars_HW2F_new, Swapt_Cube2F, resnorm_HW2F_new] = Calib_HW2F_cube([x' LB' UB'], DF, FWD,Pricing_Date, 0, 0, Cube_Expiry,Cube_Tenor,Cube_Strike, 0, 0,0,Cube_Price_Physical/10^4,0,0,Cube_Omega);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% G2++ calibration time dep. volatility %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calibration upper & lower bound
LB_HW2F = [0.1; 0.001; 0.0001; 0.001; -1]';
UB_HW2F = [1.5; 1; 0.25; 0.25; -0.5]';
% guessing points
Pars_HW2F_Gamma = [[Pars_HW2F_new zeros(5,1)];[unique(ATM_Swap_T) ones(length(unique(Cube_Expiry)),1)]];

LB_Gamma = [LB_HW2F' ; zeros(length(unique(Cube_Expiry)),1)];
UB_Gamma = [UB_HW2F' ; 10*ones(length(unique(Cube_Expiry)),1)];

[Pars_HW2F_Gamma_new, Swapt_Cube2F_Gamma, resnorm_HW2F_Gamma_new] = Calib_HW2F_Gamma_cube([Pars_HW2F_Gamma LB_Gamma UB_Gamma] ,DF,FWD,...
                                                                                           Pricing_Date,0,0,Cube_Expiry,Cube_Tenor,Cube_Strike,0,0,0,Cube_Price_Physical/10^4,0,0,Cube_Omega);