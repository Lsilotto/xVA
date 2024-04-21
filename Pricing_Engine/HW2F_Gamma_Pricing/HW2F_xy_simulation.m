function [r_out_d,r_out_f,x,y,Phi_d,Phi_f] = HW2F_xy_simulation(Parameters,DF_disc,DF_fwd,T,n,M,time_grid,flag_seed) 
%{
% The function simulates the short rates for G2++ model with time dependent volatility
%
%%% input:
%
% Parameters = matrix n_expiries+1 by 5 cosi costruita:
%              - row1 calibrated parameters
%                col1 from 2:end expiries of market swaptions used for calibration
%              - col2 from 2:end Gamma calibrated values 
%
% DF_disc    = martrix nx2: 
%               - col1 n time to maturities in days, 
%               - col2 discount factors (discounting curve)
%
% DF_fwd     = martrix nx2: 
%               - col1 n time to maturities in days, 
%               - col2 discount factors (forwarding curve)
%
% T          = time to set T-forward measure (in year fraction)
%
% n          = number of time grid points
%
% M          = number of MC scenarios
%
% time_grid  = simulation time grid vector in year fraction
%
%%% output:
%
% r_out_d    = matrix M by n containing the paths of simulated short rates
%              (discounting)
%
% r_out_f    = matrix M by n containing the paths of simulated short rates
%              (forwarding)
%
% x          = matrix M by n containing the paths of process x
%
% y          = matrix M by n containing the paths of process y
%
% Phi_d      = matrix M by n containing the paths of phi (discounting)
%
% Phi_f      = matrix M by n containing the paths of phi (forwarding)
%}

%%% Compute zero rates from discount factors 
% Discounting curve
tmp1 = -(1./(DF_disc(:,1)./365)).*log(DF_disc(:,2));
B_rate_d = [DF_disc(:,1) tmp1];
% Forwarding curve
tmp2 = -(1./(DF_fwd(:,1)./365)).*log(DF_fwd(:,2));
B_rate_f = [DF_fwd(:,1) tmp2];

%%% Merge of term structures and interpolation
Curve_Pillars = unique([B_rate_d(:,1);B_rate_f(:,1)]);
B_rate_d_temp  = interp1([0; B_rate_d(:,1)],[B_rate_d(1,2); B_rate_d(:,2)],Curve_Pillars,'spline');
B_rate_f_temp  = interp1([0; B_rate_f(:,1)],[B_rate_f(1,2); B_rate_f(:,2)],Curve_Pillars,'spline');

%%% Compute instantaneous forward rates
% f_x(t,T) = z_x(t,T)+ d/dT{z_x(t,T)}tau_x(t,T)
% pag 15 Everything... Bianchetti
B_rate_d_int   = interp1(Curve_Pillars,B_rate_d_temp,time_grid,'spline');
B_rate_f_int   = interp1(Curve_Pillars,B_rate_f_temp,time_grid,'spline');

dt = T/n; 
fwd_rate_d = B_rate_d_int+time_grid.*(diff([B_rate_d_int(1),B_rate_d_int])./dt);
fwd_rate_f = B_rate_f_int+time_grid.*(diff([B_rate_f_int(1),B_rate_f_int])./dt);

%%% G2++ parameters initialization
Params = Parameters(1,:);
a = Parameters(1,1);
Sigma = Parameters(1,2);
b = Parameters(1,3);
eta = Parameters(1,4);
rho = Parameters(1,5);
T_Gamma = Parameters(2:end,1);
Gamma   = Parameters(2:end,2);
Gamma_intrpl = interp1(T_Gamma,Gamma,time_grid,'next');
Gamma_intrpl(isnan(Gamma_intrpl)) = Gamma(1);

%%% Compute Phi(t)  
[r0_d,Phi_d] = Phi_T_Gamma(Params,fwd_rate_d,Gamma_intrpl,time_grid);
[r0_f,Phi_f] = Phi_T_Gamma(Params,fwd_rate_f,Gamma_intrpl,time_grid);

%%% Simulation initialization
r_out_d       = zeros(M,length(time_grid));
r_out_d(:,1)  = r0_d;
r_out_f       = zeros(M,length(time_grid));
r_out_f(:,1)  = r0_f;
x             = zeros(M,length(time_grid));
y             = zeros(M,length(time_grid));

%%% Simulation
for i=1:length(time_grid)-1
    var_x    = VarX_HW2F_Gamma([[a Sigma];[T_Gamma Gamma]],time_grid(i),time_grid(i+1)); 
    var_y    = VarY_HW2F_Gamma([[b eta];[T_Gamma Gamma]],time_grid(i),time_grid(i+1));
    cov_xy   = rho*VarX_HW2F_Gamma([[(a+b)/2 sqrt(Sigma*eta)];[T_Gamma Gamma]],time_grid(i),time_grid(i+1));
    
    cov_matrix    = [var_x,cov_xy;cov_xy,var_y];    
% Compute Mx and My
[Mx,My]  = MT_xy_HW2F_Gamma(Parameters,T,time_grid(i),time_grid(i+1));

% Random numbers
if flag_seed == 1
    rng(i) % fixed seed
end
Temp = mvnrnd([0,0],cov_matrix,M/2);
V = [Temp;-Temp]; 


% Compute paths x and y
x(:,i+1) =  x(:,i)*exp(-a*(time_grid(i+1)-time_grid(i))) + Mx + V(:,1); 
y(:,i+1) =  y(:,i)*exp(-b*(time_grid(i+1)-time_grid(i))) + My + V(:,2);

% measure Q
% x(:,i+1) =  x(:,i)*exp(-a*(time_grid(i+1)-time_grid(i))) + V(:,1); 
% y(:,i+1) =  y(:,i)*exp(-b*(time_grid(i+1)-time_grid(i))) + V(:,2); 

% Compute short rates
r_out_d(:,i+1) = x(:,i+1) + y(:,i+1) + Phi_d(1,i+1);
r_out_f(:,i+1) = x(:,i+1) + y(:,i+1) + Phi_f(1,i+1);
end



