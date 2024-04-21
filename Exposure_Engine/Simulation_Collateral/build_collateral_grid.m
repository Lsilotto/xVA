function[collateral_grid,collateral_grid_yf,x_grid,y_grid] = build_collateral_grid(exposure_grid,MPOR_dd,x,y,xy_time_grid)
%{
% The function builds the time grid used for simulating the collateral
% starting from the time grid used for simulating the exposure.
%
%%% input:
%
% exposure_grid = time grid used for simulating the exposure (Matlab date format)
%
% MPOR_dd       = Margin Period of Risk in days
%
% x             = matrix n_path by length(xy_time_grid) -> simulated G2++ x
%                 process
%
% y             = matrix n_path by length(xy_time_grid) -> simulated G2++ y
%                 process
%
% xy_time_grid  = x and y processes daily simulation time grid
%
% instruments_CF_dates (optional) = instruments cash flow dates 
%
%%% output:
%
% collateral_grid    = time grid used for simulating the collateral  (Matlab date format)
%
% collateral_grid_yf = time grid used for simulating the collateral  (in year fraction)
%
% x_grid             = G2++ x process in correspondence of collateral simulation dates
%
% y_grid             = G2++ y process in correspondence of collateral simulation dates
%}


collateral_grid    = exposure_grid - MPOR_dd;

pricing_date       = exposure_grid(1);
collateral_grid_yf = (collateral_grid - pricing_date)/365;

% x and y processes in correspondence of collateral simulation dates
if verLessThan ('matlab','9.4') == 0
    for i = 1:length(collateral_grid_yf)
        if collateral_grid_yf(i) >= 0
            index_sim(i)    = find(round(collateral_grid_yf(i),4) == round(xy_time_grid,4));
        else
            index_sim(i)    = 1;
        end
    end
else
    for i = 1:length(collateral_grid_yf)
        if collateral_grid_yf(i) >= 0
            index_sim(i)    = find(round(collateral_grid_yf(i)*10000)/10000 == round(xy_time_grid*10000)/10000); % workaround per versioni matlab < 2018a 
        else
            index_sim(i)    = 1;
        end
    end
end

x_grid                  = x(:,index_sim);
y_grid                  = y(:,index_sim);


