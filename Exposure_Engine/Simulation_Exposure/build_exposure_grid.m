function[exposure_grid,exposure_grid_yf,x_grid,y_grid] = build_exposure_grid(pricing_date,T,grid_basis,grid_interval,x,y,xy_time_grid,instrument_CF_dates)
%{
% The function builds the time grid used for simulating the exposure.
% Starting from the pricing date (today) the time points are added according to the basis (day,month)
% and the granularity (interval) untill the instrument (or portfolio) maturity.
% It is possible to add additional time points using the last argument of
% the function (optional).
%
%%% input:
%
% pricing_date  = today (Matlab date format)
%
% T             = last simulation point (Matlab date format)
%
% grid_basis    = grid basis: 'month' -> monthly time grid ,'day' -> daily time grid  
%            
% grid_interval = distance between time points
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
% exposure_grid    = time grid used for simulating the exposure  (Matlab date format)
%
% exposure_grid_yf = time grid used for simulating the exposure  (in year fraction)
%
% x_grid           = G2++ x process in correspondence of exposure simulation dates
%
% y_grid           = G2++ y process in correspondence of exposure simulation dates
%}


% starting from the pricing date time points are untill maturity
exposure_grid(1)  = pricing_date;
i = 1;
while exposure_grid(i) <= addtodate(T,-grid_interval,grid_basis)
    exposure_grid(i+1)  = addtodate(exposure_grid(i),grid_interval,grid_basis);
    i = i+1;
end
% ensuring the last date is included in the time grid
exposure_grid     = unique([exposure_grid T]);

% Joint grid construction: CF_dates + 1 day are added 
if nargin == 8
    exposure_grid = sort(unique([exposure_grid (instrument_CF_dates+1)]),'ascend');
end

exposure_grid_yf  = (exposure_grid-pricing_date)/365;

% x and y processes in correspondence of exposure simulation dates
if verLessThan ('matlab','9.4') == 0
    for i=1:length(exposure_grid_yf)
        index_sim(i)                = find(round(exposure_grid_yf(i),4) == round(xy_time_grid,4));
    end
else
    for i=1:length(exposure_grid_yf)
        index_sim(i)                = find(round(exposure_grid_yf(i)*10000)/10000 == round(xy_time_grid*10000)/10000); % workaround per versioni matlab < 2018a 
    end
end

x_grid                          = x(:,index_sim);
y_grid                          = y(:,index_sim);


