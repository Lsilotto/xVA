function[EIM,EIM_ub,EIM_lb] = find_metrics_IM(initial_margin,n_sigma)

n_sim            = size(initial_margin,1);
time_grid_points = size(initial_margin,2);

for t = 1:time_grid_points

    EIM(t)    = 1/n_sim*sum(initial_margin(:,t));
    sigma     = sqrt((1/(n_sim-1))*sum((initial_margin(:,t)-EIM(t)).^2));
    MC_error  = sigma/sqrt(n_sim);
    EIM_ub(t) = EIM(t) + n_sigma*MC_error;
    EIM_lb(t) = EIM(t) - n_sigma*MC_error;
        
end