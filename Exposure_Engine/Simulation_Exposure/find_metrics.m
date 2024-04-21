function[EPE,EPE_lb,EPE_ub,ENE,ENE_lb,ENE_ub,EE,PFE_90,PFE_10] = find_metrics(exposure,n_sigma)

tmp              = struct2cell(exposure);
time_grid_points = size(tmp{1},2);
n_sim            = size(tmp{1},1);

for t = 1:time_grid_points
    
    % Uncollateralized
    EPE.uncollateralized(t)            = 1/n_sim*sum(max(exposure.uncollateralized(:,t),0));
    sigma_EPE                          = sqrt((1/(n_sim-1))*sum((max(exposure.uncollateralized(:,t),0)-EPE.uncollateralized(t)).^2));
    MC_error_EPE                       = sigma_EPE/sqrt(n_sim);
    EPE_ub.uncollateralized(t)         = EPE.uncollateralized(t) + n_sigma*MC_error_EPE;
    EPE_lb.uncollateralized(t)         = EPE.uncollateralized(t) - n_sigma*MC_error_EPE;
    
    ENE.uncollateralized(t)            = 1/n_sim*sum(min(exposure.uncollateralized(:,t),0));
    sigma_ENE                          = sqrt((1/(n_sim-1))*sum((min(exposure.uncollateralized(:,t),0)-ENE.uncollateralized(t)).^2));
    MC_error_ENE                       = sigma_ENE/sqrt(n_sim);
    ENE_ub.uncollateralized(t)         = ENE.uncollateralized(t) + n_sigma*MC_error_ENE;
    ENE_lb.uncollateralized(t)         = ENE.uncollateralized(t) - n_sigma*MC_error_ENE;
    
    EE.uncollateralized(t)             = 1/n_sim*sum(exposure.uncollateralized(:,t));
    
    PFE_90.uncollateralized(t)         = prctile(exposure.uncollateralized(:,t),90);
    PFE_10.uncollateralized(t)         = prctile(exposure.uncollateralized(:,t),10); 
    
    % wVM
    EPE.wVM(t)            = 1/n_sim*sum(max(exposure.wVM(:,t),0));
    sigma_EPE             = sqrt((1/(n_sim-1))*sum((max(exposure.wVM(:,t),0)-EPE.wVM(t)).^2));
    MC_error_EPE          = sigma_EPE/sqrt(n_sim);
    EPE_ub.wVM(t)         = EPE.wVM(t) + n_sigma*MC_error_EPE;
    EPE_lb.wVM(t)         = EPE.wVM(t) - n_sigma*MC_error_EPE;
    
    ENE.wVM(t)            = 1/n_sim*sum(min(exposure.wVM(:,t),0));
    sigma_ENE             = sqrt((1/(n_sim-1))*sum((min(exposure.wVM(:,t),0)-ENE.wVM(t)).^2));
    MC_error_ENE          = sigma_ENE/sqrt(n_sim);
    ENE_ub.wVM(t)         = ENE.wVM(t) + n_sigma*MC_error_ENE;
    ENE_lb.wVM(t)         = ENE.wVM(t) - n_sigma*MC_error_ENE;
    
    EE.wVM(t)             = 1/n_sim*sum(exposure.wVM(:,t));
    
    PFE_90.wVM(t)         = prctile(exposure.wVM(:,t),90);
    PFE_10.wVM(t)         = prctile(exposure.wVM(:,t),10); 
    
    % wIM
    EPE.wIM(t)            = 1/n_sim*sum(max(exposure.wIM(:,t),0));
    sigma_EPE             = sqrt((1/(n_sim-1))*sum((max(exposure.wIM(:,t),0)-EPE.wIM(t)).^2));
    MC_error_EPE          = sigma_EPE/sqrt(n_sim);
    EPE_ub.wIM(t)         = EPE.wIM(t) + n_sigma*MC_error_EPE;
    EPE_lb.wIM(t)         = EPE.wIM(t) - n_sigma*MC_error_EPE;
    
    ENE.wIM(t)            = 1/n_sim*sum(min(exposure.wIM(:,t),0));
    sigma_ENE             = sqrt((1/(n_sim-1))*sum((min(exposure.wIM(:,t),0)-ENE.wIM(t)).^2));
    MC_error_ENE          = sigma_ENE/sqrt(n_sim);
    ENE_ub.wIM(t)         = ENE.wIM(t) + n_sigma*MC_error_ENE;
    ENE_lb.wIM(t)         = ENE.wIM(t) - n_sigma*MC_error_ENE;
    
    EE.wIM(t)             = 1/n_sim*sum(exposure.wIM(:,t));
    
    PFE_90.wIM(t)         = prctile(exposure.wIM(:,t),90);
    PFE_10.wIM(t)         = prctile(exposure.wIM(:,t),10); 
    
    % wVM + IM
    EPE.wVM_IM(t)         = 1/n_sim*sum(max(exposure.wVM_IM(:,t),0));
    sigma_EPE             = sqrt((1/(n_sim-1))*sum((max(exposure.wVM_IM(:,t),0)-EPE.wVM_IM(t)).^2));
    MC_error_EPE          = sigma_EPE/sqrt(n_sim);
    EPE_ub.wVM_IM(t)      = EPE.wVM_IM(t) + n_sigma*MC_error_EPE;
    EPE_lb.wVM_IM(t)      = EPE.wVM_IM(t) - n_sigma*MC_error_EPE;
    
    ENE.wVM_IM(t)         = 1/n_sim*sum(min(exposure.wVM_IM(:,t),0));
    sigma_ENE             = sqrt((1/(n_sim-1))*sum((min(exposure.wVM_IM(:,t),0)-ENE.wVM_IM(t)).^2));
    MC_error_ENE          = sigma_ENE/sqrt(n_sim);
    ENE_ub.wVM_IM(t)      = ENE.wVM_IM(t) + n_sigma*MC_error_ENE;
    ENE_lb.wVM_IM(t)      = ENE.wVM_IM(t) - n_sigma*MC_error_ENE;
    
    EE.wVM_IM(t)          = 1/n_sim*sum(exposure.wVM_IM(:,t));
    
    PFE_90.wVM_IM(t)      = prctile(exposure.wVM_IM(:,t),90);
    PFE_10.wVM_IM(t)      = prctile(exposure.wVM_IM(:,t),10); 
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
