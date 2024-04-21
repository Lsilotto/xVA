function[exposure] = find_exposure(ptf_uncollateralized_exposure,ptf_variation_margin,initial_margin)

n_time_steps = length(ptf_uncollateralized_exposure);
n_sim        = length(ptf_uncollateralized_exposure{1});

ptf_variation_margin{end} = ptf_variation_margin{end}.*0;
initial_margin(:,end)     = initial_margin(:,end).*0;

exposure.uncollateralized = cell2mat(ptf_uncollateralized_exposure);
variation_margin          = cell2mat(ptf_variation_margin);


for t = 1:n_time_steps
    
    for n = 1:n_sim
        % Exposure with only VM
        exposure.wVM(n,t) = exposure.uncollateralized(n,t) - variation_margin(n,t);
        
        % Exposure with only VM IM:
        if exposure.uncollateralized(n,t) > 0
            exposure.wIM(n,t) = max(exposure.uncollateralized(n,t) - initial_margin(n,t),0);
        elseif exposure.uncollateralized(n,t) < 0
            exposure.wIM(n,t) = min(exposure.uncollateralized(n,t) + initial_margin(n,t),0);
        else
            exposure.wIM(n,t) = 0;
        end
        
        % Exposure with VM and IM
        if exposure.wVM(n,t) > 0
            exposure.wVM_IM(n,t) = max(exposure.wVM(n,t) - initial_margin(n,t),0); % initial_margin > 0 (includes convexity component)
        elseif exposure.wVM(n,t) < 0
            exposure.wVM_IM(n,t) = min(exposure.wVM(n,t) + initial_margin(n,t),0); % initial_margin > 0 (includes convexity component)
        else
            exposure.wVM_IM(n,t) = 0;
        end
        
    end
end








