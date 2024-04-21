function[variation_margin] = collateral_simulation(threshold,minimum_transfer_amount,netting_set_exposure,collateral_account,DF_disc_simulated)
%
% The function simulates the evolution of the collateral account
%
% threshold                   = CSA threshold 
%
% minimum_transfer_amount     = CSA minimum transfer amount 
%
% netting_set_exposure        = matrix n_sim by n_time_steps containing exposure values 
%
% exposure_evaluation_dates   = exposure simulation dates (Matlab date format) 
%
%%% output:
%
% variation_margin
%

collateral_account = collateral_account./DF_disc_simulated;
    tmp1_ind1                       = netting_set_exposure + threshold < 0;      
    tmp2_ind1                       = collateral_account < 0; 
    % see formula 13.8 Brigo,Morini,Pallavicini
    indicatore1                     = abs(tmp1_ind1.*(netting_set_exposure + threshold) - tmp2_ind1.*collateral_account) > minimum_transfer_amount; 

    tmp1_ind2                       = netting_set_exposure - threshold > 0;      
    tmp2_ind2                       = collateral_account > 0; 
    % see formula 13.9 Brigo,Morini,Pallavicini
    indicatore2                     = abs(tmp1_ind2.*(netting_set_exposure - threshold) - tmp2_ind2.*collateral_account) > minimum_transfer_amount; 

    % see formula 13.10 Brigo,Morini,Pallavicini
    variation_margin                = collateral_account + indicatore1.*(min(netting_set_exposure + threshold,0) - min(collateral_account,0))...
                                    + indicatore2.*(max(netting_set_exposure - threshold,0) - max(collateral_account,0));
end
