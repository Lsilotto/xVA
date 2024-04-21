function[initial_margin] = IM_k_mta(threshold,minimum_transfer_amount,IM)


    

    indicatore1                     = IM - threshold > minimum_transfer_amount; 
 

    % formula 13.10 Brigo,Morini,Pallavicini
    initial_margin                  = indicatore1.*(max(IM - threshold,0));
                                    
end
