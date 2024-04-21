function[CVA,DVA,CVA_grid,DVA_grid] = find_cva_dva(EPE,ENE,RR_ctp,DP_ctp_grid,SP_ctp_grid,RR_inv,DP_inv_grid,SP_inv_grid)

tmp_EPE                   = struct2cell(EPE);
tmp_ENE                   = struct2cell(ENE);

% Uncollateralized
CVA_grid.uncollateralized = -(1-RR_ctp)*[0; diff(DP_ctp_grid)].*SP_inv_grid.*tmp_EPE{1}'; 
CVA.uncollateralized      = sum(CVA_grid.uncollateralized);
DVA_grid.uncollateralized = -(1-RR_inv)*[0; diff(DP_inv_grid)].*SP_ctp_grid.*tmp_ENE{1}'; 
DVA.uncollateralized      = sum(DVA_grid.uncollateralized);

% w VM
CVA_grid.wVM              = -(1-RR_ctp)*[0; diff(DP_ctp_grid)].*SP_inv_grid.*tmp_EPE{2}'; 
CVA.wVM                   = sum(CVA_grid.wVM);
DVA_grid.wVM              = -(1-RR_inv)*[0; diff(DP_inv_grid)].*SP_ctp_grid.*tmp_ENE{2}'; 
DVA.wVM                   = sum(DVA_grid.wVM);

% w IM
CVA_grid.wIM              = -(1-RR_ctp)*[0; diff(DP_ctp_grid)].*SP_inv_grid.*tmp_EPE{3}'; 
CVA.wIM                   = sum(CVA_grid.wIM);
DVA_grid.wIM              = -(1-RR_inv)*[0; diff(DP_inv_grid)].*SP_ctp_grid.*tmp_ENE{3}'; 
DVA.wIM                   = sum(DVA_grid.wIM);

% w VM + IM
CVA_grid.wVM_IM           = -(1-RR_ctp)*[0; diff(DP_ctp_grid)].*SP_inv_grid.*tmp_EPE{4}'; 
CVA.wVM_IM                = sum(CVA_grid.wVM_IM);
DVA_grid.wVM_IM           = -(1-RR_inv)*[0; diff(DP_inv_grid)].*SP_ctp_grid.*tmp_ENE{4}'; 
DVA.wVM_IM                = sum(DVA_grid.wVM_IM);