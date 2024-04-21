function[XVA,XVA_grid] = find_fmva(EXP,FUND,SP_1,SP_2)

if nargin == 2 % MVA

    % Uncollateralized
    XVA_grid.uncollateralized = - EXP.*FUND;
    XVA.uncollateralized      = sum(XVA_grid.uncollateralized);

else % FVA

    tmp_EXP                   = struct2cell(EXP);

    % Uncollateralized
    XVA_grid.uncollateralized = -SP_1.*SP_2.*tmp_EXP{1}'.*FUND';
    XVA.uncollateralized      = sum(XVA_grid.uncollateralized);

    % w VM
    XVA_grid.wVM              = -SP_1.*SP_2.*tmp_EXP{2}'.*FUND';
    XVA.wVM                   = sum(XVA_grid.wVM);

    % w IM
    XVA_grid.wIM              = -SP_1.*SP_2.*tmp_EXP{3}'.*FUND';
    XVA.wIM                   = sum(XVA_grid.wIM);

    % w VM + IM
    XVA_grid.wVM_IM           = -SP_1.*SP_2.*tmp_EXP{4}'.*FUND';
    XVA.wVM_IM                = sum(XVA_grid.wVM_IM);
end
end
