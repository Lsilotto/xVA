function AA=Diff_PriceHW2F_Gamma_cube(x,Pars_in,T_GammaI,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,Cube_Omega,F_comp,W,dv,perc);





%addpath(genpath('../'));
tau=(DF(2:end,1))/365;
Rf=-log(DF(2:end,2))./tau;
Rf=[0;Rf];

tauf=(FWD(2:end,1))/365;
Rfwd=-log(FWD(2:end,2))./tauf;
Rfwd=[0;Rfwd];

% mettere vettore
%PayRec=0;
Method=1;
Pars_Std=Pars_in(1,:);
T_Gamma=Pars_in(2:end,1);
Gamma=Pars_in(2:end,2);
index_Gamma=find(T_GammaI==T_Gamma);

% nexp=size(Expiry,1);
% nten_Swap=size(Tenor_Swap,1);

nexp_Cube=size(Expiry_Cube,1);

[yy,mm,dd]=datevec(Pricing_Date);

AA=zeros(nexp_Cube,1);

Gamma(index_Gamma)=x;
Pars=[Pars_Std;[T_Gamma Gamma ones(length(unique(Expiry_Cube)),3)]];

% Pricing Swaption ATM
dcbasis_Fix=6;
dcbasis_Float=2;


for i=1:nexp_Cube
    Expiry_Swaption=busdate(datenum([yy mm+Expiry_Cube(i) dd])-1,1,dv);
    Start_DateX=busdate(busdate(Expiry_Swaption,1,dv),1,dv);
    
    Coupon_Index_S=Swap_Date_Generator(Start_DateX,Tenor_Cube(i),12,6,dv);
    % Parte Fix
    Swap_Dates=Coupon_Index_S(Coupon_Index_S(:,2)==1);
    Swap_Dates=[Expiry_Swaption;Start_DateX;Swap_Dates];
    % Serve anche l'expiry per il pricing corretto
    % tau_fixx=(Swap_Dates-Pricing_Date)/365;
    % DF_M=exp(-interp1(Pricing_Date+DF(:,1),Rf,Swap_Dates,'linear').*tau_fixx);
    % DFwd_M=exp(-interp1(Pricing_Date+FWD(:,1),Rfwd,Swap_Dates,'linear').*tau_fixx);
    Strike=Strike_Cube(i);
    Swap_Dates=Swap_Dates(2:end);
    
    % Parte Float
    Swap_Dates_Float=Coupon_Index_S(Coupon_Index_S(:,3)==1);
    Swap_Dates_Float=[Expiry_Swaption;Start_DateX;Swap_Dates_Float];
    % Serve anche l'expiry per il pricing corretto
    tau_float=(Swap_Dates_Float-Pricing_Date)/365;
    DF_FLOAT=exp(-interp1(Pricing_Date+DF(:,1),Rf,Swap_Dates_Float,'linear').*tau_float);
    Fwd_FLOAT=exp(-interp1(Pricing_Date+FWD(:,1),Rfwd,Swap_Dates_Float,'linear').*tau_float);
    Swap_Dates_Float=Swap_Dates_Float(2:end);
    
    [Price_P,Price_R]=Price_Swaption_HW2F_Gamma_MULTI(Pricing_Date,Expiry_Swaption,Pars,Swap_Dates,Swap_Dates_Float,DF_FLOAT,Fwd_FLOAT,...
        Strike,dcbasis_Fix,0,Method);
    
    if Cube_Omega(i) == 1
        Price = Price_P;
    else
        Price = Price_R;
    end
    discount_expiry = exp(-interp1(Pricing_Date+DF(:,1),Rf,Expiry_Swaption,'linear').*(Expiry_Swaption-Pricing_Date)/365);
    AA(i)=Price/discount_expiry;
    
    
end


if perc==1
    AA=((AA+0.001)./(F_comp+0.001)-1).*W/sum(W);
elseif perc==0
    AA=(AA-F_comp).*W/sum(W);
else
    AA=AA-F_comp;
end

