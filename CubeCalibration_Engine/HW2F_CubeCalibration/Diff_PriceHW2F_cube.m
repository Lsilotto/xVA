function AA=Diff_PriceHW2F_cube(Pars,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Cube_Strike,ATM_ATM,Cube_Omega,F_comp,W,dv,perc);

Pars'

tau=(DF(2:end,1))/365;
Rf=-log(DF(2:end,2))./tau;
Rf=[0;Rf];

tauf=(FWD(2:end,1))/365;
Rfwd=-log(FWD(2:end,2))./tauf;
Rfwd=[0;Rfwd];

Method=1;
Pars=Pars';

nexp_Cube=size(Expiry_Cube,1);

Tenor_Float=6;
Tenor_Fix=12;
[yy mm dd]=datevec(Pricing_Date);

AA=zeros(nexp_Cube,1);

% Pricing Swaption ATM
dcbasis_Fix=6;
dcbasis_Float=2;

Method=1;
%PayRec=0;
for i=1:nexp_Cube
    Expiry_Swaption=busdate(datenum([yy mm+Expiry_Cube(i) dd])-1,1,dv);
    Start_DateX=busdate(busdate(Expiry_Swaption,1,dv),1,dv);
    
    Coupon_Index_S=Swap_Date_Generator(Start_DateX,Tenor_Cube(i),12,6,dv);
    % Fix side
    Swap_Dates=Coupon_Index_S(Coupon_Index_S(:,2)==1);
    Swap_Dates=[Expiry_Swaption;Start_DateX;Swap_Dates];
    Strike=Cube_Strike(i);
    Swap_Dates=Swap_Dates(2:end);
    
    % Float side
    Swap_Dates_Float=Coupon_Index_S(Coupon_Index_S(:,3)==1);
    Swap_Dates_Float=[Expiry_Swaption;Start_DateX;Swap_Dates_Float];
    tau_float=(Swap_Dates_Float-Pricing_Date)/365;
    DF_FLOAT=exp(-interp1(Pricing_Date+DF(:,1),Rf,Swap_Dates_Float,'linear').*tau_float);
    Fwd_FLOAT=exp(-interp1(Pricing_Date+FWD(:,1),Rfwd,Swap_Dates_Float,'linear').*tau_float);
    Swap_Dates_Float=Swap_Dates_Float(2:end);
    
    [Price_P,Price_R]=Price_Swaption_HW2F_MULTI(Pricing_Date,Expiry_Swaption,Pars,Swap_Dates,Swap_Dates_Float,DF_FLOAT,Fwd_FLOAT,...
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
    AA=((AA+0.001)./(F_comp+0.001)-1);
elseif perc==0
    AA=(AA-F_comp);
else
    AA=AA-F_comp;
end

