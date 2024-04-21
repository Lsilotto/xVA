function AA=Diff_PriceHW2F(Pars,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,ATM_Cube,F_comp,W,dv,perc);

Pars'

tau=(DF(2:end,1))/365;
Rf=-log(DF(2:end,2))./tau;
Rf=[0;Rf];

tauf=(FWD(2:end,1))/365;
Rfwd=-log(FWD(2:end,2))./tauf;
Rfwd=[0;Rfwd];

PayRec=0;
Method=1;
Pars=Pars';
nexp=size(Expiry,1);
nten_Swap=size(Tenor_Swap,1);
Tenor_Float=6;
Tenor_Fix=12;
[yy mm dd]=datevec(Pricing_Date);

AA=zeros(nexp*(nten_Swap),1);

% Pricing Swaption ATM
dcbasis_Fix=6;
dcbasis_Float=2;

Method=1;
PayRec=0;
for i=1:nexp
    Expiry_Swaption=busdate(datenum([yy mm+Expiry(i) dd])-1,1,dv);
    Start_DateX=busdate(busdate(Expiry_Swaption,1,dv),1,dv);
	for j=1:nten_Swap
        Coupon_Index_S=Swap_Date_Generator(Start_DateX,Tenor_Swap(j),12,6,dv);
        % Fix side
        Swap_Dates=Coupon_Index_S(Coupon_Index_S(:,2)==1);
        Swap_Dates=[Expiry_Swaption;Start_DateX;Swap_Dates];
        Strike=ATM_ATM(i,j);
        Swap_Dates=Swap_Dates(2:end);

        % Float side
        Swap_Dates_Float=Coupon_Index_S(Coupon_Index_S(:,3)==1);
        Swap_Dates_Float=[Expiry_Swaption;Start_DateX;Swap_Dates_Float];
        tau_float=(Swap_Dates_Float-Pricing_Date)/365;
        DF_FLOAT=exp(-interp1(Pricing_Date+DF(:,1),Rf,Swap_Dates_Float,'linear').*tau_float);
        Fwd_FLOAT=exp(-interp1(Pricing_Date+FWD(:,1),Rfwd,Swap_Dates_Float,'linear').*tau_float);
        Swap_Dates_Float=Swap_Dates_Float(2:end);
  
        [Price_P,Price_R]=Price_Swaption_HW2F_MULTI(Pricing_Date,Expiry_Swaption,Pars,Swap_Dates,Swap_Dates_Float,DF_FLOAT,Fwd_FLOAT,...
        Strike,dcbasis_Fix,PayRec,Method);

        discount_expiry = exp(-interp1(Pricing_Date+DF(:,1),Rf,Expiry_Swaption,'linear').*(Expiry_Swaption-Pricing_Date)/365);
        AA(1+(i-1)+(j-1)*nexp,1)=(Price_P+Price_R)/discount_expiry;
	end
end


if perc==1
    AA=((AA+0.001)./(F_comp+0.001)-1);
elseif perc==0
    AA=(AA-F_comp);
else
    AA=AA-F_comp;
end

