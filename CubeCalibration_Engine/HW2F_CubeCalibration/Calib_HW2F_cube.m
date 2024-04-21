function [Pars_HW2F_new, Swapt_Cube,resnorm]=Calib_HW2F_cube(Pars,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,ATM_Cube,Swapt_Black,Swapt_Cube_Black,ATM_Weight,Cube_Weight,Cube_Omega);

[d1 d2 d3] = textread('TargetDates.txt','%s %s %s');
dv = strcat(d1,d2,d3);
dv=datenum(dv);


F_comp = Swapt_Cube_Black;


x=Pars(:,1);
LB=Pars(:,2);
UB=Pars(:,3);
Pars_old=Pars;

options=optimset('TolX',1e-15, 'TolFun',1e-10, 'Display','iter','MaxFunEvals',10000,'MaxIter',30,'UseParallel','never');
perc=1;
[x,resnorm] = lsqnonlin(@(x) Diff_PriceHW2F_cube(x,DF,FWD,Pricing_Date,0,0,Expiry_Cube,Tenor_Cube,Strike_Cube,0,Cube_Omega,F_comp,0,dv,perc),x,LB,UB,options);
Pars_HW2F_new=x;

AA=Diff_PriceHW2F_cube(x,DF,FWD,Pricing_Date,0,0,Expiry_Cube,Tenor_Cube,Strike_Cube,0,Cube_Omega,F_comp,0,dv,2);
AA=AA+F_comp;

Swapt_Cube = AA;

