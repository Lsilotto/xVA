function [Pars_HW2F_new, Swapt_ATM, resnorm]=Calib_HW2F(Pars,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,ATM_Cube,Swapt_Black,Swapt_Cube_Black,ATM_Weight,Cube_Weight);

[d1 d2 d3] = textread('TargetDates.txt','%s %s %s'); 
dv = strcat(d1,d2,d3); 
dv=datenum(dv);

[N,M]=size(Swapt_Black);
F_comp=zeros(N*M,1);
for i=1:M
F_comp(1+(i-1)*N:i*N,1)=Swapt_Black(:,i);
end

x=Pars(:,1);
LB=Pars(:,2);
UB=Pars(:,3);
Pars_old=Pars;

options=optimset('TolX',1e-15, 'TolFun',1e-10, 'Display','iter','MaxFunEvals',10000,'MaxIter',30,'UseParallel','never');
perc=1;
[x, resnorm] = lsqnonlin(@(x) Diff_PriceHW2F(x,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,ATM_Cube,F_comp,0,dv,perc),x,LB,UB,options);
Pars_HW2F_new=x;

AA=Diff_PriceHW2F(x,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,ATM_Cube,F_comp,0,dv,2);
AA=AA+F_comp;

Swapt_ATM=zeros(N,M);
for i=1:M
Swapt_ATM(:,i)=AA(1+(i-1)*N:i*N,1);



end
