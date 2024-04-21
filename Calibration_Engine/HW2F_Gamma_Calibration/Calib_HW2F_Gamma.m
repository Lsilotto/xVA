function [Pars_HW2F_new, Swapt_ATM, Swapt_Cube, resnorm]=Calib_HW2F_Gamma(Pars,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Expiry_Cube,Tenor_Cube,Strike_Cube,ATM_ATM,ATM_Cube,Swapt_Black,Swapt_Cube_Black,ATM_Weight,Cube_Weight);

[d1,d2,d3] = textread('TargetDates.txt','%s %s %s'); 
dv = strcat(d1,d2,d3); 
dv=datenum(dv);

[N,M]=size(Swapt_Black);
for i=1:M
F_comp(1+(i-1)*N:i*N,1)=Swapt_Black(:,i);
end

Pars_Std=Pars(1:5,1)';
T_Gamma=Pars(6:end,1);
Gamma=Pars(6:end,2);
LGamma=Pars(6:end,3);
UGamma=Pars(6:end,4);
NT=size(T_Gamma,1);

nexp=size(Expiry,1);
nten_Swap=size(Tenor_Swap,1);


for i=1:NT
% weights only at i Expiry
ATM_WeightX=ones(nexp,nten_Swap);
ATM_WeightX(:,:)=0;
ATM_WeightX(i,:)=1;

for j=1:M
W(1+(j-1)*N:j*N,1)=ATM_WeightX(:,j);
end


x=Gamma(i,1);
LB=LGamma(i,1);
UB=UGamma(i,1);

options= optimset('TolX',1e-8, 'TolFun',1e-8,'Display','iter','MaxFunEvals',10000,'MaxIter',15);
perc=1;
Pars_old=[Pars_Std;[T_Gamma Gamma ones(NT,3)]];
i
[T_Gamma Gamma]
[x,resnorm(:,i)] = lsqnonlin(@(x) Diff_PriceHW2F_Gamma(x,Pars_old,T_Gamma(i),DF,FWD,Pricing_Date,Expiry,Tenor_Swap,0,0,0,ATM_ATM,0,F_comp,W,dv,perc),x,LB,UB,options);
Gamma(i,1)=x;
i
[T_Gamma Gamma]
end

Pars_HW2F_new=[Pars_Std;[T_Gamma Gamma ones(NT,3)]];

AA=Diff_PriceHW2F_Gamma(Gamma(end,1),Pars_HW2F_new,T_Gamma(end,1),DF,FWD,Pricing_Date,Expiry,Tenor_Swap,0,0,0,ATM_ATM,0,F_comp,W,dv,2);
AA=AA+F_comp;


Swapt_ATM=zeros(N,M);
for i=1:M
Swapt_ATM(:,i)=AA(1+(i-1)*N:i*N,1);
end

Swapt_Cube=0;

