function [Pars_HW2F_new, Swapt_Cube, resnorm]=Calib_HW2F_Gamma_cube(Pars,DF,FWD,Pricing_Date,Expiry,Tenor_Swap,Cube_Expiry,Cube_Tenor,Cube_Strike,ATM_ATM,ATM_Cube,Swapt_Black,Swapt_Cube_Black,ATM_Weight,Cube_Weight,Cube_Omega);

% Carico Vacanze Calendario Target
[d1,d2,d3] = textread('TargetDates.txt','%s %s %s');
dv = strcat(d1,d2,d3);
dv=datenum(dv);
% caricato il file delle ferie

% [N,M]=size(Swapt_Black);
% for i=1:M
% F_comp(1+(i-1)*N:i*N,1)=Swapt_Black(:,i);
% end

F_comp = Swapt_Cube_Black;

Pars_Std=Pars(1:5,1)';
T_Gamma=Pars(6:end,1);
Gamma=Pars(6:end,2);
LGamma=Pars(6:end,3);
UGamma=Pars(6:end,4);
NT=size(T_Gamma,1);


[yy,mm,dd]  = datevec(Pricing_Date);
Expiry_Swpt = busdate(datenum([yy*ones(length(Cube_Expiry),1) mm+Cube_Expiry dd*ones(length(Cube_Expiry),1)])-1,1,dv);
Tmat = (Expiry_Swpt-Pricing_Date)/365;

% nexp=size(Expiry,1);
% nten_Swap=size(Tenor_Swap,1);


for i=1:NT
    W = zeros(length(Cube_Expiry),1);    
    % weights only at i Expiry
    idx    = find(Tmat == T_Gamma(i));
    W(idx) = 1;
    
    x=Gamma(i,1);
    LB=LGamma(i,1);
    UB=UGamma(i,1);
    
    options= optimset('TolX',1e-8, 'TolFun',1e-8,'Display','iter','MaxFunEvals',10000,'MaxIter',15);
    perc=1;
    Pars_old=[Pars_Std;[T_Gamma Gamma ones(NT,3)]];
    i
    [T_Gamma Gamma]
    [x,resnorm(:,i)] = lsqnonlin(@(x) Diff_PriceHW2F_Gamma_cube(x,Pars_old,T_Gamma(i),DF,FWD,Pricing_Date,0,0,Cube_Expiry,Cube_Tenor,Cube_Strike,0,Cube_Omega,F_comp,W,dv,perc),x,LB,UB,options);
    Gamma(i,1)=x;
    i
    [T_Gamma Gamma]
end

Pars_HW2F_new=[Pars_Std;[T_Gamma Gamma ones(NT,3)]];

AA=Diff_PriceHW2F_Gamma_cube(Gamma(end,1),Pars_HW2F_new,T_Gamma(end,1),DF,FWD,Pricing_Date,0,0,Cube_Expiry,Cube_Tenor,Cube_Strike,0,Cube_Omega,F_comp,W,dv,2);
AA=AA+F_comp;

Swapt_Cube=AA;

