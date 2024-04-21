function [Price_P,Price_R]=Price_Swaption_HW2F_Gamma_MULTI(Pricing_Date,Expiry_Date,Pars,Swap_Dates,Swap_Dates_Float,DF_FLOAT,Fwd_FLOAT,Strike,dcbasis_fix,PayRec,Method)


% size
nK=size(Strike,1); % par swap rate
nT=size(Swap_Dates_Float,1)-1; % number of cash flows

% Parameters
a=Pars(1,1);
sigma=Pars(1,2);
b=Pars(1,3);
eta=Pars(1,4);
rho=Pars(1,5);
Ti=Pars(2:end,1); % expiries
Gamma_i=Pars(2:end,2); % gamma
nGamma=size(Ti,1);

phid = DF_FLOAT(2:end)./DF_FLOAT(1:end-1);
phif = Fwd_FLOAT(2:end)./Fwd_FLOAT(1:end-1);

psi = phid./phif;

taux=yearfrac(Swap_Dates(1:end-1),Swap_Dates(2:end),dcbasis_fix);
tau=[ ];
for i=1:size(Swap_Dates,1)-1
    tau=[tau;0;Strike*taux(i)];
end

psi = psi(2:end);
c = tau+1;
c(1:end-1) = c(1:end-1) - psi(2:end);
c = c/psi(1);

% times
t=(Expiry_Date-Pricing_Date)/365;
T=(Swap_Dates_Float-Pricing_Date)/365;
% t expiry, T swap cash flows

AtT=AtTX_HW2F_Gamma(Pars,t,T);
AtT=DF_FLOAT(2:end,1)./DF_FLOAT(1,1).*AtT;
Ba=BtT_HW2F(a,t,T);
Bb=BtT_HW2F(b,t,T);

% pag 147 Brigo-Mercurio
Var_X=VarX_HW2F_Gamma([[a sigma];[Ti Gamma_i]],0,t); 
Var_Y=VarY_HW2F_Gamma([[b eta];[Ti Gamma_i]],0,t);
Covar_X_Y=VarX_HW2F_Gamma([[(a+b)/2 sqrt(sigma*eta)];[Ti Gamma_i]],0,t);


sigx=sqrt(Var_X);
sigy=sqrt(Var_Y);
rhoxy=rho*Covar_X_Y./(sigx*sigy);

[Mx,My]=MT_xy_HW2_Gamma(Pars,t,0,t);
mux=-Mx;
muy=-My;

Price_P=zeros(1,nK);
Price_R=Price_P;

if Method==1
% Numerical integration (rectangular rule)
nsigma = 4;
nsteps = 10;
liminf = mux - sqrt(2)*nsigma*sigx;
limsup = mux + sqrt(2)*nsigma*sigx;
dx  = (limsup-liminf)/(nsteps);
ii=(0:1:nsteps);
x=liminf+ii*dx;
[X, Y]=HW2F_SwaptionIntegrand(x,mux,muy,sigx,sigy,rhoxy,AtT,Ba,Bb,c,PayRec);
if PayRec==0
Price_P=sum(X,1)*dx.*DF_FLOAT(1,1);
Price_R=sum(Y,1)*dx.*DF_FLOAT(1,1);
else
% Y = 0
index=find(PayRec==1);
Price_P(index)=sum(X(:,index),1)*DF_FLOAT(1,1)*dx;
index1=find(PayRec==-1);
Price_R(index1)=sum(X(:,index1),1)*DF_FLOAT(1,1)*dx;
end
elseif Method==2
nsigma=10;
tolf=1.e-10;
liminf = mux - nsigma*sigx;
limsup = mux + nsigma*sigx;
Q=quadv(@HW2F_SwaptionIntegrand,liminf,limsup,tolf,[],mux,muy,sigx,sigy,rhoxy,AtT,Ba,Bb,c,PayRec);

QQ = @(y)quadv(@HW2F_SwaptionIntegrand,liminf,limsup,y,[],mux,muy,sigx,sigy,rhoxy,AtT,Ba,Bb,c,PayRec);

for i = 1:10
tic
value(i) = QQ(10^-i);
time(i) = toc;
end




yyaxis left
plot(value)
yyaxis right
plot(time)





Price_P=Q(:,1).*DF_FLOAT(1,1);
Price_R=Q(:,2).*DF_FLOAT(1,1);
end


Price_P = Price_P*psi(1);
Price_R = Price_R*psi(1);
