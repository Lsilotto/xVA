function [fx,fy] = HW2F_SwaptionIntegrand( x,mux,muy,sigx,sigy,rhoxy,ATt,BaTt,BbTt,c,PayRec)

% Integrand function for swaption pricing
% Input: x row vector with x values
%	mux,muy: momens of x. y
%	sigx,sigy vola
%	rhoxy correl
%	ATt: vector nT+1 (with swap start date)
%	BatT BbtT coefficient
% 	c= vector of cash flows dimension NT by nStk
% PayRec= row vector nK : +1 payer, -1 Receiver
% the other is calculated via put_call parity
% if PayRec=0--> computed both payer and receiver

% pag 158 Brigo-Mercurio

nx=size(x,2);
nT=size(c,1);
nK=size(c,2);

ybar=zeros(nx,nK);
fvale=ybar;
flagge=ybar;
c=[ones(1,nK);c];
for j=1:nK
    y0=0;
    y=0;
    for i=1:nx
        y0=y;
        xx=x(i);
        [y,fval,exitflag]  = fzero(@(y) HW2F_Annuity(y,xx,c(:,j),ATt,BaTt,BbTt), y0);
        ybar(i,j) =y;
        fvale(i,j)=fval;
        flagge(i,j)=exitflag;
    end
end
%
nexc=size(find(flagge~=1),1)/(nK*nx); % percentage of grid points for which the solution was not found
err=max(max(abs(fvale)));
xbar=repmat(x',1,nK);

redrhoxy = sqrt(1-rhoxy*rhoxy);
h1x      = ( (ybar -muy)/sigy -rhoxy*(xbar-mux)/sigx )/redrhoxy;
h2x      = repmat(h1x,[1 1 nT+1]) + permute(repmat(BbTt,[1 nx nK]),[2 3 1])*sigy*redrhoxy; % compute h2 for all x and all t_i
lambdax  = permute(repmat((c.*repmat(ATt,1,nK)),[1 1 nx]),[3 2 1]).*exp(-permute(repmat(BaTt,[1 nx nK]),[2 3 1]).*repmat(x',[1 nK nT+1]));
kappax   = -permute(repmat(BbTt,[1 nx nK]),[2 3 1]).*(muy-0.5*redrhoxy*redrhoxy*sigy*sigy*permute(repmat(BbTt,[1 nx nK]),[2 3 1])+rhoxy*sigy*(repmat(x',[1 nK nT+1])-mux)/sigx);
NORMDIST=normpdf(x',mux,sigx);
NORMDIST=repmat(NORMDIST,1,nK);

if PayRec==0
    fx = normcdf(-h2x(:,:,1)).*lambdax(:,:,1).*exp(kappax(:,:,1))-sum(normcdf(-h2x(:,:,2:end)).*lambdax(:,:,2:end).*exp(kappax(:,:,2:end)),3);
    fy=-normcdf(h2x(:,:,1)).*lambdax(:,:,1).*exp(kappax(:,:,1))+sum(normcdf(h2x(:,:,2:end)).*lambdax(:,:,2:end).*exp(kappax(:,:,2:end)),3);
    fx=fx.*NORMDIST; % Normdist fisrt factor of integrand
    fy=fy.*NORMDIST;
else
    PayRec=-repmat(PayRec,[nx 1 nT+1]);
    fx=normcdf(PayRec(:,:,1).*h2x(:,:,1)).*lambdax(:,:,1).*exp(kappax(:,:,1))-sum(normcdf(PayRec(:,:,2:end).*h2x(:,:,2:end)).*lambdax(:,:,2:end).*exp(kappax(:,:,2:end)),3);
    fx=-PayRec(:,:,1).*NORMDIST.*fx;
    fy=zeros(nx,nK);
end
