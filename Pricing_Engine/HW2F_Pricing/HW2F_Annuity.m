function diff=HW2F_Annuity(y,x,C,A,Ba,Bb)

% Function for Jamshidian decomposition for swaptions HW 2F
% Start_Date:Expiry+2 
% c=cash flows
% A,B: 1st element Start Date, following elements swap cash flow dates

% pag 158 Brigo-Mercurio
nK=size(C,2);

diff=sum(C(2:end,:).*A(2:end,:).*exp(-Ba(2:end,:).*x-Bb(2:end,:).*y),1)-C(1,:).*A(1,:).*exp(-Ba(1,:).*x-Bb(1,:).*y);
