function Coupon_Index=Swap_Date_Generator(Start_Date,Maturity,Tenor_Fix,Tenor_Float,dv,Convention)


if nargin == 5
Convention='modfol';
end

% following 'follow'
% modified following 'modfol'
% preceding 'preced'
% modified preceding 'modprec'
% actual 'actual'

% Input
% Start_Date
% Maturity in months
% Tenor_Fix
% Tenor_Float
% dv holidays

% Output
% Coupon_Index: 3 cols: col1 all dates, col2= 1 if fixed leg CF date, col3= 1 if float leg CF date
Start_DateX=Bus_Adj(Start_Date,Convention,dv);
[yy mm dd] = datevec(Start_DateX);

ncoup_Fix=round(Maturity/Tenor_Fix);
vec_Fix=[repmat(yy,ncoup_Fix+1,1) mm+(0:Tenor_Fix:ncoup_Fix*Tenor_Fix)' repmat(dd,ncoup_Fix+1,1)]; 
vec_Fix1=datevec(datenum(vec_Fix));

for i=1:ncoup_Fix
if vec_Fix(i,3)~=vec_Fix1(i,3)
  vec_Fix(i,:)=[yy mm+(i-1)*Tenor_Fix dd-1];
  vec_Fix1(i,:)=datevec(datenum(vec_Fix(i,:)));
  
  if vec_Fix(i,3)~=vec_Fix1(i,3)
  vec_Fix(i,:)=[yy mm+(i-1)*Tenor_Fix dd-2];
  vec_Fix1(i,:)=datevec(datenum(vec_Fix(i,:)));
  
  if vec_Fix(i,3)~=vec_Fix1(i,3)
  vec_Fix(i,:)=[yy mm+(i-1)*Tenor_Fix dd-3];
  vec_Fix1(i,:)=datevec(datenum(vec_Fix(i,:)));
  end
  
  end
end
end

Date_Fix=datenum(vec_Fix);
% Schedule gen - if Start End of month - all end of months
Date_Fix=Bus_Adj(Date_Fix,Convention,dv); % check holiday

% Date Float coupon
ncoup_Float=round(Maturity/Tenor_Float);
vec_Float=[repmat(yy,ncoup_Float+1,1) mm+(0:Tenor_Float:ncoup_Float*Tenor_Float)' repmat(dd,ncoup_Float+1,1)];
vec_Float1=datevec(datenum(vec_Float));
for i=1:ncoup_Float
if vec_Float(i,3)~=vec_Float1(i,3)
  vec_Float(i,:)=[yy mm+(i-1)*Tenor_Float dd-1];
  vec_Float1(i,:)=datevec(datenum(vec_Float(i,:)));
  
  if vec_Float(i,3)~=vec_Float1(i,3)
  vec_Float(i,:)=[yy mm+(i-1)*Tenor_Float dd-2];
  vec_Float1(i,:)=datevec(datenum(vec_Float(i,:)));
  
  if vec_Float(i,3)~=vec_Float1(i,3)
  vec_Float(i,:)=[yy mm+(i-1)*Tenor_Float dd-3];
  vec_Float1(i,:)=datevec(datenum(vec_Float(i,:)));
  end
  
  end
end
end
Date_Float=datenum(vec_Float);
Date_Float=Bus_Adj(Date_Float,Convention,dv);
% Date Floating Leg
Dates_Tot=unique([Date_Fix;Date_Float]);

ndates=size(Dates_Tot,1);
Coupon_Index=[Dates_Tot zeros(ndates,2)]; 
for i=1:ndates
index=find(Coupon_Index(i,1)==Date_Fix);
if isempty(index)==0
Coupon_Index(i,2)=1;
end

index=find(Coupon_Index(i,1)==Date_Float);
if isempty(index)==0
Coupon_Index(i,3)=1;
end
end

Coupon_Index(1,2:end)=0;

