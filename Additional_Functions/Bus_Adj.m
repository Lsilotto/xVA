function Date_Out=Bus_Adj(Date_In,Convention,dv)

% dv : holidays
% Convention
% following 'follow'
% modified following 'modfol'
% preceding 'preced'
% modified preceding 'modprec'
% actual 'actual'
% Date_in: date input

N_Dates=size(Date_In,1);
Date_Out=zeros(N_Dates,1);

for i=1:N_Dates

if isbusday(Date_In(i))
Date_Out(i)=Date_In(i);
else
	if Convention=='follow'
	Date_Out(i)=busdate(Date_In(i),1,dv);
	elseif Convention=='modfol'
	Date_Out(i)=busdate(Date_In(i),1,dv);
	[yy mm dd]=datevec(Date_In(i));
	[yy1 mm1 dd1]=datevec(Date_Out(i));
		if mm1 ~= mm
		% changed month
		Date_Out(i)=busdate(Date_In(i),-1,dv);
		end
	elseif Convention=='preced'	
	Date_Out(i)=busdate(Date_In(i),-1,dv);
	elseif Convention=='modpre'
	Date_Out(i)=busdate(Date_In(i),-1,dv);
	[yy mm dd]=datevec(Date_In(i));
	[yy1 mm1 dd1]=datevec(Date_Out(i));
		if mm1 ~= mm
		% changed month
		Date_Out(i)=busdate(Date_In(i),1,dv);
		end
	end
end

end

