function[r_tT] = zcb2r(P_tT,tau_tT)

if verLessThan ('matlab','9.4') == 0
    r_tT = -log(P_tT)./tau_tT;
else
    r_tT = bsxfun(@rdivide,-log(P_tT),tau_tT);
end

r_tT(isnan(r_tT)) = 0;