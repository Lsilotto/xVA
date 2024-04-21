function S_t=survival_probability(terms,hazard_rates)

terms=terms/365;
S_t=[];
hazard_terms=[0; hazard_rates(:,1)]/365;
hazard_rates(:,1)=[];
for i=1:size(terms,1)
    n=find(hazard_terms<=terms(i),1,'last');
    if n==length(hazard_terms)
        n=n-1;
    end
    tn=[0; cumsum(diff(hazard_terms).*hazard_rates)];
    S_t(i)=exp(-tn(n)-(terms(i)-hazard_terms(n,1))*hazard_rates(n));
end