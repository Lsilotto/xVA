function[r_shock] = r_shock(r)

n_pillar = length(r);

for i = 1:n_pillar
    r_shock(i,:) = [r(1:i-1) r(i)+0.0001 r(i+1:end)];
end
