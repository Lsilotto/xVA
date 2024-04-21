function[x,y] = HW2F_fixedSeed_xy_simulation(Parameters,T,time_grid,random_numbers)

a      = Parameters(1,1);
b      = Parameters(1,3);
n_path = size(random_numbers{1},1)*2;

x      = zeros(n_path,length(time_grid));
y      = zeros(n_path,length(time_grid));

for i=1:length(time_grid)-1
    
    V  = zeros(n_path,2);
    AA = repelem(random_numbers{i},2,1);
    random_numbers{i} = [];
    
    V(mod([1:size(V,1)],2) == 1, :) = AA(mod([1:size(V,1)],2) == 1, :);
    V(mod([1:size(V,1)],2) == 0, :) = -AA(mod([1:size(V,1)],2) == 0, :);
    
    [Mx,My]  = MT_xy_HW2F_Gamma(Parameters,T,time_grid(i),time_grid(i+1));
        
    x(:,i+1) =  x(:,i)*exp(-a*(time_grid(i+1)-time_grid(i))) + Mx + V(:,1);
    y(:,i+1) =  y(:,i)*exp(-b*(time_grid(i+1)-time_grid(i))) + My + V(:,2);
    
end