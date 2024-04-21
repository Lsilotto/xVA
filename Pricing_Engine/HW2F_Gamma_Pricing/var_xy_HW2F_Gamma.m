function [var_x,var_y]=var_xy_HW2F_Gamma(Params,Gamma,dt)

a     = Params(1);
sigma = Params(2);
b     = Params(3);
eta   = Params(4);

var_x = var_HW2F_Gamma([a,sigma],Gamma,dt);
var_y = var_HW2F_Gamma([b,eta],Gamma,dt);
end