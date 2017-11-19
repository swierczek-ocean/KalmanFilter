function [x_start,x_star] = fdvar(x_0)
fun = @(x)argh(x);
options = optimoptions('lsqnonlin','SpecifyObjectiveGradient',true)
x_star = lsqnonlin(fun,x_0,options);
x_start = MRK2(x_star,dt,obsdt,jump);
end

