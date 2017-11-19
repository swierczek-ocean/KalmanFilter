function [x_start,x_star,TimeSeries] = fdvar(x_0)
global dt, global obsdt, global jump, global lb,global ub
fun = @(x)argh(x);
% options = optimoptions(@lsqnonlin,'Display','iter-detailed',...
%     'SpecifyObjectiveGradient',true);
options = optimoptions(@lsqnonlin,'Display','iter-detailed',...
    'DerivativeCheck', 'on',...
    'SpecifyObjectiveGradient',true);
[x_star,~,~,~,~,~,J] = lsqnonlin(fun,x_0,[],[],options);
[x_start,TimeSeries] = MRK2(x_star);
end

