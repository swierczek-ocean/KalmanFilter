function [x_out,TimeSeries] = MRK2(x)
global L1, global L2, global F, global n
global dt, global jump

f = @(x)(L1*x.*(L2*x)+(F-x));
TimeSeries = [x,zeros(n,jump)];

for i=1:jump
    k = 0.5*dt.*f(TimeSeries(:,i)+dt.*f(TimeSeries(:,i)));
    TimeSeries(:,i+1) = TimeSeries(:,i) + 0.5*dt.*f(TimeSeries(:,i)) + k;
end

x_out = TimeSeries(:,end);

end

