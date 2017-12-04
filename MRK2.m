function [x_out,TimeSeries] = MRK2(x,L1,L2,F,n,dt,jump)

f = @(x)(L1*x.*(L2*x)+(F-x));
TimeSeries = [x,zeros(n,jump)];

for i=1:jump
    k = 0.5*dt.*f(TimeSeries(:,i)+dt.*f(TimeSeries(:,i)));
    TimeSeries(:,i+1) = TimeSeries(:,i) + 0.5*dt.*f(TimeSeries(:,i)) + k;
end

x_out = TimeSeries(:,end);

end

