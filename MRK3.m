function [x_out,TimeSeries] = MRK3(x,L1,L2,F,n,dt,jump)

f = @(x)(L1*x.*(L2*x)+(F-x));
TimeSeries = [x,zeros(n,jump)];

for i=1:jump
    k1 = f(TimeSeries(:,i));
    k2 = f(TimeSeries(:,i)+0.5*dt.*k1);
    k3 = f(TimeSeries(:,i)+2*dt.*k2 -dt.*k1);
    TimeSeries(:,i+1) = TimeSeries(:,i) + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);
end

x_out = TimeSeries(:,end);

end

