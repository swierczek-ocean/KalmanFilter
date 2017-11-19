function [arrr,Jaco] = argh(x_0)
global H, global bcov, global bmean, global y_t
global n, global R
B = 0.5.*(bcov+bcov');
[x_out,TimeSeries] = MRK2(x_0);
arrr = [sqrtm(B)\(x_0-bmean);(H*x_out-y_t)./sqrt(R)];
Jaco = [sqrtm(B)\eye(n);H*Mlin(TimeSeries)./sqrt(R)];
end

