function [argh,Jaco] = argh(x_0)
B = 0.5.*(B+B');
[U,lambda] = eig(B);
Q = U*sqrt(lambda);
[x_out,TimeSeries] = MRK2(x,dt,obsdt,jump)
argh = [Q*(x_0-mean);(H*x_out-y_t)./sqrt(R)];
Jaco = [Q;H*Mlin(TimeSeries,dt,obsdt,jump)./sqrt(R)];
end

