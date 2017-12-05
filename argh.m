function [arrr,Jaco] = argh(x_0,H,bcov,bmean,y_t,n,R,L1,L2,F,dt,jump)

B = 0.5.*(bcov+bcov');
[x_out,TimeSeries] = MRK2(x_0,L1,L2,F,n,dt,jump);
arrr = [real(sqrtm(B))\(x_0-bmean);(H*x_out-y_t)./sqrt(R)];
Jaco = [real(sqrtm(B))\eye(n);H*Mlin3(TimeSeries,L1,L2,F,n,dt,jump)./sqrt(R)];
end

