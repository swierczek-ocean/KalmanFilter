function X = lorenz63s1(x,dt,M1,M2,M3)
f = @(x)((M1*x).*(M2*x) + M3*x);

k1 = f(x);
k2 = f(x+0.5*dt.*k1);
k3 = f(x+2*dt.*k2 -dt.*k1);
X = x + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);

end

