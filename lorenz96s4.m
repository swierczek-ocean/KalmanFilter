function X = lorenz96s4(x,dt,M,N,F)
f = @(x)(M*x.*(N*x)+(F-x));

k1 = f(x);
k2 = f(x+0.5*dt.*k1);
k3 = f(x+2*dt.*k2 -dt.*k1);
X = x + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);

end

