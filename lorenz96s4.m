function X = lorenz96s4(x,dt,M,N,F)
f = @(x)(M*x.*(N*x)+(F-x));

k1=f(x);
k2=f(x+0.5*dt.*k1);
k3=f(x+0.5*dt.*k2);
k4=f(x+dt.*k3);
X = x + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);

end

