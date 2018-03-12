function X = lorenz63s1(x,dt,M1,M2,M3)
f = @(x)((M1*x).*(M2*x) + M3*x);

X = x + dt.*f(x);
end

