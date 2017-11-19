function ensemble = ensemble_init3(dt,X_start)
tic();
num_iter = 1000/dt;

X = X_start;
f = @(x)(L1*x.*(L2*x)+(F-x));
ensemble = [];

for i=1:num_iter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
   ensemble = [ensemble,X];
end

toc()
end
