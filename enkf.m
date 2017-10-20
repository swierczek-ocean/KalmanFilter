function [] = enkfpo(dt,ensemble_size)
tic();

f = @(x)(M*x.*(N*x)+(8-x));

% First, make an initial ensemble by generating a lot of 

for i=1:init_iter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
end

% Next, create synthetic data set from initial state

for i=2:numiter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
   SynthData(:,i)=X;
end 






toc()
end

