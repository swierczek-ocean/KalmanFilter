function [M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(n,dt,t_final,F,R,obsdt)
tic();

v1 = ones(n-2,1);
v2 = ones(n-1,1);
M = diag(v2,1) - diag(v1,-2);
M(n,1)=1; M(1,n-1)=-1; M(2,n)=-1;
N=diag(v2,-1);
N(1,n)=1;
jump = ceil(obsdt/dt);

X = unifrnd(-1,1,n,1);
init_iter = ceil(100/dt);
numiter = ceil(t_final/dt);
w = floor(numiter/jump);

f = @(x)(M*x.*(N*x)+(F-x));

% First, find an initial state by running ODE for a while

for i=1:init_iter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
end

X_start=X;

% Next, create synthetic data set from initial state

H = zeros(n/2,n);

for i=1:n/2
   H(i,2*i-1)=1; 
end

SynthDataTrue = zeros(n,numiter);
SynthDataTrue(:,1)=X;
SynthDataObs = zeros(n/2,w);

for i=2:numiter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
   SynthDataTrue(:,i)=X;
end 

for i=1:floor(numiter/jump)
   SynthDataObs(:,i) = H*SynthDataTrue(:,jump*(i-1)+1) + normrnd(0,R,n/2,1); 
end

toc()
end

