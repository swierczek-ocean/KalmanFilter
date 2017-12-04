function [SynthDataTrue,SynthDataObs,X_start] = lorenz3(n,t_final,L1,L2,H,F,dt,jump,R)

X = normrnd(0,1,n,1);
init_iter = ceil(2000/dt);
numiter = ceil(t_final/dt);
w = floor(numiter/jump);

f = @(x)(L1*x.*(L2*x)+(F-x));

% First, find an initial state by running ODE for a while

for i=1:init_iter
    k1 = f(X);
    k2 = f(X+0.5*dt.*k1);
    k3 = f(X+2*dt.*k2 -dt.*k1);
    X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);
end

X_start=X;

% Next, create synthetic data set from initial state

SynthDataTrue = zeros(n,numiter);
SynthDataTrue(:,1)=X;
SynthDataObs = zeros(n/2,w);

for i=2:numiter
    k1 = f(X);
    k2 = f(X+0.5*dt.*k1);
    k3 = f(X+2*dt.*k2 -dt.*k1);
    X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);
    SynthDataTrue(:,i)=X;
end

for i=1:floor(numiter/jump)
    SynthDataObs(:,i) = H*SynthDataTrue(:,jump*(i-1)+1) + normrnd(0,R,n/2,1);
end

end

