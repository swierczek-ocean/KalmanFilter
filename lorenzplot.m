function [M,N,H,SynthDataTrue,SynthDataObs,X_start] = lorenzplot(n,dt,t_final,F,R,jump,j,FR)
tic();

v1 = ones(n-2,1);
v2 = ones(n-1,1);
M = diag(v2,1) - diag(v1,-2);
M(n,1)=1; M(1,n-1)=-1; M(2,n)=-1;
N=diag(v2,-1);
N(1,n)=1;

X = unifrnd(-1,1,n,1);
init_iter = ceil(100/dt);
numiter = ceil(t_final/dt);

f = @(x)(M*x.*(N*x)+(F-x));

% First, find an initial state by running ODE for a while

for i=1:init_iter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
end

xmin = min(X)-1;
xmax = max(X)+1;

X_start=X;

figure, set(gcf, 'Color','white')
plot3(X_start(j),X_start(j+1),X_start(j+2))
axis([xmin xmax xmin xmax xmin xmax])
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz.avi');
vidObj.Quality = 100;
vidObj.FrameRate = FR;
open(vidObj);
writeVideo(vidObj, getframe(gca));


% Next, create synthetic data set from initial state

H = zeros(n/2,n);

for i=1:n/2
   H(i,2*i-1)=1; 
end

SynthDataTrue = zeros(n,numiter);
SynthDataTrue(:,1)=X;
SynthDataObs = zeros(n/2,numiter);

for i=2:numiter
   k1=f(X);
   k2=f(X+0.5*dt.*k1);
   k3=f(X+0.5*dt.*k2);
   k4=f(X+dt.*k3);
   X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
   SynthDataTrue(:,i)=X;
   
    if mod(i,jump)==0
        plot3(X(j),X(j+1),X(j+2))
        axis([xmin xmax xmin xmax xmin xmax])
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

for i=1:numiter
   SynthDataObs(:,i) = H*SynthDataTrue(:,i) + normrnd(0,R,n/2,1); 
end

toc()
end

