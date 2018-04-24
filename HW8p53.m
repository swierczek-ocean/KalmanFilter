tic()
clc
close all
clear

%% preliminaries
colors
n = 3;
nobs = 2;
r = 5;
L = localize2(n,r);
alpha = 0.15;
sigma = 10;
beta = 8/3;
rho = 28;
dt = 0.01;
R = 1;
Rm = R*eye(2);
Q = 0.00001*dt*eye(3);
sqrtQ = sqrtm(Q);
nsteps = 2600;
Traj = zeros(3,nsteps);
Traj(:,1) = [3;3;3];
KF = zeros(3,nsteps);
Ne = 20;
M1 = [0,0,0;-1,0,0;1,0,0];
M2 = [0,0,0;0,0,1;0,1,0];
M3 = [-sigma,sigma,0;rho,-1,0;0,0,-beta];
color1 = 11;
color2 = 12;
color3 = 16;
color4 = 31;
color5 = 33;
coords1 = [-10 60 -25 25 -30 30];
coords2 = [-21 21 -3 52];
name = 'ENKFPO2';
ll = 30;
%%

%% true state

for ii=2:nsteps
    Traj(:,ii) = lorenz63s4(Traj(:,ii-1),dt,M1,M2,M3)+sqrtQ*randn(3,1);
end
%%

%% observations
H = [1,0,0;0,0,1];
Obs = H*Traj +R.*randn(2,nsteps);
%%

%% ensemble + simulation + movie
Ens = randn(3,Ne);

for jj=2:nsteps
    Ens = lorenz63s4(Ens,dt,M1,M2,M3) + sqrt(dt).*randn(3,Ne);
    if(mod(jj,10)==0)
        [Ens,mu_a,spread] = LETKF(Ens,dt,n,Ne,H,M1,M2,M3,alpha,Obs(:,jj),L,Rm,nobs,sqrtQ);
    end
    KF(:,jj) = mean(Ens,2);
end
%%

%% plot

Error = Traj - KF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

% lorenz63plots2(Traj,KF,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
lorenz63plots5(Traj,KF,Obs,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
%%

toc()