tic()
clc
close all
clear

%% preliminaries
colors
n = 3;
r = 4;
alpha = 0.1;
sigma = 10;
beta = 8/3;
rho = 28;
dt = 0.01;
Q = 0.000001*dt;
R = 1;
nsteps = 20000;
Traj = zeros(3,nsteps);
KF = zeros(3,nsteps);
Ne = 20;
M1 = [0,0,0;-1,0,0;1,0,0];
M2 = [0,0,0;0,0,1;0,1,0];
M3 = [-sigma,sigma,0;rho,-1,0;0,0,-beta];
color1 = 11;
color2 = 12;
color3 = 22;
color4 = 11;
color5 = 12;
coords1 = [-10 60 -25 25 -30 30];
coords2 = [-25 25 -10 60];
name = 'ENKFPO1';
ll = 25;
%%

%% true state

for ii=2:nsteps
    Traj(:,ii) = lorenz63s4(Traj(:,ii-1),dt,M1,M2,M3)+sqrt(dt).*randn(3,1);
end
%%


%% observations
H = [1,0,0;0,0,1];
Obs = H*Traj +R.*randn(2,nsteps);
%%

%% ensemble + simulation
Ens = randn(3,Ne);

for jj=2:nsteps
    [Ens,mu_a,P_a] = ENKFPO2(Ens,dt,n,Ne,H,R,r,alpha,Obs(:,jj),M1,M2,M3,Q);
    KF(:,jj) = mean(Ens,2);
end
%%


%% plot

Error = Traj - KF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

lorenz63plots2(Traj,KF,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
%lorenz63plots3(Traj,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
%%

toc()