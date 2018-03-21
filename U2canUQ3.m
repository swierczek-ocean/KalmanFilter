tic()
clc
close all
clear

%% preliminaries
colors
dt = 0.01;
jump = 10;
R = 1;
spinsteps = 100000;
nsteps = 20000;
n = 40;
nexp = 20;
nobs = floor(n/2);
Traj = zeros(n,nsteps);
Ne = 300;
color1 = 13;
color2 = 9;
color3 = 6;
color4 = 8;
[M,N,H] = makeL96matrices(n);
F = 8;
alpha = 0;
r = 5;
%L = localize2(n,r);
L = ones(n,n);
Rm = R*eye(nobs);
EnKFestSR = zeros(n,nsteps/jump);
EnKFestNR = zeros(n,nsteps/jump);
EnKFestRR = zeros(n,nsteps/jump);
Error = zeros(nexp,nsteps/jump);
RMSE = zeros(nexp,1);
xx=1;
yy=2;
zz=3;
%%

%% initial spin up
Wallaby = randn(n,1);

for ii=2:spinsteps
    Wallaby = lorenz96s4(Wallaby,dt,M,N,F);
end
%%

%% true state
Traj(:,1) = Wallaby;
for ii=2:nsteps
    Traj(:,ii) = lorenz96s4(Traj(:,ii-1),dt,M,N,F);
end
%%

%% observations
Obs = H*Traj +sqrt(R).*randn(nobs,nsteps);
%%


    
%% ensemble + simulation
Ens = ensemble_init(dt,Ne,M,N,F,Wallaby);
EnKFestSR(:,1) = mean(Ens,2);
obssteps = [1];

for jj=2:(nsteps/jump)
    [Ens,mu_a,spread] = ENKFSR(Ens,dt,jump,n,Ne,H,M,N,F,alpha,Obs(:,(jj-1)*jump+1),L,Rm,nobs);
    EnKFestSR(:,jj) = mu_a;
    obssteps = [obssteps,(jj-1)*jump+1];
end
%%

Ens = ensemble_init(dt,Ne,M,N,F,Wallaby);
EnKFestNR(:,1) = mean(Ens,2);

for jj=2:(nsteps/jump)
    [Ens,mu_a,spread] = ENKFSRU(Ens,dt,jump,n,Ne,H,M,N,F,alpha,Obs(:,(jj-1)*jump+1),L,Rm,nobs);
    EnKFestNR(:,jj) = mu_a;
end
%%

Ne=100;
Ens = ensemble_init(dt,Ne,M,N,F,Wallaby);
EnKFestRR(:,1) = mean(Ens,2);

for jj=2:(nsteps/jump)
    [Ens,mu_a,spread] = ENKFSRU2(Ens,dt,jump,n,Ne,H,M,N,F,alpha,Obs(:,(jj-1)*jump+1),L,Rm,nobs);
    EnKFestRR(:,jj) = mu_a;
end
%%


%% RMSE
Error_temp = Traj(:,obssteps) - EnKFestSR;
Error = sqrt(sum(Error_temp.^2))./sqrt(n);
st = 20;
fprintf('SR RMSE %g\n',mean(Error(st:end)))
%%

%% RMSE
Error_temp = Traj(:,obssteps) - EnKFestNR;
Error = sqrt(sum(Error_temp.^2))./sqrt(n);
st = 20;
fprintf('NR RMSE %g\n',mean(Error(st:end)))
%%

%% RMSE
Error_temp = Traj(:,obssteps) - EnKFestRR;
Error = sqrt(sum(Error_temp.^2))./sqrt(n);
st = 20;
fprintf('RR RMSE %g\n',mean(Error(st:end)))
%%



toc()