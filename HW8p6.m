tic()
clc
close all
clear

%% preliminaries
colors
dt = 0.01;
Q = 0.1*dt;
R = 1;
spinsteps = 100000;
nsteps = 20000;
n = 10;
Traj = zeros(n,nsteps);
PF = zeros(n,nsteps);
Ne = 400;
w = zeros(1,Ne);
color1 = 11;
color2 = 9;
color3 = 22;
color4 = 8;
coords1 = [-9 12 -9 12 -9 12];
name = 'SPF3';
ll = 25;
xx = 1;
yy = 2;
zz = 3;
[M,N,H] = makeL96matrices(n);
F = 8;
szh = size(H,1);
szq = size(H,2);
QM = Q.*eye(szq);
sqrtQ = sqrtm(QM);
RM = R.*eye(szh);
%%

%% initial spin up
Wallaby = randn(n,1);

for ii=2:nsteps
    Wallaby = lorenz96s4(Wallaby,dt,M,N,F)+sqrtQ*randn(n,1);
end
%%

%% true state
Traj(:,1) = Wallaby;
for ii=2:nsteps
    Traj(:,ii) = lorenz96s4(Traj(:,ii-1),dt,M,N,F)+sqrtQ*randn(n,1);
end
%%

%% observations
Obs = H*Traj +R.*randn(n/2,nsteps);
%%

%% ensemble + simulation
Ens = ensemble_init(dt,Ne,M,N,F,Wallaby);

for jj=2:nsteps
    Ens = lorenz96s4(Ens,dt,M,N,F) + sqrtQ*randn(n,Ne);
    for kk=1:Ne
        w(kk) = 0.5.*(Obs(:,jj)-H*Ens(:,kk))'*(eye(n/2)./R)*(Obs(:,jj)-H*Ens(:,kk));
    end
    W = normalizeweights(w);
    Ens = resamplingmmo(W,Ens,Ne,n);
    PF(:,jj) = mean(Ens,2);
end
%%


%% plot

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

%lorenz96plots2(Traj,color4,coords1,name,nsteps,xx,yy,zz)
%lorenz96plots(Traj,PF,color1,color2,color3,coords1,name,ll,nsteps,xx,yy,zz)
%%

toc()