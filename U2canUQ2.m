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
Ne = 100;
color1 = 13;
color2 = 9;
color3 = 6;
[M,N,H] = makeL96matrices(n);
F = 8;
alpha = 0;
r = 5;
%L = localize2(n,r);
L = ones(n,n);
Rm = R*eye(nobs);
EnKFest = zeros(n,nsteps/jump);
spreadv = zeros(1,nsteps/jump);
Error = zeros(nexp,nsteps/jump);
RMSE = zeros(nexp,1);
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

for zz=1:nexp
    zz
%% ensemble + simulation
Ens = ensemble_init(dt,Ne,M,N,F,Wallaby);
EnKFest(:,1) = mean(Ens,2);
obssteps = [1];

for jj=2:(nsteps/jump)
    [Ens,mu_a,spread] = ENKFSRU(Ens,dt,jump,n,Ne,H,M,N,F,alpha,Obs(:,(jj-1)*jump+1),L,Rm,nobs);
    EnKFest(:,jj) = mu_a;
    spreadv(jj) = spread;
    obssteps = [obssteps,(jj-1)*jump+1];
end
%%

%% RMSE
Error_temp = Traj(:,obssteps) - EnKFest;
Error(zz,:) = sqrt(sum(Error_temp.^2))./sqrt(n);
st = 20;
RMSE(zz) = mean(Error(zz,st:end));
fprintf('RMSE %g\n',RMSE(zz))
fprintf('Spread %g\n',mean(spreadv(st:end)))

%%
end

%% plots
sk = 10;
figure
h2 = plot(sk:(nsteps/jump),Error(1,sk:end),'Color',Color(:,color3),'Linewidth',1.5);
hold on
for ll=2:nexp
   plot(sk:(nsteps/jump),Error(ll,sk:end),'Color',Color(:,color3),'Linewidth',1.5) 
end
h1 = plot(sk:(nsteps/jump),spreadv(sk:end),'Color',Color(:,color1),'Linewidth',2.5);
title('SqEnKF RR RMSE and spread')
xlabel('time steps')
ylabel('error')
legend([h1(1),h2(1)],'spread','RMSE')
print('SqEnKF_RR_RMSE','-djpeg')
hold off

%%


toc()