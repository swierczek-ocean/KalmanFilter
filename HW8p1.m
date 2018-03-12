tic()
clc
close all
clear all

%% preliminaries
colors
sigma = 10;
beta = 8/3;
rho = 28;
dt = 0.01;
Q = .0001;
R = .0001;
nsteps = 10000;
Traj = zeros(3,nsteps);
PF = zeros(3,nsteps);
Ne = 100;
%%

%% simulation and plots

for ii=2:nsteps
    Traj(:,ii) = lorenz63s(Traj(:,ii-1),dt,sigma,beta,rho,Q);
end

figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,9),'Linewidth',1.5)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
print('3D_trajectory','-djpeg')

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,8),'Linewidth',1.5)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
print('2D_trajectory','-djpeg')
%%

%% observations
H = [1,0,0;0,0,1];
Obs = H*Traj +R.*randn(2,nsteps);
%%

%% ensemble
Ens = zeros(3,Ne);

for jj=2:nsteps
    Ens = lorenz63s(Ens,dt,sigma,beta,rho,Q);
    w = 0.5.*(Obs(:,jj)-H*Ens)'*(Obs(:,jj)-H*Ens);
    W = normalizeweights(w);
    Ens = resamplingmmo(W,Ens,Ne,3);
    PF(:,jj) = mean(Ens')';
end

%%

%% plot
figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,8),'Linewidth',1.5)
hold on
plot3(PF(3,:),PF(1,:),PF(2,:),'*','Color',Color(:,9),'MarkerSize',3)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
print('3D_trajectory_plus_SPF','-djpeg')
hold off

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,8),'Linewidth',1.5)
hold on
plot(PF(1,:),PF(3,:),'*','Color',Color(:,9),'MarkerSize',3)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
print('2D_trajectory_plus_SPF','-djpeg')
hold off

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);

figure
plot(Error,'Color',Color(:,22),'Linewidth',3)
xlabel('time step')
ylabel('error')
print('SPF_error','-djpeg')
%%













toc()