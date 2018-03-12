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
Q = 1;
R = 1;
nsteps = 20000;
Traj = zeros(3,nsteps);
KF = zeros(3,nsteps);
Ne = 20;
M1 = [0,0,0;-1,0,0;1,0,0];
M2 = [0,0,0;0,0,1;0,1,0];
M3 = [-sigma,sigma,0;rho,-1,0;0,0,-beta];
%%

%% simulation and plots

for ii=2:nsteps
    Traj(:,ii) = lorenz63s4(Traj(:,ii-1),dt,M1,M2,M3)+sqrt(dt)*sqrt(Q).*randn(3,1);
end

figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,9),'Linewidth',1.5)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
print('3D_trajectory_ENKFPO','-djpeg')

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.5)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
print('2D_trajectory_ENKFPO','-djpeg')
%%


%% observations
H = [1,0,0;0,0,1];
Obs = H*Traj +R.*randn(2,nsteps);
%%

%% ensemble + simulation + movie
Ens = randn(3,Ne);
szh = size(H,1);
szq = size(H,2);
QM = Q.*eye(szq);
RM = R.*eye(szh);
Z = eye(szh)/(H*QM*H'+R);
K = QM*H'*Z;

for jj=2:nsteps
    [Ens,mu_a,P_a] = ENKFPO2(Ens,dt,n,Ne,H,R,r,alpha,Obs(:,jj),M1,M2,M3,Q);
    KF(:,jj) = mean(Ens,2);
end


%%


%% plot
figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot3(KF(3,:),KF(1,:),KF(2,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
print('3D_trajectory_plus_ENKFPO','-djpeg')
hold off

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot(KF(1,:),KF(3,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
print('2D_trajectory_plus_ENKFPO','-djpeg')
hold off

Error = Traj - KF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,22),'MarkerSize',3)
xlabel('time step')
ylabel('error')
print('ENKFPO_error_2','-djpeg')

%%

toc()