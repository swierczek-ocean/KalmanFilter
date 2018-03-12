tic()
clc
close all
clear

%% preliminaries
colors
sigma = 10;
beta = 8/3;
rho = 28;
dt = 0.01;
Q = 1;
R = 1;
jump = 10;
nsteps = 20000;
Traj = zeros(3,nsteps);
PF = zeros(3,nsteps);
Ne = 400;
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
print('3D_trajectory_SPF2','-djpeg')

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.5)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
print('2D_trajectory_SPF2','-djpeg')
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
w = zeros(1,Ne);

for jj=2:nsteps
    Ens = lorenz63s4(Ens,dt,M1,M2,M3) + sqrt(dt)*sqrt(Q).*randn(3,Ne);
    if(mod(jj,10)==0)
        for mm=1:Ne
            w(mm) = 0.5.*(Obs(:,jj)-H*Ens(:,mm))'*Z*(Obs(:,jj)-H*Ens(:,mm));
        end
        W = normalizeweights(w);
        Ens = resamplingmmo(W,Ens,Ne,3);
    end
    PF(:,jj) = mean(Ens,2);
end

%%

%% plot
figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot3(PF(3,:),PF(1,:),PF(2,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
print('3D_trajectory_plus_SPF2','-djpeg')
hold off

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot(PF(1,:),PF(3,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
print('2D_trajectory_plus_SPF2','-djpeg')
hold off

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,22),'MarkerSize',3)
xlabel('time step')
ylabel('error')
print('SPF_error_4','-djpeg')

%%

toc()