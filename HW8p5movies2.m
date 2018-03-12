tic()
clc
close all
clear

%% preliminaries
colors
n = 3;
r = 4;
alpha = 0.15;
sigma = 10;
beta = 8/3;
rho = 28;
dt = 0.01;
Q = 1;
R = 1;
jump = 10;
nsteps = 6000;
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

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.5)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
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

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot3(KF(3,1),KF(1,1),KF(2,1),'*','Color',Color(:,9),'MarkerSize',2)
axis([-10 60 -25 25 -30 30])
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_3d_plus_ENKF2_whole.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:nsteps
    Ens = lorenz63s4(Ens,dt,M1,M2,M3) + sqrt(dt)*sqrt(Q).*randn(3,Ne);
    if(mod(jj,10)==0)
        [Ens,mu_a,P_a] = ENKFPO2(Ens,dt,n,Ne,H,R,r,alpha,Obs(:,jj),M1,M2,M3,Q);
    end
    KF(:,jj) = mean(Ens,2);
    if(mod(jj,10)==0)
        plot3(Traj(3,1:jj),Traj(1,1:jj),Traj(2,1:jj),'Color',Color(:,11),'Linewidth',1.6)
        hold on
        plot3(KF(3,1:jj),KF(1,1:jj),KF(2,1:jj),'*','Color',Color(:,9),'MarkerSize',3)
        axis([-10 60 -25 25 -30 30])
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

%%

Error = Traj - KF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,22),'MarkerSize',3)
xlabel('time step')
ylabel('error')

%% ensemble + simulation + movie
Ens = randn(3,Ne);

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot(KF(1,1),KF(3,1),'*','Color',Color(:,9),'MarkerSize',2)
axis([-25 25 -10 60])
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_2d_plus_ENKF2_whole.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:nsteps
    Ens = lorenz63s4(Ens,dt,M1,M2,M3) + sqrt(dt)*sqrt(Q).*randn(3,Ne);
    if(mod(jj,10)==0)
        [Ens,mu_a,P_a] = ENKFPO2(Ens,dt,n,Ne,H,R,r,alpha,Obs(:,jj),M1,M2,M3,Q);
    end
    KF(:,jj) = mean(Ens,2);
    if(mod(jj,10)==0)
        plot(Traj(1,1:jj),Traj(3,1:jj),'Color',Color(:,11),'Linewidth',1.6)
        hold on
        plot(KF(1,1:jj),KF(3,1:jj),'*','Color',Color(:,9),'MarkerSize',3)
        axis([-25 25 -10 60])
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

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
hold off

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot(KF(1,:),KF(3,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
hold off

Error = Traj - KF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,22),'MarkerSize',3)
xlabel('time step')
ylabel('error')

%%

toc()