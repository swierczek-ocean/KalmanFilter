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
Q = 1;
R = 1;
nsteps = 7500;
Traj = zeros(3,nsteps);
PF = zeros(3,nsteps);
PF(:,1) = randn(3,1);
Ne = 20;
w = zeros(1,Ne);
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
W = ones(Ne);

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot3(PF(3,1),PF(1,1),PF(2,1),'*','Color',Color(:,9),'MarkerSize',2)
axis([-10 60 -25 25 -30 30])
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_3d_plus_SPF_whole.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:nsteps
    Ens = lorenz63s(Ens,dt,sigma,beta,rho,Q);
    for kk=1:Ne
        w(kk) = 0.5.*(Obs(:,jj)-H*Ens(:,kk))'*(eye(2)./R)*(Obs(:,jj)-H*Ens(:,kk));
    end
    W = normalizeweights(w);
    Ens = resamplingmmo(W,Ens,Ne,3);
    PF(:,jj) = mean(Ens,2);
    if(mod(jj,10)==0)
        plot3(Traj(3,1:jj),Traj(1,1:jj),Traj(2,1:jj),'Color',Color(:,11),'Linewidth',1.3)
        hold on
        plot3(PF(3,1:jj),PF(1,1:jj),PF(2,1:jj),'*','Color',Color(:,9),'MarkerSize',2)
        axis([-10 60 -25 25 -30 30])
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

%%

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,22),'MarkerSize',3)
xlabel('time step')
ylabel('error')



%% ensemble + simulation + movie
Ens = randn(3,Ne);
W = ones(Ne);

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot(PF(1,1),PF(3,1),'*','Color',Color(:,9),'MarkerSize',2)
axis([-25 25 -10 60])
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_2d_plus_SPF_whole.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:nsteps
    Ens = lorenz63s(Ens,dt,sigma,beta,rho,Q);
    for kk=1:Ne
        w(kk) = 0.5.*(Obs(:,jj)-H*Ens(:,kk))'*(eye(2)./R)*(Obs(:,jj)-H*Ens(:,kk));
    end
    W = normalizeweights(w);
    Ens = resamplingmmo(W,Ens,Ne,3);
    PF(:,jj) = mean(Ens,2);
    if(mod(jj,10)==0)
        plot(Traj(1,1:jj),Traj(3,1:jj),'Color',Color(:,11),'Linewidth',1.3)
        hold on
        plot(PF(1,1:jj),PF(3,1:jj),'*','Color',Color(:,9),'MarkerSize',2)
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
plot3(PF(3,:),PF(1,:),PF(2,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
hold off

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,11),'Linewidth',1.3)
hold on
plot(PF(1,:),PF(3,:),'*','Color',Color(:,9),'MarkerSize',2)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
hold off

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,22),'MarkerSize',3)
xlabel('time step')
ylabel('error')


%%

toc()