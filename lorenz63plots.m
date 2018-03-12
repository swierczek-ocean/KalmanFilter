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
Q = 0.0001;
R = 1;
nsteps = 10000;
Traj = zeros(3,nsteps);

%%

%% simulation and plot

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,9),'Linewidth',1.3)
axis([-10 60 -25 25 -30 30])
xlabel('z')
ylabel('x')
zlabel('y')
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_3d.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 60;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for ii=2:nsteps
    Traj(:,ii) = lorenz63s(Traj(:,ii-1),dt,sigma,beta,rho,Q);
    if(mod(ii,8)==0)
        plot3(Traj(3,1:ii),Traj(1,1:ii),Traj(2,1:ii),'Color',Color(:,9),'Linewidth',1.3)
        axis([-10 60 -25 25 -30 30])
        xlabel('z')
        ylabel('x')
        zlabel('y')
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);


figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,8),'Linewidth',1.5)
axis([-25 25 -10 60])
xlabel('x')
ylabel('z')
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_2d.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 60;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for ii=2:nsteps
    Traj(:,ii) = lorenz63s(Traj(:,ii-1),dt,sigma,beta,rho,Q);
    if(mod(ii,8)==0)
        plot(Traj(1,1:ii),Traj(3,1:ii),'Color',Color(:,8),'Linewidth',1.3)
        axis([-25 25 -10 60])
        xlabel('x')
        ylabel('z')
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

%%





















toc()