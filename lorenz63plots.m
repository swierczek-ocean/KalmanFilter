tic()
% clc
% close all


%% preliminaries
colors
sigma = 10;
beta = 8/3;
rho = 28;
dt = 0.01;
Q = 0.0001;
R = 1;
M1 = [0,0,0;-1,0,0;1,0,0];
M2 = [0,0,0;0,0,1;0,1,0];
M3 = [-sigma,sigma,0;rho,-1,0;0,0,-beta];
nsteps = 10000;
Traj = zeros(3,nsteps);
Traj(:,1) = 2*ones(3,1);

Traj2 = zeros(3,nsteps);
Traj2(:,1) = 2*ones(3,1);
Traj2(1,1) = Traj2(1,1)+0.0000000001;

for ii=2:nsteps
    Traj(:,ii) = lorenz63s4(Traj(:,ii-1),dt,M1,M2,M3);
end

for ii=2:nsteps
    Traj2(:,ii) = lorenz63s4(Traj2(:,ii-1),dt,M1,M2,M3);
end
%%

%% simulation and plot

% figure, set(gcf, 'Color','white')
% plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,9),'Linewidth',1.3)
% hold on
% axis([-10 60 -25 25 -30 30])
% hold off
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% nFrames = 471;
% vidObj = VideoWriter('lorenz63_3d.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 30;
% open(vidObj);
% writeVideo(vidObj, getframe(gca));
% 
% for ii=2:nsteps
%     Traj(:,ii) = lorenz63s4(Traj(:,ii-1),dt,M1,M2,M3);
%     if(mod(ii,8)==0)
%         plot3(Traj(3,1:ii),Traj(1,1:ii),Traj(2,1:ii),'Color',Color(:,9),'Linewidth',1.3)
%         hold on
%         axis([-10 60 -25 25 -30 30])
%         hold off
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
% 
% close(vidObj);


% figure, set(gcf, 'Color','white')
% plot(Traj(1,1),Traj(3,1),'Color',Color(:,8),'Linewidth',1.5)
% axis([-25 25 -10 60])
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% nFrames = 471;
% vidObj = VideoWriter('lorenz63_2d.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 30;
% open(vidObj);
% writeVideo(vidObj, getframe(gca));
% 
% for ii=2:nsteps
%     % Traj(:,ii) = lorenz63s4(Traj(:,ii-1),dt,M1,M2,M3);
%     if(mod(ii,8)==0)
%         plot(Traj(1,1:ii),Traj(3,1:ii),'Color',Color(:,8),'Linewidth',1.3)
%         axis([-25 25 -10 60])
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
% 
% close(vidObj);



figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,31),'Linewidth',3.5)
hold on
plot(Traj2(1,1),Traj2(3,1),'Color',Color(:,33),'Linewidth',1.5)
axis([-25 25 -10 60])
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz63_2d_split.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for ii=2:25
    if(mod(ii,8)==0)
        plot(Traj(1,1:ii),Traj(3,1:ii),'Color',Color(:,31),'Linewidth',3.3)
        hold on
        plot(Traj2(1,1:ii),Traj2(3,1:ii),'Color',Color(:,33),'Linewidth',1.3)
        axis([-25 25 -10 60])
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for ii=26:nsteps
    if(mod(ii,8)==0)
        plot(Traj(1,ii-24:ii),Traj(3,ii-24:ii),'Color',Color(:,31),'Linewidth',3.3)
        hold on
        plot(Traj2(1,ii-24:ii),Traj2(3,ii-24:ii),'Color',Color(:,33),'Linewidth',1.3)
        axis([-25 25 -10 60])
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);


















toc()