function lorenz63plots2(Traj,PF,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
colors

%% Error plot

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,color3),'MarkerSize',3)
xlabel('time step')
ylabel('error')
print([name,'_error_1'],'-djpeg')
%%

%% plots of trajectory in 2D and 3D

figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,color1),'Linewidth',1.5)
axis(coords1)
xlabel('z')
ylabel('x')
zlabel('y')
print([name,'_3D_trajectory'],'-djpeg')

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,color4),'Linewidth',1.5)
axis(coords2)
xlabel('x')
ylabel('z')
print([name,'_2D_trajectory'],'-djpeg')
%%

%% 3D movie 1

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(PF(3,1),PF(1,1),PF(2,1),'*','Color',Color(:,color2),'MarkerSize',2)
axis(coords1)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz63_3d.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:(ll+1)
    if(mod(jj,10)==0)
        plot3(Traj(3,1:jj),Traj(1,1:jj),Traj(2,1:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(PF(3,1:jj),PF(1,1:jj),PF(2,1:jj),'*','Color',Color(:,color2),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for jj=(ll+2):nsteps
    if(mod(jj,10)==0)
        plot3(Traj(3,jj-ll:jj),Traj(1,jj-ll:jj),Traj(2,jj-ll:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(PF(3,jj-ll:jj),PF(1,jj-ll:jj),PF(2,jj-ll:jj),'*','Color',Color(:,color2),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%

%% 2D movie 1

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,color4),'Linewidth',1.3)
hold on
plot(PF(1,1),PF(3,1),'*','Color',Color(:,color5),'MarkerSize',2)
axis(coords2)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz63_2d.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:(ll+1)
    if(mod(jj,10)==0)
        plot(Traj(1,1:jj),Traj(3,1:jj),'Color',Color(:,color4),'Linewidth',1.6)
        hold on
        plot(PF(1,1:jj),PF(3,1:jj),'*','Color',Color(:,color5),'MarkerSize',3)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for jj=(ll+2):nsteps
    if(mod(jj,10)==0)
        plot(Traj(1,jj-ll:jj),Traj(3,jj-ll:jj),'Color',Color(:,color4),'Linewidth',1.6)
        hold on
        plot(PF(1,jj-ll:jj),PF(3,jj-ll:jj),'*','Color',Color(:,color5),'MarkerSize',3)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

%%

%% plots of assimilation
figure
plot3(Traj(3,:),Traj(1,:),Traj(2,:),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(PF(3,:),PF(1,:),PF(2,:),'*','Color',Color(:,color2),'MarkerSize',2)
axis(coords1)
xlabel('z')
ylabel('x')
zlabel('y')
print([name,'3D_trajectory'],'-djpeg')
hold off

figure
plot(Traj(1,:),Traj(3,:),'Color',Color(:,color4),'Linewidth',1.3)
hold on
plot(PF(1,:),PF(3,:),'*','Color',Color(:,color5),'MarkerSize',2)
axis(coords2)
xlabel('x')
ylabel('z')
print([name,'2D_trajectory'],'-djpeg')
hold off
%% 

%% 3D movie 2

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(PF(3,1),PF(1,1),PF(2,1),'*','Color',Color(:,color2),'MarkerSize',2)
axis(coords1)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz63_3d_whole.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

nsteps = floor(0.25*nsteps);

for jj=2:nsteps
    if(mod(jj,10)==0)
        plot3(Traj(3,1:jj),Traj(1,1:jj),Traj(2,1:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(PF(3,1:jj),PF(1,1:jj),PF(2,1:jj),'*','Color',Color(:,color2),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%

%% 2D movie 2

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,color4),'Linewidth',1.3)
hold on
plot(PF(1,1),PF(3,1),'*','Color',Color(:,color5),'MarkerSize',2)
axis(coords2)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz63_2d_whole.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:nsteps
    if(mod(jj,10)==0)
        plot(Traj(1,1:jj),Traj(3,1:jj),'Color',Color(:,color4),'Linewidth',1.6)
        hold on
        plot(PF(1,1:jj),PF(3,1:jj),'*','Color',Color(:,color5),'MarkerSize',3)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%


end

