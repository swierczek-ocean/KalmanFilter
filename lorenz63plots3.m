function lorenz63plots3(Traj,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
colors


%% 3D movie 1

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,color1),'Linewidth',1.3)
axis(coords1)
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
        axis(coords1)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for jj=(ll+2):nsteps
    if(mod(jj,10)==0)
        plot3(Traj(3,jj-ll:jj),Traj(1,jj-ll:jj),Traj(2,jj-ll:jj),'Color',Color(:,color1),'Linewidth',1.6)
        axis(coords1)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%

%% 2D movie 1

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,color4),'Linewidth',1.3)
axis(coords2)
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
        axis(coords2)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for jj=(ll+2):nsteps
    if(mod(jj,10)==0)
        plot(Traj(1,jj-ll:jj),Traj(3,jj-ll:jj),'Color',Color(:,color4),'Linewidth',1.6)
        axis(coords2)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

%%


%% 3D movie 2

figure, set(gcf, 'Color','white')
plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,color1),'Linewidth',1.3)
axis(coords1)
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
        axis(coords1)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%

%% 2D movie 2

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,color4),'Linewidth',1.3)
axis(coords2)
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
        axis(coords2)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%


end

