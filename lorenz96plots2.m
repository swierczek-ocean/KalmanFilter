function lorenz96plots2(Traj,color1,coords1,name,nsteps,xx,yy,zz)
colors

%% 3D movie 2

figure, set(gcf, 'Color','white')
plot3(Traj(xx,1),Traj(yy,1),Traj(zz,1),'Color',Color(:,color1),'Linewidth',1.3)
axis(coords1)
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz96_3d_whole.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

nsteps = floor(0.25*nsteps);

for jj=2:nsteps
    if(mod(jj,10)==0)
        plot3(Traj(xx,1:jj),Traj(yy,1:jj),Traj(zz,1:jj),'Color',Color(:,color1),'Linewidth',1.3)
        axis(coords1)
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%


end



