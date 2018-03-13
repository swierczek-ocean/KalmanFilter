function lorenz96plots(Traj,PF,color1,color2,color3,coords1,name,ll,nsteps,xx,yy,zz)
colors

%% plots of trajectory in 2D and 3D

figure
plot3(Traj(xx,:),Traj(yy,:),Traj(zz,:),'Color',Color(:,color1),'Linewidth',1.5)
axis(coords1)
xlabel('x')
ylabel('y')
zlabel('z')
print([name,'_L96_3D_trajectory'],'-djpeg')

%%

%% 3D movie 1

figure, set(gcf, 'Color','white')
plot3(Traj(xx,1),Traj(yy,1),Traj(zz,1),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(PF(xx,1),PF(yy,1),PF(zz,1),'*','Color',Color(:,color2),'MarkerSize',2)
axis(coords1)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz96_3d.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:(ll+1)
    if(mod(jj,10)==0)
        plot3(Traj(xx,1:jj),Traj(yy,1:jj),Traj(zz,1:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(PF(xx,1:jj),PF(yy,1:jj),PF(zz,1:jj),'*','Color',Color(:,color2),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for jj=(ll+2):nsteps
    if(mod(jj,10)==0)
        plot3(Traj(xx,jj-ll:jj),Traj(yy,jj-ll:jj),Traj(zz,jj-ll:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(PF(xx,jj-ll:jj),PF(yy,jj-ll:jj),PF(zz,jj-ll:jj),'*','Color',Color(:,color2),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%


%% Error plot

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'*','Color',Color(:,color3),'MarkerSize',3)
xlabel('time step')
ylabel('error')
print([name,'_L96_error_1'],'-djpeg')
%%

%% plots of assimilation
figure
plot3(Traj(xx,:),Traj(yy,:),Traj(zz,:),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(PF(xx,:),PF(yy,:),PF(zz,:),'*','Color',Color(:,color2),'MarkerSize',2)
axis(coords1)
xlabel('x')
ylabel('y')
zlabel('z')
print([name,'_L963D_trajectory'],'-djpeg')
hold off

%% 

%% 3D movie 2

figure, set(gcf, 'Color','white')
plot3(Traj(xx,1),Traj(yy,1),Traj(zz,1),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(PF(xx,1),PF(yy,1),PF(zz,1),'*','Color',Color(:,color2),'MarkerSize',2)
axis(coords1)
hold off
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
        hold on
        plot3(PF(xx,1:jj),PF(yy,1:jj),PF(zz,1:jj),'*','Color',Color(:,color2),'MarkerSize',2)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%


end

