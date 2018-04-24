function lorenz63plots5(Traj,PF,Obs,color1,color2,color3,color4,color5,coords1,coords2,name,ll,nsteps)
colors

%% Error plot

Error = Traj - PF;
Error = sqrt(sum(Error.^2))./sqrt(3);
st = floor(0.25*nsteps);
average_RMSE = mean(Error(st:end))

figure
plot(Error,'.','Color',Color(:,color3),'MarkerSize',6)
title('RMSE with DA')
xlabel('time step')
ylabel('error')
print([name,'_error_1'],'-djpeg')
%%


[n,m] = size(PF);
for iii=1:floor(m/2)
   PF(:,2*iii-1)=PF(:,2*iii); 
end

%% 3D movie 1

% figure, set(gcf, 'Color','white')
% plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,color1),'Linewidth',1.3)
% hold on
% plot3(PF(3,1),PF(1,1),PF(2,1),'*','Color',Color(:,color2),'MarkerSize',6)
% axis(coords1)
% hold off
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% nFrames = 471;
% vidObj = VideoWriter([name,'_lorenz63_3d.avi']);
% vidObj.Quality = 100;
% vidObj.FrameRate = 30;
% open(vidObj);
% writeVideo(vidObj, getframe(gca));
% 
% for jj=2:(ll+1)
%     if(mod(jj,10)==0)
%         plot3(Traj(3,1:jj),Traj(1,1:jj),Traj(2,1:jj),'Color',Color(:,color1),'Linewidth',1.6)
%         hold on
%         plot3(PF(3,1:jj),PF(1,1:jj),PF(2,1:jj),'.','Color',Color(:,color2),'MarkerSize',6)
%         axis(coords1)
%         hold off
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
% 
% for jj=(ll+2):nsteps
%     if(mod(jj,10)==0)
%         plot3(Traj(3,jj-ll:jj),Traj(1,jj-ll:jj),Traj(2,jj-ll:jj),'Color',Color(:,color1),'Linewidth',1.6)
%         hold on
%         plot3(PF(3,jj-ll:jj),PF(1,jj-ll:jj),PF(2,jj-ll:jj),'.','Color',Color(:,color2),'MarkerSize',6)
%         axis(coords1)
%         hold off
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
% 
% close(vidObj);
%%

%% 2D movie 1
xx=20;
yy=10;
zz=1;
obssz = 13;
kfsz = 11;

figure, set(gcf, 'Color','white')
plot(Traj(1,1),Traj(3,1),'Color',Color(:,color4),'Linewidth',1.3)
hold on
plot(PF(1,1),PF(3,1),'.','Color',Color(:,color5),'MarkerSize',6)
axis(coords2)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter([name,'_lorenz63_2d.avi']);
vidObj.Quality = 100;
vidObj.FrameRate = 26;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:(ll+1)
    if(mod(jj,10)==0)
        zz = yy;
        yy = xx;
        xx = jj;
        plot(Traj(1,1:jj+1),Traj(3,1:jj+1),'Color',Color(:,color4),'Linewidth',3.5)
        hold on
        plot(PF(1,1:2:jj),PF(3,1:2:jj),'.','Color',Color(:,color5),'MarkerSize',kfsz)
        plot(Obs(1,jj),Obs(2,jj),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,yy),Obs(2,yy),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,zz),Obs(2,zz),'.','Color',Color(:,color3),'MarkerSize',obssz)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    else
        plot(Traj(1,1:jj+1),Traj(3,1:jj+1),'Color',Color(:,color4),'Linewidth',3.5)
        hold on
        plot(PF(1,1:2:jj),PF(3,1:2:jj),'.','Color',Color(:,color5),'MarkerSize',kfsz)
        plot(Obs(1,xx),Obs(2,xx),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,yy),Obs(2,yy),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,zz),Obs(2,zz),'.','Color',Color(:,color3),'MarkerSize',obssz)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));        
    end
end

for jj=(ll+2):(nsteps-2)
    if(mod(jj,10)==0)
        zz = yy;
        yy = xx;
        xx = jj;
        plot(Traj(1,jj-ll:jj+1),Traj(3,jj-ll:jj+1),'Color',Color(:,color4),'Linewidth',3.5)
        hold on
        plot(PF(1,jj-ll:2:jj),PF(3,jj-ll:2:jj),'.','Color',Color(:,color5),'MarkerSize',kfsz)
        plot(Obs(1,jj),Obs(2,jj),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,yy),Obs(2,yy),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,zz),Obs(2,zz),'.','Color',Color(:,color3),'MarkerSize',obssz)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    else
        plot(Traj(1,jj-ll:jj+1),Traj(3,jj-ll:jj+1),'Color',Color(:,color4),'Linewidth',3.5)
        hold on
        plot(PF(1,jj-ll:2:jj),PF(3,jj-ll:2:jj),'.','Color',Color(:,color5),'MarkerSize',kfsz)
        plot(Obs(1,xx),Obs(2,xx),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,yy),Obs(2,yy),'.','Color',Color(:,color3),'MarkerSize',obssz)
        plot(Obs(1,zz),Obs(2,zz),'.','Color',Color(:,color3),'MarkerSize',obssz)
        axis(coords2)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));        
        
    end
end

close(vidObj);

%%



%% 3D movie 2
% 
% figure, set(gcf, 'Color','white')
% plot3(Traj(3,1),Traj(1,1),Traj(2,1),'Color',Color(:,color1),'Linewidth',1.3)
% hold on
% plot3(PF(3,1),PF(1,1),PF(2,1),'.','Color',Color(:,color2),'MarkerSize',5)
% axis(coords1)
% hold off
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% nFrames = 471;
% vidObj = VideoWriter([name,'_lorenz63_3d_whole.avi']);
% vidObj.Quality = 100;
% vidObj.FrameRate = 30;
% open(vidObj);
% writeVideo(vidObj, getframe(gca));
% 
% nsteps = floor(0.25*nsteps);
% 
% for jj=2:nsteps
%     if(mod(jj,10)==0)
%         plot3(Traj(3,1:jj),Traj(1,1:jj),Traj(2,1:jj),'Color',Color(:,color1),'Linewidth',1.6)
%         hold on
%         plot3(PF(3,1:jj),PF(1,1:jj),PF(2,1:jj),'.','Color',Color(:,color2),'MarkerSize',5)
%         axis(coords1)
%         hold off
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
% 
% close(vidObj);
% %%
% 
% %% 2D movie 2
% 
% figure, set(gcf, 'Color','white')
% plot(Traj(1,1),Traj(3,1),'Color',Color(:,color4),'Linewidth',1.3)
% hold on
% plot(PF(1,1),PF(3,1),'.','Color',Color(:,color5),'MarkerSize',5)
% axis(coords2)
% hold off
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% nFrames = 471;
% vidObj = VideoWriter([name,'_lorenz63_2d_whole.avi']);
% vidObj.Quality = 100;
% vidObj.FrameRate = 12;
% open(vidObj);
% writeVideo(vidObj, getframe(gca));
% 
% for jj=2:nsteps
%     if(mod(jj,2)==0)
%         plot(Traj(1,1:jj),Traj(3,1:jj),'Color',Color(:,color4),'Linewidth',1.6)
%         hold on
%         plot(PF(1,1:jj),PF(3,1:jj),'.','Color',Color(:,color5),'MarkerSize',5)
%         axis(coords2)
%         hold off
%         drawnow()
%         writeVideo(vidObj, getframe(gca));
%     end
% end
% 
% close(vidObj);
%%


end



