%% 3D movie 1
color1 = 15;
color2 = 4;
color3 = 11;
color4 = 11;
ll = 69;
lll = floor(ll/10);
coords1 = [-9 12 -9 12 -9 12];

figure, set(gcf, 'Color','white')
plot3(Traj(xx,1),Traj(yy,1),Traj(zz,1),'Color',Color(:,color1),'Linewidth',1.3)
hold on
plot3(EnKFestSR(xx,1),EnKFestSR(yy,1),EnKFestSR(zz,1),'*','Color',Color(:,color2),'MarkerSize',2)
plot3(EnKFestNR(xx,1),EnKFestNR(yy,1),EnKFestNR(zz,1),'*','Color',Color(:,color3),'MarkerSize',2)
%plot3(EnKFestRR(xx,1),EnKFestRR(yy,1),EnKFestRR(zz,1),'*','Color',Color(:,color4),'MarkerSize',2)
axis(coords1)
hold off
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('ENKF_lorenz96_3d.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 17;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for jj=2:(ll+15)
    if(mod(jj,10)==0)
        jjj = floor(jj/10);
        plot3(Traj(xx,1:jj),Traj(yy,1:jj),Traj(zz,1:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(EnKFestSR(xx,1:jjj),EnKFestSR(yy,1:jjj),EnKFestSR(zz,1:jjj),'*','Color',Color(:,color2),'MarkerSize',3)
        plot3(EnKFestNR(xx,1:jjj),EnKFestNR(yy,1:jjj),EnKFestNR(zz,1:jjj),'*','Color',Color(:,color3),'MarkerSize',3)
        %plot3(EnKFestRR(xx,1:jjj),EnKFestRR(yy,1:jjj),EnKFestRR(zz,1:jjj),'*','Color',Color(:,color4),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

for jj=(ll+16):nsteps
    if(mod(jj,10)==0)
        jjj = floor(jj/10);
        plot3(Traj(xx,jj-ll:jj),Traj(yy,jj-ll:jj),Traj(zz,jj-ll:jj),'Color',Color(:,color1),'Linewidth',1.6)
        hold on
        plot3(EnKFestSR(xx,jjj-lll:jjj),EnKFestSR(yy,jjj-lll:jjj),EnKFestSR(zz,jjj-lll:jjj),'*','Color',Color(:,color2),'MarkerSize',3)
        plot3(EnKFestNR(xx,jjj-lll:jjj),EnKFestNR(yy,jjj-lll:jjj),EnKFestNR(zz,jjj-lll:jjj),'*','Color',Color(:,color3),'MarkerSize',3)
        %plot3(EnKFestRR(xx,jjj-lll:jjj),EnKFestRR(yy,jjj-lll:jjj),EnKFestRR(zz,jjj-lll:jjj),'*','Color',Color(:,color4),'MarkerSize',3)
        axis(coords1)
        hold off
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);
%%