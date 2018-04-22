colors

% t = 0:0.01:25;
% b = 0.25;
% k = 4;
% B = b/(2*sqrt(k-b^2/4));
% X = exp(-b/2.*t).*(cos(sqrt(k-b^2/4).*t) + B.*sin(sqrt(k-b^2/4).*t));
% 
% plot(t,X,'LineWidth',4)
% xlabel('time')
% ylabel('X(t)')
% title('m=1  b=0.25  k=4  x(0)=1  dx/dt(0)=0')


t = 0:0.5:500;
X = 15.*sin(pi.*t./12-3)+75;
Y = -10.*sin(pi.*t./12-3.5)+20;

% yyaxis left
% h1=plot(t,X,'LineWidth',2.5,'Color',Color(:,11));
% axis([-2 52 55 95])
% xlabel('time (hours)')
% ylabel('degrees Fahrenheit')
% 
% yyaxis right
% h2=plot(t,Y,'LineWidth',2.5,'Color',Color(:,16))
% axis([-2 52 5 35])
% ylabel('%')
% legend([h1(1),h2(1)],'temperature','relative humidity','Location','northwest')
% title('temp & humidity at a single point')


% plot3(Y,t,X,'LineWidth',3,'Color',Color(:,11))
% axis([5 35 -1 51 55 95])
% ylabel('time (hours)')
% zlabel('degrees Fahrenheit')
% xlabel('% relative humidity')
% title('temp & humidity at a single point')


%% 3D movie 2

fig = figure, set(gcf, 'Color','white')
plot(X(1),Y(1),'Color',Color(:,11),'Linewidth',1.3)
axis([55 95 5 35])
xlabel('degrees Fahrenheit')
ylabel('% relative humidity')
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('temphumid.avi');
vidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
writeVideo(vidObj, getframe(fig));

nsteps = length(t);

for jj=2:nsteps
    plot(X(jj-1:jj),Y(jj-1:jj),'Color',Color(:,11),'Linewidth',1.3)
    hold on
    axis([55 95 5 35])
    xlabel('degrees Fahrenheit')
    ylabel('% relative humidity')
    hold off
    drawnow()
    writeVideo(vidObj, getframe(fig));
end

close(vidObj);
%%





