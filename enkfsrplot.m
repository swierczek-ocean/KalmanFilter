function [RMSE,spread] = enkfsrplot(dt,ensemble,M,N,H,t_final,R,Y,T,jump,X_start,j)
tic();

f = @(x)(M*x.*(N*x)+(8-x));
[n,ne]=size(ensemble);
numiter = ceil(t_final/dt);
[m,q]=size(Y);
X = ensemble;
RMSE = [];
spread = [];
Rm = R*eye(m);
time = [1:1:numiter-1];
xmin = min(ensemble(j,:))-1;
xmax = max(ensemble(j,:))+1;
ymin = min(ensemble(j+1,:))-1;
ymax = max(ensemble(j+1,:))+1;
zmin = min(ensemble(j+2,:))-1;
zmax = max(ensemble(j+2,:))+1;

figure, set(gcf, 'Color','white')
plot3(X_start(j),X_start(j+1),X_start(j+2))
axis([xmin xmax ymin ymax zmin zmax])
set(gca, 'nextplot','replacechildren', 'Visible','off');

nFrames = 471;
vidObj = VideoWriter('lorenz.avi');
vidObj.Quality = 100;
vidObj.FrameRate = FR;
open(vidObj);
writeVideo(vidObj, getframe(gca));

for i=1:(numiter-1)
    x_f = f(X);
    mu_f = (1/ne).*transpose(sum(transpose(x_f)));
    X_f = (x_f - mu_f).*(1/sqrt(ne-1));
    P_f = X_f*transpose(X_f);
    K = P_f*H'/(H*P_f*H' + Rm);
    Y_tilde = Y(:,i)+normrnd(0,R,m,ne);
    X = x_f + K*(Y_tilde-H*x_f);
    mu_a = mu_f + K*(Y(:,i)-H*mu_f);
    P_a = (eye(n)-K*H)*P_f;
    error = mu_a-T(:,i+1);
    RMSE = [RMSE,sqrt((1/n).*transpose(error)*error)];
    spread = [spread,sqrt(trace(P_a)/n)];
    
    if mod(i,jump)==0
        plot3(X_start(j),X_start(j+1),X_start(j+2))
        axis([xmin xmax ymin ymax zmin zmax])
        drawnow()
        writeVideo(vidObj, getframe(gca));
    end
end

close(vidObj);

figure
plot(time,RMSE,'*','MarkerSize',5,'Color','red')
hold on
plot(time,spread,'o','MarkerSize',5,'Color','blue')
title('Errors')
xlabel('time')
legend('root mean square error','spread')
print('Errors','-djpeg')

toc()
end

