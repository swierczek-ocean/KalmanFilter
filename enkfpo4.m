function [ARMSE,aspread,X_a,mu_a,P_a] = enkfpo4(ensemble,Y,T,r,alpha,spy)
global L1, global L2, global H, global R, global dt, global jump
global n, global ne, 

colors
f = @(x)(L1*x.*(L2*x)+(8-x));
[m,q]=size(T);
[m,a]=size(Y);
X = ensemble;
RMSE = [];
spread = [];
Rm = R*eye(m);
time = [1:1:a-1];
counter = 0;
spyvec = zeros(1,a-1);
indices = [];

for i=1:(q-1)
    k = f(X) + f(X+dt.*f(X)); 
    X = X + 0.5*dt.*k;              % forecast ensemble
    
    if(mod(i,jump)==0)
        counter = counter + 1;
        mu_f = (1/ne).*transpose(sum(transpose(X)));        % forecast mean
        x_f = mu_f + sqrt(1+alpha).*(X-mu_f);               % ensemble inflation
        X_f = (x_f - mu_f).*(1/sqrt(ne-1));                 % forecast perturbations
        P_f = X_f*transpose(X_f);                           % forecast covariance
        L = localize2(P_f,r);                               % creating localization matrix L
        P_f = L.*P_f;                                       % localization
        K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
        Y_tilde = Y(:,counter+1)+normrnd(0,R,m,ne);         % analysis perturbations
        X_a = X + K*(Y_tilde-H*X);                          % analysis ensemble
        mu_a = mu_f + K*(Y(:,counter+1)-H*mu_f);            % analysis mean
        P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
        error = mu_a-T(:,jump*counter+1);
        RMSE = [RMSE,sqrt((1/n).*transpose(error)*error)];
        spread = [spread,sqrt(trace(P_a)/n)];
        spyvec(counter) = mu_a(spy);
        indices = [indices,i];
        X = X_a;
    end
end

[m,q] = size(RMSE);
z = floor(q/2);
ARMSE = (1/(q-1-z))*sum(transpose(RMSE(z:q-1)));
aspread = (1/(q-1-z))*sum(transpose(spread(z:q-1)));

figure
h1=plot(time,RMSE,'*','MarkerSize',7,'Color',Color(:,9));
hold on
h2=plot(time,spread,'o','MarkerSize',7,'Color',Color(:,8));
title('EnKF Perturbed Obs Errors')
xlabel('time')
legend([h1(1),h2(1)],'root mean square error','spread')
print(['ErrorsPO_r=',num2str(r),'_alpha=',num2str(alpha)],'-djpeg')
hold off

figure
h1=plot(time,spyvec,'*','MarkerSize',7,'Color',Color(:,7));
hold on
h2=plot(time,T(spy,indices),'o','MarkerSize',7,'Color',Color(:,14));
title(['True Data vs. Kalman Estimate in coordinate ', num2str(spy)])
xlabel('time')
legend([h1(1),h2(1)],'Kalman','True')
print('KalmanvsTrue','-djpeg')
hold off

end

