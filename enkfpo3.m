function [ARMSE,aspread,J] = enkfpo3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold,r,alpha,spy)
tic();

f = @(x)(M*x.*(N*x)+(8-x));
[n,ne]=size(ensemble);
[m,q]=size(T);
[m,a]=size(Y);
X = ensemble;
RMSE = [];
spread = [];
Rm = R*eye(m);
time = [1:1:a-1];
counter = 0;
% spyvec = zeros(1,a-1);
% indices = [];

for i=1:(q-1)
    k1=f(X);
    k2=f(X+0.5*dt.*k1);
    k3=f(X+0.5*dt.*k2);
    k4=f(X+dt.*k3);
    X = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);              % forecast ensemble
    
    if(mod(i,jump)==0)
        counter = counter + 1;
        mu_f = (1/ne).*transpose(sum(transpose(X)));        % forecast mean
        x_f = mu_f + sqrt(1+alpha).*(X-mu_f);               % ensemble inflation
        X_f = (x_f - mu_f).*(1/sqrt(ne-1));                 % forecast perturbations
        P_f = X_f*transpose(X_f);                           % forecast covariance
        J = P_f;
        L = localize(P_f,threshold,r);                      % creating localization matrix L
        P_f = L.*P_f;                                       % localization
        K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
        Y_tilde = Y(:,counter+1)+normrnd(0,R,m,ne);         % analysis perturbations
        X = X + K*(Y_tilde-H*X);                            % analysis ensemble
        mu_a = mu_f + K*(Y(:,counter+1)-H*mu_f);            % analysis mean
        P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
        error = mu_a-T(:,jump*counter+1);
        RMSE = [RMSE,sqrt((1/n).*transpose(error)*error)];
        spread = [spread,sqrt(trace(P_a)/n)];
%         spyvec(counter) = mu_a(spy);
%         indices = [indices,i];
    end
end

[m,q] = size(RMSE);
z = floor(q/2);
ARMSE = (1/(q-1-z))*sum(transpose(RMSE(z:q-1)));
aspread = (1/(q-1-z))*sum(transpose(spread(z:q-1)));

figure
plot(time,RMSE,'*','MarkerSize',5,'Color','red')
hold on
plot(time,spread,'o','MarkerSize',5,'Color','blue')
title('EnKF Perturbed Obs Errors')
xlabel('time')
legend('root mean square error','spread')
print('ErrorsPO','-djpeg')

% figure
% plot(time,spyvec,'*','MarkerSize',5,'Color','red')
% hold on
% plot(time,T(spy,indices),'o','MarkerSize',5,'Color','blue')
% title(['True Data vs. Kalman Estimate in coordinate ', num2str(spy)])
% xlabel('time')
% legend('Kalman','True')
% print('KalmanvsTrue','-djpeg')

% toc();
end

