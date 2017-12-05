function [ARMSE,aspread,X_a,mu_a,P_a] = enkfpo5(ensemble,Y,T,r,alpha,spy,L1,L2,H,R,dt,jump,n,ne)

colors
f = @(x)(L1*x.*(L2*x)+(8-x));
[m,a]=size(Y);
X = ensemble;
RMSE = zeros(a,1);
spread = zeros(a,1);
Rm = R*eye(m);
time = [1:1:a];
spyvec = zeros(1,a-1);
L = localize2(n,r);

for i=1:a
    mu_f = mean(X,2);                                   % forecast mean
    X = mu_f + sqrt(1+alpha).*(X-mu_f);                 % ensemble inflation
    P_f = cov(X');                                      % forecast covariance
    P_f = L.*P_f;                                       % localization
    K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
    Y_tilde = Y(:,i)+normrnd(0,R,m,ne);                 % analysis perturbations
    X_a = X + K*(Y_tilde-H*X);                          % analysis ensemble
    mu_a = mu_f + K*(Y(:,i)-H*mu_f);                    % analysis mean
    P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
    error = mu_a-T(:,jump*(i-1)+1);
    RMSE(i) = sqrt((1/n).*transpose(error)*error);
    spread(i) = sqrt(trace(P_a)/n);
    spyvec(i) = mu_a(spy);
    X = X_a;
    for j=1:jump
        k1 = f(X);
        k2 = f(X+0.5*dt.*k1);
        k3 = f(X+2*dt.*k2 -dt.*k1);
        X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);      % forecast ensemble
    end
    
end

ARMSE = mean(RMSE(20:end));
aspread = mean(spread(20:end));

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
h2=plot(time,T(spy,time),'o','MarkerSize',7,'Color',Color(:,14));
title(['True Data vs. Kalman Estimate in coordinate ', num2str(spy)])
xlabel('time')
legend([h1(1),h2(1)],'Kalman','True')
print('KalmanvsTrue','-djpeg')
hold off

end
