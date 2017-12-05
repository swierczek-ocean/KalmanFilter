function [ARMSE,aspread,X_a,mu_a,P_a] = enkfpo4(ensemble,Y,T,r,alpha,spy,L1,L2,H,R,dt,jump,n,ne)

colors
f = @(x)(L1*x.*(L2*x)+(8-x));
[m,q]=size(T);
[m,a]=size(Y);
X = ensemble;
RMSE = [];
spread = [];
Rm = R*eye(m);
time = 1:a;
counter = 0;
spyvec = zeros(1,a-1);
indices = [];
L = localize2(n,r);

for i=1:q
%     k1 = f(X);
%     k2 = f(X+0.5*dt.*k1);
%     k3 = f(X+2*dt.*k2 -dt.*k1);
%     X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);      % forecast ensemble
    
if(mod(i,jump)==1)
    counter = counter + 1;
    mu_f = mean(X,2);                                   % forecast mean
    P_f = cov(X');                                      % forecast covariance
    P_f = (1+alpha).*L.*P_f;                            % localization
    K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
    Y_tilde = Y(:,counter)+normrnd(0,R,m,ne);           % analysis perturbations
    X_a = X + K*(Y_tilde-H*X);                          % analysis ensemble
    mu_a = mu_f + K*(Y(:,counter)-H*mu_f);              % analysis mean
    P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
    error = mu_a-T(:,jump*(counter-1)+1);
    RMSE = [RMSE,sqrt((1/n).*transpose(error)*error)];
    spread = [spread,sqrt(trace(P_a)/n)];
    spyvec(counter) = mu_a(spy);
    indices = [indices,i];
    X = X_a;
    
end
k1 = f(X);
k2 = f(X+0.5*dt.*k1);
k3 = f(X+2*dt.*k2 -dt.*k1);
X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);      % forecast ensemble
    
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

