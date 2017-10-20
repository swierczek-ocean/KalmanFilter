function [ARMSE,aspread] = enkfpo(dt,ensemble,M,N,H,t_final,R,Y,T,jump)
tic();

f = @(x)(M*x.*(N*x)+(8-x));
[n,ne]=size(ensemble);
numiter = ceil(t_final/dt);
[m,q]=size(Y);
X = ensemble;
RMSE = [];
spread = [];
Rm = R*eye(m);
time = [1:1:q-1];

for i=1:(q-1)
    k1=f(X);
    k2=f(X+0.5*dt.*k1);
    k3=f(X+0.5*dt.*k2);
    k4=f(X+dt.*k3);
    x_f = X + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);
    mu_f = (1/ne).*transpose(sum(transpose(x_f)));
    X_f = (x_f - mu_f).*(1/sqrt(ne-1));
    P_f = X_f*transpose(X_f);
    K = P_f*(H')/(H*P_f*H' + Rm);
    Y_tilde = Y(:,i)+normrnd(0,R,m,ne);
    X = x_f + K*(Y_tilde-H*x_f);
    mu_a = mu_f + K*(Y(:,i)-H*mu_f);
    P_a = (eye(n)-K*H)*P_f;
    error = mu_a-T(:,jump*(i-1)+1);
    RMSE = [RMSE,sqrt((1/n).*transpose(error)*error)];
    spread = [spread,sqrt(trace(P_a)/n)];
end

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

toc()
end

