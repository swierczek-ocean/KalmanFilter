global n, global dt, global obsdt, global ne, global jump, global R, global bg
bg=3000;
n=40;
dt=0.01;
obsdt=0.2;
ne=20;
jump = ceil(obsdt/dt);
R=1;
t_final=25;
r=4;
alpha=0.14;
spy=17;
spl=100;
Rm = R*eye(20);

global L1, global L2, global H, global F
F = 8;
[L1,L2,H] = prelim(n);

[SynthDataTrue,SynthDataObs,X_start] = lorenz2(n,t_final);

T = SynthDataTrue(:,1:spl*jump);
Y = SynthDataObs(:,1:spl);

T2 = SynthDataTrue(:,spl*jump+1:end);
Y2 = SynthDataObs(:,spl+1:end);

j1 = size(T2,2);
j2 = size(Y2,2);
RMSE = zeros(1,j2);
spread = zeros(1,j2);

ensemble = ensemble_init4(X_start);
[ARMSE,aspread,X_a,mu_a,P_a] = enkfpo4(ensemble,Y,T,r,alpha,spy);

ne=19;
X_a=X_a(:,1:ne);

global bcov, global bmean, global y_t
bmean=mu_a;
bcov=P_a;
spyvec = zeros(1,j2);
indices = [];

for i=1:j2
    y_t = Y2(:,i);
    [new_mu,special_x,TimeSeries] = fdvar(bmean);
    for l=1:ne
        X_a(:,l)=fdvar(X_a(:,l));
    end
    Xtemp = [new_mu,X_a];
    mu_f = (1/(ne+1)).*transpose(sum(transpose(Xtemp)));% forecast mean
    X_f = (Xtemp - new_mu).*(1/sqrt(ne-1));             % forecast perturbations
    P_f = X_f*transpose(X_f);                           % forecast covariance
    L = localize2(P_f,r);                               % creating localization matrix L
    P_f = (1+alpha).*L.*P_f;                            % localization
    K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
    Y_tilde = y_t+[zeros(20,1),normrnd(0,R,20,ne)];     % analysis perturbations
    X_a = Xtemp + K*(Y_tilde-H*Xtemp);                  % analysis ensemble
    X_a = X_a(:,1:ne);
    bmean = (1/(ne+1)).*transpose(sum(transpose(X_a))); % analysis mean
    bcov = (eye(n)-K*H)*P_f;                            % analysis covariance
    gamma = new_mu-bmean;
    X_a = X_a + gamma;
    bmean = new_mu;
    error = bmean-T2(jump*i+1);
    RMSE(i) = sqrt((1/n).*transpose(error)*error);
    spread(i) = sqrt(trace(bcov)/n);
    spyvec(i) = bmean(spy);
    indices = [indices,i*jump];
end

figure
h1=plot(1:j2,RMSE,'*','MarkerSize',7,'Color',Color(:,9));
hold on
h2=plot(1:j2,spread,'o','MarkerSize',7,'Color',Color(:,8));
title('EnKF Perturbed Obs Errors')
xlabel('time')
legend([h1(1),h2(1)],'root mean square error','spread')
print(['ErrorsPO_r=',num2str(r),'_alpha=',num2str(alpha)],'-djpeg')
hold off

figure
h1=plot(1:j2,spyvec,'*','MarkerSize',7,'Color',Color(:,7));
hold on
h2=plot(1:j2,T(spy,indices),'o','MarkerSize',7,'Color',Color(:,14));
title(['True Data vs. Kalman Estimate in coordinate ', num2str(spy)])
xlabel('time')
legend([h1(1),h2(1)],'Kalman','True')
print('KalmanvsTrue','-djpeg')
hold off




