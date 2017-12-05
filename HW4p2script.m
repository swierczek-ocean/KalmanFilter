colors
bg=3000;
n=40;
dt=0.01;
obsdt=0.2;
ne=20;
jump = ceil(obsdt/dt);
R=1;
t_final=100;
r=4;
alpha=0.1;
spy=16;
spl=100;
Rm = R*eye(ne);
F = 8;
[L1,L2,H] = prelim(n);

[SynthDataTrue,SynthDataObs,X_start] = lorenz3(n,t_final,L1,L2,H,F,dt,jump,R);

T = SynthDataTrue(:,1:spl*jump);
Y = SynthDataObs(:,1:spl);

T2 = SynthDataTrue(:,spl*jump+1:end);
Y2 = SynthDataObs(:,spl+1:end);

j1 = size(T2,2);
j2 = size(Y2,2);
RMSE = zeros(1,j2-1);
spread = zeros(1,j2-1);

ensemble = ensemble_init4(X_start,L1,L2,F,dt,ne,n);
[ARMSE,aspread,X_a,mu_a,P_a] = enkfpo4(ensemble,Y,T,r,alpha,spy,L1,L2,H,R,dt,jump,n,ne);
ARMSE
aspread

X=X_a;

bmean=mu_a;
x0=bmean;
bcov=P_a;
spyvec = zeros(1,j2-1);
indices = zeros(1,j2-1);
L = localize2(n,r);

for i=1:j2-1
    y_t = Y2(:,i+1);
    [xT,x0] = fdvar(x0,dt,jump,L1,L2,H,bcov,bmean,y_t,n,R,F);
    [X,bmean,bcov] = ENKFPO(X,dt,jump,n,ne,H,R,L1,L2,F,r,alpha,y_t);
    x0=xT;
    X = X-bmean+xT;
    error = xT-T2(:,jump*i+1);
    RMSE(i) = sqrt((1/n).*transpose(error)*error);
    spread(i) = sqrt(trace(bcov)/n);
    spyvec(i) = xT(spy);
    indices(i) = i*jump+1;
end

figure
h1=plot(1:j2-1,RMSE,'*','MarkerSize',7,'Color',Color(:,9));
hold on
h2=plot(1:j2-1,spread,'o','MarkerSize',7,'Color',Color(:,8));
title('4DVar Errors')
xlabel('time')
legend([h1(1),h2(1)],'root mean square error','spread')
print(['ErrorsPO_r=',num2str(r),'_alpha=',num2str(alpha)],'-djpeg')
hold off

figure
h1=plot(1:j2-1,spyvec,'*','MarkerSize',7,'Color',Color(:,7));
hold on
h2=plot(1:j2-1,T2(spy,indices),'o','MarkerSize',7,'Color',Color(:,14));
title(['True Data vs. 4DVar Estimate in coordinate ', num2str(spy)])
xlabel('time')
legend([h1(1),h2(1)],'Kalman','True')
print('KalmanvsTrue','-djpeg')
hold off


average_RMSE = mean(RMSE(30:end))

