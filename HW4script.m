global n, global dt, global obsdt, global ne, global jump, global R
n=40;
dt=0.01;
obsdt=0.2;
ne=20;
jump = ceil(obsdt/dt);
R=1;
t_final=10;

global L1, global L2, global H, global F
F = 8;
[L1,L2,H] = prelim(n);

[SynthDataTrue,SynthDataObs,X_start] = lorenz2(n,dt,t_final,R,jump);

ensemble = ensemble_init3(dt,X_start);

bmean = transpose(mean(transpose(ensemble)));
bcov = (ensemble-bmean)*transpose(ensemble-bmean)./(ne-1);

x_0 = mean;
global y_t
y_t = SynthDataObs(:,2);
[x_start,x_star] = fdvar(x_0);

figure
plot(1:40,x_star,'Color',Color(:,16),'Linewidth',5)
hold on
plot(1:40,bmean,'Color',Color(:,9),'Linewidth',5)
plot(1:40,SynthDataTrue(:,1),'Color',Color(:,8),'Linewidth',5)
title('True vs. Prior vs. Posterior')
xlabel('time')
legend('posterior','prior','true')
print('prior_vs_posterior,'-djpeg')


figure
plot(1:40,x_start,'Color',Color(:,16),'Linewidth',5)
hold on
plot(1:40,MRK2(bmean,dt,obsdt,jump),'Color',Color(:,9),'Linewidth',5)
plot(1:40,SynthDataTrue(:,2),'Color',Color(:,8),'Linewidth',5)
title('True vs. Prior vs. Posterior')
xlabel('time')
legend('posterior','prior','true')
print('prior_vs_posterior,'-djpeg')



