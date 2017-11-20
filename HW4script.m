global n, global dt, global obsdt, global ne, global jump, global R, global bg
bg=3000;
n=40;
dt=0.01;
obsdt=0.2;
ne=20;
jump = ceil(obsdt/dt);
R=1;
t_final=25;
SampleSize=100;

run(colors)

global L1, global L2, global H, global F
F = 8;
[L1,L2,H] = prelim(n);

tic()
[SynthDataTrue,SynthDataObs,X_start] = lorenz2(n,t_final);
toc()
tic()
ensemble = ensemble_init3(X_start);
toc()
size(ensemble)

global bmean
global bcov
bmean = transpose(mean(transpose(ensemble)));
bcov = (ensemble-bmean)*transpose(ensemble-bmean)./(bg-1);

x_0 = bmean;
global y_t
y_t = SynthDataObs(:,2);
[x_start,x_star,TimeSeries] = fdvar(x_0);

[arrr,Jaco] = argh(x_star);
M = Mlin(TimeSeries);
Bcov_star = (2.*Jaco'*Jaco)\eye(n);

px0=zeros(n,SampleSize);
px0y=zeros(n,SampleSize);
pxT=zeros(n,SampleSize);
pxTy=zeros(n,SampleSize);

Bcov_star = 0.5.*(Bcov_star + Bcov_star');
P = (M*bcov*M');
P = 0.5.*(P+P');

for i=1:SampleSize
    px0(:,i)=bmean+sqrtm(bcov)*randn(40,1);
    px0y(:,i)=x_star + sqrtm(Bcov_star)*randn(40,1);
    pxT(:,i)=MRK2(bmean) + sqrtm(P)*randn(40,1);
    pxTy(:,i)=x_start+sqrtm(M*Bcov_star\M')*randn(40,1);
end


figure
h1=plot(1:40,px0','Color',Color(:,19),'Linewidth',1.5);
hold on
h2=plot(1:40,px0y','Color',Color(:,15),'Linewidth',1.5);
h3=plot(1:40,SynthDataTrue(:,1),'Color',Color(:,9),'Linewidth',3.5);
title('True vs. Prior vs. Posterior')
xlabel('dimension')
legend([h1(1),h2(1),h3(1)],'prior','posterior','true')
print('prior_vs_posterior at time 0','-djpeg')


figure
h1=plot(1:40,pxT','Color',Color(:,19),'Linewidth',2);
hold on
h2=plot(1:40,pxTy','Linewidth',2,'Color',Color(:,15));
h3=plot(1:40,SynthDataTrue(:,21),'Color',Color(:,9),'Linewidth',3.5);
h4=plot(1:2:39,SynthDataObs(:,2),'*','MarkerSize',10,'Color','y');
title('True vs. M(Prior) vs. M(Posterior) vs. Obs at time t')
xlabel('dimension')
legend([h1(1),h2(1),h3(1),h4(1)],'prior','posterior','true','obs')
print('prior_vs_posterior_vs_obs at time t','-djpeg')


RMSE_Background0=sqrt(1/n*sum((SynthDataTrue(:,1)-bmean).^2))
RMSE_x0=sqrt(1/n*sum((SynthDataTrue(:,1)-x_star).^2))
RMSE_BackgroundT=sqrt(1/n*sum((SynthDataTrue(:,21)-MRK2(bmean)).^2))
RMSE_xT=sqrt(1/n*sum((SynthDataTrue(:,21)-x_start).^2))
