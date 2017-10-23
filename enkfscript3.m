dt = 0.01;
t_final = 10;
R = 1;
obsdt = 0.1;
r = ;
alpha = ;


[M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(400,dt,t_final,8,R,obsdt);


T = SynthDataTrue;
Y = SynthDataObs;
size(Y);
size(T);

ne=20;
ensemble = ensemble_init(dt,ne,M,N,8,X_start);


[ARMSE1,aspread2] = enkfpo3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r,alpha);


[ARMSE1,aspread2] = enkfsr3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r,alpha);


[ARMSE1,aspread2] = enkfsr4(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r,alpha);


[ARMSE1,aspread2] = enkfsr5(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r,alpha);




















