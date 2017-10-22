dt = 0.01;
t_final = 10;
R = 1;
obsdt = 0.1;
threshold = 0.2;
r = 2;
alpha = 0.2;
set(groot, 'DefaultFigureVisible', 'on');

[M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(40,dt,t_final,8,R,obsdt);


T = SynthDataTrue;
Y = SynthDataObs;
size(Y);
size(T);

ne=100;
ensemble = ensemble_init(dt,ne,M,N,8,X_start);
W = randomrotation(ne);


[ARMSE,aspread] = enkfpo3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold,r,alpha,37);

[ARMSE,aspread] = enkfsr3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold,r,alpha,37)

[ARMSE,aspread] = enkfsr4(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold,r,alpha,37)

[ARMSE,aspread] = enkfsr5(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold,r,alpha,37)