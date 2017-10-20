dt = 0.01;
t_final = 10;
R = 1;
obsdt = 0.1;

[M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(40,dt,t_final,8,R,obsdt);


T = SynthDataTrue;
Y = SynthDataObs;
size(Y);
size(T);

ne=400;
[ensemble] = ensemble_init(dt,ne,M,N,8,X_start);
W = randomrotation(ne);

%[ARMSE,aspread] = enkfpo(obsdt,ensemble,M,N,H,t_final,R,Y,T,jump)

[ARMSE,aspread] = enkfsr2(dt,ensemble,M,N,H,t_final,R,Y,T,jump,W)

%[RMSE,spread] = enkfpoplot(dt,ensemble,M,N,H,t_final,R,Y,T,20,X_start,1);

%[M,N,H,SynthDataTrue,SynthDataObs,X_start] = lorenzplot(40,dt,t_final,8,R,15,1,30);