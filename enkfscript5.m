dt = 0.01;
t_final = 70;
R = 1;
obsdt = 0.1;
spy=3;


Y = load('Obs40.txt');
Y = transpose(Y);
size(Y);

[M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(40,dt,t_final,8,R,obsdt);

T = SynthDataTrue;
size(T);

ne=20;

ensemble = ensemble_init(dt,ne,M,N,8,X_start);

enkfpo3r(dt,ensemble,M,N,H,t_final,R,Y,jump,4.8,0.08,spy)

enkfsr3r(dt,ensemble,M,N,H,t_final,R,Y,jump,5,0.28,spy)

enkfsr4r(dt,ensemble,M,N,H,t_final,R,Y,jump,4,0.24,spy)

enkfsr5r(dt,ensemble,M,N,H,t_final,R,Y,jump,3.6,0.27,spy)









