n=40;
dt=0.01;
obsdt=0.2;

[L1,L2,H] = prelim(n);
[SynthDataTrue,SynthDataObs,X_start,jump] = lorenz2(n,dt,t_final,F,R,obsdt,L1,L2,H)





