dt = 0.01;
t_final = 10;
R = 1;
obsdt = 0.1;
r = 0.5:0.5:4;
alpha = 0:0.05:1;
set(groot, 'DefaultFigureVisible', 'off');

[M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(40,dt,t_final,8,R,obsdt);


T = SynthDataTrue;
Y = SynthDataObs;
size(Y);
size(T);

ne=20;
ensemble = ensemble_init(dt,ne,M,N,8,X_start);
W = randomrotation(ne);


[o,rsz] = size(r);
[o,asz] = size(alpha);
POsearche = ones(rsz,asz);
SRRR1searche = ones(rsz,asz);
SRRRsearche = ones(rsz,asz);
SRsearche = ones(rsz,asz);
POsearchs = ones(rsz,asz);
SRRR1searchs = ones(rsz,asz);
SRRRsearchs = ones(rsz,asz);
SRsearchs = ones(rsz,asz);

for i=1:rsz
    for j=1:asz


            [ARMSE,aspread] = enkfpo3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r(i),alpha(j));
            POsearche(i,j) = ARMSE;
            POsearchs(i,j) = aspread;
            
            [ARMSE,aspread] = enkfsr3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r(i),alpha(j));
            SRRR1searche(i,j) = ARMSE;
            SRRR1searchs(i,j) = aspread;
            
            [ARMSE,aspread] = enkfsr4(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r(i),alpha(j));
            SRRRsearche(i,j) = ARMSE;
            SRRRsearchs(i,j) = aspread;
            
            [ARMSE,aspread] = enkfsr5(dt,ensemble,M,N,H,t_final,R,Y,T,jump,r(i),alpha(j));
            SRsearche(i,j) = ARMSE;
            SRsearchs(i,j) = aspread;

    end
end

[minPO,I1] = min(POsearche(:));
[i1,j1] = ind2sub([rsz,asz],I1);
[minSRRR1,I2] = min(SRRR1searche(:));
[i2,j2] = ind2sub([rsz,asz],I2);
[minSRRR,I3] = min(SRRRsearche(:));
[i3,j3] = ind2sub([rsz,asz],I3);
[minSR,I4] = min(SRsearche(:));
[i4,j4] = ind2sub([rsz,asz],I4);

Result = [minPO,r(i1),alpha(j1),POsearchs(i1,j1)];
Result = [Result;minSRRR1,r(i2),alpha(j2),SRRR1searchs(i2,j2)];
Result = [Result;minSRRR,r(i3),alpha(j3),SRRRsearchs(i3,j3)];
Result = [Result;minSR,r(i4),alpha(j4),SRsearchs(i4,j4)];

save POsearche
save SRRR1searche
save SRRRsearche
save SRsearche
save POsearchs
save SRRR1searchs
save SRRRsearchs
save SRsearchs
save Result