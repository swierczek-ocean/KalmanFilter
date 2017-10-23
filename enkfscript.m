dt = 0.01;
t_final = 10;
R = 1;
obsdt = 0.1;
threshold = [0.05:0.05:0.3];
r = [0.1:0.05:5];
alpha = [0:0.05:3];
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
[o,tsz] = size(threshold);
POsearche = zeros(rsz,asz,tsz);
SRRR1searche = zeros(rsz,asz,tsz);
SRRRsearche = zeros(rsz,asz,tsz);
SRsearche = zeros(rsz,asz,tsz);
POsearchs = zeros(rsz,asz,tsz);
SRRR1searchs = zeros(rsz,asz,tsz);
SRRRsearchs = zeros(rsz,asz,tsz);
SRsearchs = zeros(rsz,asz,tsz);

for i=1:rsz
    for j=1:asz
        for k=1:tsz

            [ARMSE,aspread] = enkfpo3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            POsearche(i,j,k) = ARMSE;
            POsearchs(i,j,k) = aspread;
            
            [ARMSE,aspread] = enkfsr3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            SRRR1searche(i,j,k) = ARMSE;
            SRRR1searchs(i,j,k) = aspread;
            
            [ARMSE,aspread] = enkfsr4(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            SRRRsearche(i,j,k) = ARMSE;
            SRRRsearchs(i,j,k) = aspread;
            
            [ARMSE,aspread] = enkfsr5(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            SRsearche(i,j,k) = ARMSE;
            SRsearchs(i,j,k) = aspread;

        end
    end
end

[minPO,I1] = min(POsearche(:));
[i1,j1,k1] = ind2sub([rsz,asz,tsz],I1);
[minSRRR1,I2] = min(SRRR1searche(:));
[i2,j2,k2] = ind2sub([rsz,asz,tsz],I2);
[minSRRR,I3] = min(SRRRsearche(:));
[i3,j3,k3] = ind2sub([rsz,asz,tsz],I3);
[minSR,I4] = min(SRsearche(:));
[i4,j4,k4] = ind2sub([rsz,asz,tsz],I4);

Result = [minPO,r(i1),alpha(j1),threshold(k1),POsearchs(i1,j1,k1)];
Result = [Result;minSRRR1,r(i2),alpha(j2),threshold(k2),SRRR1searchs(i2,j2,k2)];
Result = [Result;minSRRR,r(i3),alpha(j3),threshold(k3),SRRRsearchs(i3,j3,k3)];
Result = [Result;minSR,r(i4),alpha(j4),threshold(k4),SRsearchs(i4,j4,k4)];

save POsearch
save SRRR1search
save SRRRsearch
save SRsearch
save Result