dt = 0.01;
t_final = 10;
R = 1;
obsdt = 0.1;
threshold = [0.05,0.1];
r = [1:0.5:2];
alpha = [0:0.1:0.2];
set(groot, 'DefaultFigureVisible', 'off');

[M,N,H,SynthDataTrue,SynthDataObs,X_start,jump] = lorenz(40,dt,t_final,8,R,obsdt);


T = SynthDataTrue;
Y = SynthDataObs;
size(Y);
size(T);

ne=40;
ensemble = ensemble_init(dt,ne,M,N,8,X_start);
W = randomrotation(ne);


[o,rsz] = size(r);
[o,asz] = size(alpha);
[o,tsz] = size(threshold);
POsearch = zeros(rsz,asz,tsz);
SRRR1search = zeros(rsz,asz,tsz);
SRRRsearch = zeros(rsz,asz,tsz);
SRsearch = zeros(rsz,asz,tsz);

for i=1:rsz
    for j=1:asz
        for k=1:tsz

            [ARMSE,aspread] = enkfpo3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            POsearch(i,j,k) = ARMSE;
            
            [ARMSE,aspread] = enkfsr3(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            SRRR1search(i,j,k) = ARMSE;
            
            [ARMSE,aspread] = enkfsr4(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            SRRRsearch(i,j,k) = ARMSE;
            
            [ARMSE,aspread] = enkfsr5(dt,ensemble,M,N,H,t_final,R,Y,T,jump,threshold(k),r(i),alpha(j));
            SRsearch(i,j,k) = ARMSE;

        end
    end
end

[minPO,I1] = min(POsearch(:));
[i1,j1,k1] = ind2sub([rsz,asz,tsz],I1);
[minSRRR1,I2] = min(SRRR1search(:));
[i2,j2,k2] = ind2sub([rsz,asz,tsz],I2);
[minSRRR,I3] = min(SRRRsearch(:));
[i3,j3,k3] = ind2sub([rsz,asz,tsz],I3);
[minSR,I4] = min(SRsearch(:));
[i4,j4,k4] = ind2sub([rsz,asz,tsz],I4);

Result = [minPO,r(i1),alpha(j1),threshold(k1)];
Result = [Result;minSRRR1,r(i2),alpha(j2),threshold(k2)];
Result = [Result;minSRRR,r(i3),alpha(j3),threshold(k3)];
Result = [Result;minSR,r(i4),alpha(j4),threshold(k4)];

save POsearch
save SRRR1search
save SRRRsearch
save SRsearch
save Result