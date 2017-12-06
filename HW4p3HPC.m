tic();
set(groot, 'DefaultFigureVisible', 'off');
colors
bg=3000;
n=40;
dt=0.01;
obsdt=0.2;
ne=20;
jump = ceil(obsdt/dt);
R=1;
t_final=60;
spy=8;
spl=250;
Rm = R*eye(ne);
F = 8;
[L1,L2,H] = prelim(n);

[SynthDataTrue,SynthDataObs,X_start] = lorenz3(n,t_final,L1,L2,H,F,dt,jump,R);

T = SynthDataTrue(:,1:spl*jump);
Y = SynthDataObs(:,1:spl);

T2 = SynthDataTrue(:,spl*jump+1:end);
Y2 = SynthDataObs(:,spl+1:end);

j1 = size(T2,2);
j2 = size(Y2,2);
RMSE = zeros(1,j2-1);
spread = zeros(1,j2-1);


r=2.6:0.2:4.6;
alpha=0:0.025:0.2;

rsz = max(size(r));
asz = max(size(alpha));

tsz = rsz*asz;
summary3 = zeros(tsz,3);
counter = 1;

for kk=1:rsz
    for jj=1:asz
        
        ensemble = ensemble_init4(X_start,L1,L2,F,dt,ne,n);
        [ARMSE,aspread,X_a,mu_a,P_a] = enkfpo4(ensemble,Y,T,r(kk),alpha(jj),spy,L1,L2,H,R,dt,jump,n,ne);
        ARMSE
        aspread
        
        X=X_a;
        
        bmean=mu_a;
        bcov=P_a;
        spyvec = zeros(1,j2-1);
        indices = zeros(1,j2-1);
        L = localize2(n,r(kk));
        
        for i=1:j2-1
            y_t = Y2(:,i+1);
            [xT,x0] = fdvar(bmean,dt,jump,L1,L2,H,bcov,bmean,y_t,n,R,F);
            X(:,1) = xT;
            error = xT-T2(:,jump*i+1);
            RMSE(i) = sqrt((1/n).*transpose(error)*error);
            bcov = 0.5.*(bcov +bcov');
            bcov = sqrtm(bcov);
            for j=2:ne
                PertMean = bmean +bcov*normrnd(0,1,n,1);
                PertObs = y_t + eye(n/2)*normrnd(0,1,n/2,1);
                [xT,x0] = fdvar(PertMean,dt,jump,L1,L2,H,bcov,bmean,PertObs,n,R,F);
                X(:,j) = xT;
            end
            bmean = X(:,1);
            EnMean = mean(X,2);
            X = EnMean + sqrt(1+alpha(jj)).*(X-EnMean);
            bcov = cov(X');
            bcov = L.*bcov;
            spread(i) = sqrt(trace(bcov)/n);
%             spyvec(i) = xT(spy);
%             indices(i) = i*jump+1;
        end
        
        average_RMSE = mean(RMSE(floor(end/2):end))
        
        summary3(counter,:) = [average_RMSE,r(kk),alpha(jj)];
        counter = counter + 1;
    end
end

save summary3

toc()