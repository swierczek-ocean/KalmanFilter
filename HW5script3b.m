tic()
colors
eps = 0.1;
Ne = 100000;
nexp = 100;
n = 10:10:1000;
sz = size(n,2);
sigsq = 1+eps;
sig = sqrt(sigsq);
rho = zeros(sz,1);
rhoTemp = zeros(nexp,1);


for ii=1:sz   
    for jj=1:nexp
        x = normrnd(0,sig,Ne,n(ii));
        W = zeros(Ne,1);
        Wsq = zeros(Ne,1);
        for kk=1:Ne
            W(kk) = multn(x(kk,:),0,1)/multn(x(kk,:),0,sigsq);
            Wsq(kk) = W(kk)^2;
        end
        X = sum(Wsq)/Ne;
        Y = (sum(W)/Ne)^2;
        rhoTemp(jj) = X/Y;
    end
    rho(ii) = mean(rhoTemp);
end

plot(n,rho,'Color',Color(:,19),'Linewidth',2)
title('n vs. rho')
xlabel('n')
ylabel('rho')
print('prior_vs_posterior_vs_obs at time t','-djpeg')
output = [n',rho];
save output

toc()