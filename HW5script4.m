tic()
colors
Ne = 100000;
n = 80;
ub = 4500;
nexp = 1;
rhoTemp = zeros(nexp,1);

funf = @(x)(1*(x>=4));
funp = @(x)(0.7*exp(-(x.^2)/2)/sqrt(2*pi) + 0.3*exp(-(x-4).^2/2)/sqrt(2*pi));
P = [normrnd(0,1,0.7*Ne,1); normrnd(4,1,0.3*Ne,1)];

for ll=1:nexp

% proposal #1: sum of Gaussians with slightly wider support than p, rho ~ 1.12
% funq = @(x)(0.7*exp(-(x.^2)/4)/sqrt(4*pi) + 0.3*exp(-(x-4).^2/4)/sqrt(4*pi));
% X = [normrnd(0,sqrt(2),0.7*Ne,1); normrnd(4,sqrt(2),0.3*Ne,1)];

% proposal #2: single Gaussian with mean 2 and large tails, rho ~ 1.52
% funq = @(x)(exp(-(x-2).^2/12.5)/sqrt(12.5*pi));
% X = normrnd(2,2.5,Ne,1);

% proposal #3: Gaussian that has support away from where p is large, rho ~ 2000
% funq = @(x)(exp(-(x+2).^2/2)/sqrt(2*pi));
% X = normrnd(-2,1,Ne,1);

% proposal #4: Gaussian that has support away from where p is large, rho ~ 350
% funq = @(x)(exp(-(x+1).^2/4)/sqrt(4*pi));
% X = normrnd(-1,sqrt(2),Ne,1);

% proposal #5: Gaussian that has support away from where p is large, rho ~ 6000
% funq = @(x)(exp(-(x-2).^2/.08)/sqrt(.08*pi));
% X = normrnd(2,.2,Ne,1);

% proposal #6: Uniform random variables on (-4,8), rho ~ 1.99
% funq = @(x)(1/12);
% X = unifrnd(-4,8,Ne,1);



W = funp(X)./funq(X);
rhonum = sum(W.^2)/Ne;
rhoden = (sum(W)/Ne)^2;
rhoTemp(ll) = rhonum/rhoden;

end

rho = mean(rhoTemp)

x = X;
What = W./(sum(W));
Whatsum = zeros(Ne,1);

for ii=1:Ne
    Whatsum(ii) = sum(What(1:ii));    
end

U = unifrnd(0,1,Ne,1);
U = sort(U);

for jj=1:Ne
    kk=find(Whatsum>U(jj),1);
    x(jj) = X(kk);    
end

figure()
hist(P,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,8);
%h.EdgeColor = 'w';
axis([-4 8 0 ub])
title('target distribution')
xlabel('x')
ylabel('count')

figure()
hist(X,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,7);
%h.EdgeColor = 'w';
axis([-4 8 0 ub])
title('proposal distribution')
xlabel('x')
ylabel('count')

figure()
hist(x,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,9);
%h.EdgeColor = 'w';
axis([-4 8 0 ub])
title('resampled ensemble')
xlabel('x')
ylabel('count')



toc()
