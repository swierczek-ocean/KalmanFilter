tic()

%% preliminaries
Ne = 100000;            % number of samples
y = 2.5;                % observation from actual people
std = 0.015;            % guess for Gaussian proposal standard deviation, 0.015 is *good*
mu = 1.357208808;       % guess for Gaussian proposal mean
nexp = 30;               % number of experiments
rhoTemp = zeros(nexp,1);        % vector of rhos for each experiment
W = zeros(Ne,1);                % initialize empty weight vector
%%

%% for histograms
colors                  % my own color palette for plotting
n = 80;                 % number of bins
nbins = 4500;           % vertical scale
lb = 1.2;               % x axis lower bound
ub = 1.5;               % x axis upper bound
%%

%% for plotting p(x)
Z = unifrnd(lb,ub,Ne,1);         % uniform random variables for plotting target p(x)
Z = sort(Z);
funp = @(x)exp(-0.5.*((x-1).^2 + ((y-x.^3)./0.1).^2));  % p(x)
P = funp(Z);
C = quadgk(funp,-Inf,Inf);       % scaling constant
P = P./C;
%%

for ll=1:nexp
    X = std.*randn(Ne,1) + mu;      % Ne draws from our proposal distribution
    for ii=1:Ne                         
        
        lump = -0.5*((X(ii)-1)^2 + ((y-X(ii)^3)/0.1) -((X(ii)-mu)^2/(2*std^2)));
        W(ii) = sqrt(2*pi)*std*exp(lump)/C;  % calculate weights
    end
    rhonum = sum(W.^2)/Ne;              % numerator for rho calculation
    rhoden = (sum(W)/Ne)^2;             % denominator for rho calculation
    rhoTemp(ll) = rhonum/rhoden;
end

rhoTemp
rho = mean(rhoTemp)
Neff = floor(Ne/rho)

x = X;                                  % little x will be resampled ensemble
What = W./(sum(W));                     % W hat are the self normalized weights
Whatsum = zeros(Ne,1);                  

for ii=1:Ne
    Whatsum(ii) = sum(What(1:ii));      % Whatsum is a cumulative sum of W hat
end

U = unifrnd(0,1,Ne,1);
U = sort(U);

for jj=1:Ne                             % performs resampling algorithm
    kk=find(Whatsum>=U(jj),1);           
    x(jj) = X(kk);    
end

figure()
plot(Z,P,'Color',Color(:,19),'Linewidth',2)
axis([lb ub 0 25])
title('target distribution')
xlabel('x')
ylabel('p(x)')

figure()
hist(X,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,7);
axis([lb ub 0 nbins])
title('proposal distribution')
xlabel('x')
ylabel('count')

figure()
hist(x,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,9);
axis([lb ub 0 nbins])
title('resampled ensemble')
xlabel('x')
ylabel('count')


toc()