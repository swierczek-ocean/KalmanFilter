tic()

%% preliminaries
Ne = 100000;            % number of samples
y = 2.5;                % observation from actual people
W = zeros(Ne,1);        % initialize empty weight vector
funF = @(x)(0.5*((x-1)^2 + ((y-x^3)/0.1)^2));           % F(x)
funp = @(x)exp(-0.5.*((x-1).^2 + ((y-x.^3)./0.1).^2));  % p(x)
%%

%% optimization
[mu,phi] = fminsearch(funF,1.3);         % minimization
dx = 0.001;                              % evaluate Hessian using finite difference grid
sigsq = (funF(mu-3*dx)/90-3*funF(mu-2*dx)/20+1.5*funF(mu-dx)-49*funF(mu)/18 ...
    +1.5*funF(mu+dx)-3*funF(mu+2*dx)/20+funF(mu+3*dx)/90)/(dx^2);
stdv = sqrt(1/sigsq);
%%

%% for histograms
colors                      % my own color palette for plotting
n = 80;                     % number of bins
nbins = ceil(Ne/21);        % vertical scale
lb = 1.2;                   % x axis lower bound
ub = 1.5;                   % x axis upper bound
%%

%% for plotting p(x)
Z = unifrnd(lb,ub,Ne,1);         % uniform random variables for plotting target p(x)
Z = sort(Z);
P = funp(Z);
C = quadgk(funp,-Inf,Inf);       % scaling constant
P = P./C;
%%

%% calculating weights and effective sample size
X = stdv.*randn(Ne,1) + mu;         % Ne draws from our proposal distribution
for ii=1:Ne
    lump = -0.5*((X(ii)-1)^2 + ((y-X(ii)^3)/0.1)^2 -((X(ii)-mu)^2/(2*stdv^2)));
    W(ii) = exp(lump);              % calculate weights
end
rhonum = sum(W.^2)/Ne;              % numerator for rho calculation
rhoden = (sum(W)/Ne)^2;             % denominator for rho calculation
rho = rhonum/rhoden                 % rho
Neff = floor(Ne/rho)                % effective sample size
%%

%% resampling
x = X;                              % little x will be resampled ensemble
What = W./(sum(W));                 % W hat are the self normalized weights
Whatsum = zeros(Ne,1);                  

for ii=1:Ne
    Whatsum(ii) = sum(What(1:ii));  % Whatsum is a cumulative sum of W hat
end

U = unifrnd(0,1,Ne,1);
U = sort(U);

for jj=1:Ne                         % performs resampling algorithm
    kk=find(Whatsum>=U(jj),1);           
    x(jj) = X(kk);    
end
%%

%% plots
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
%%

toc()