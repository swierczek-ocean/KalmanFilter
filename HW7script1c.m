tic()

%% preliminaries
Ne = 1e4;             % number of samples
y = 2.5;              % observation from actual people
W = zeros(Ne,1);      % initialize empty weight vector
X = zeros(Ne,1);      % initialize empty sample vector
funF = @(x)(0.5*((x-1)^2 + ((y-x^3)/0.1)^2));           % F(x)
funp = @(x)exp(-0.5.*((x-1).^2 + ((y-x.^3)./0.1).^2));  % p(x)
%%

%% optimization
[mu,phi] = fminsearch(funF,1.3);         % minimization
dx = 0.001;                              % evaluate Hessian using finite difference grid
sigsq = (funF(mu-3*dx)/90-3*funF(mu-2*dx)/20+1.5*funF(mu-dx)-49*funF(mu)/18 ...
    +1.5*funF(mu+dx)-3*funF(mu+2*dx)/20+funF(mu+3*dx)/90)/(dx^2);
L = sqrt(1/sigsq);
%%

%% for histograms
colors                  % my own color palette for plotting
n = 50;                 % number of bins
nbins = ceil(Ne/20);    % vertical scale
lb = 1.2;               % x axis lower bound
ub = 1.5;               % x axis upper bound
%%

%% for plotting p(x)
Z = unifrnd(lb,ub,Ne,1);         % uniform random variables for plotting target p(x)
Z = sort(Z);
P = funp(Z);
C = quadgk(funp,-Inf,Inf);       % scaling constant
P = P./C;
%%

%% calculating weights and effective sample size
syms lambda
for ii=1:Ne
    xi = randn;
    [l,z,w]=MyRandomMap(y,mu,L,phi,xi);
    X(ii) = z;                    % calculate sample
    W(ii) = -log(w);              % calculate weights
end
What = normalizeweights(W);

rhonum = mean(What.^2);              % numerator for rho calculation
rhoden = (mean(What))^2;             % denominator for rho calculation
rho = real(rhonum/rhoden)            % rho
Neff = floor(Ne/rho)                 % effective sample size
%%                

%% resampling
x = resampling(What,X',Ne,1);
%%

%% plot
figure()
h1 = histogram(x,'Normalization','pdf');
hold on
h2 = histogram(X,'Normalization','pdf');
h3 = plot(Z,P,'Color',Color(:,28),'Linewidth',2);
axis([lb ub 0 30])
title('random map proposal')
xlabel('x')
ylabel('y')
legend([h1(1),h2(1),h3(1)],'resampled ensemble','proposal ensemble','target distribution')
hold off
%%

toc()