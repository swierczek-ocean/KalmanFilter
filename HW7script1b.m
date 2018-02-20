
clear
close all
clc
format long g

tic()
%% preliminaries
Ne = 5e5;            % number of samples
y = 2.5;                % observation from actual people
W = zeros(Ne,1);        % initialize empty weight vector
funF2 = @(x)[(x-1)/sqrt(2); (y-x^3)/0.1];
funF = @(x)(0.5*((x-1)^2 + ((y-x^3)/0.1)^2));           % F(x)
funp = @(x)exp(-0.5.*((x-1).^2 + ((y-x.^3)./0.1).^2));  % p(x)
nu = 1:15;
nexp = size(nu,2);
%%

%% optimization
[mu,~,~,~,~,~,J] = lsqnonlin(funF,1.3);         % minimization
dx = 0.001;                              % evaluate Hessian using finite difference grid
sigsq = (funF(mu-3*dx)/90-3*funF(mu-2*dx)/20+1.5*funF(mu-dx)-49*funF(mu)/18 ...
    +1.5*funF(mu+dx)-3*funF(mu+2*dx)/20+funF(mu+3*dx)/90)/(dx^2);
stdv = 1/sqrt(sigsq);
%%

%% for histograms
colors                      % my own color palette for plotting
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
rhovec = zeros(nexp,2);
for ii=1:nexp
    X = zeros(Ne,1);
    W = zeros(Ne,1);
    for kk=1:Ne
        h = stdv*randn;
        u = norm(randn(nu(ii),1))^2;
        z = mu+h/sqrt(u/nu(ii));
        W(kk) = funF(z)-neglog_q_t(z,mu,stdv,nu(ii));
        X(kk) = z;
    end
    What = normalizeweights(W);
    
    rhonum = mean(What.^2);             % numerator for rho calculation
    rhoden = (mean(What))^2;            % denominator for rho calculation
    rhovec(ii,1) = rhonum/rhoden                 % rho
    rhovec(ii,2) = floor(Ne/rhovec(ii,1))                % effective sample size

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
title(['target vs. proposal vs. resampled ensemble for nu =',num2str(nu(ii))])
xlabel('x')
ylabel('y')
legend([h1(1),h2(1),h3(1)],'resampled ensemble','proposal ensemble','target distribution')
hold off
%%
end
toc()
