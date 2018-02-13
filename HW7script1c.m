tic()

%% preliminaries
Ne = 100;            % number of samples
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
    eqn = funF(mu + lambda*L*xi) - phi == 0.5*(xi'*xi);
    sollambda = vpasolve(eqn,lambda);
    eqnhat = funF(mu + lambda*L*xi) - phi == 0.5*(xi'*xi + dx);
    sollambdahat = vpasolve(eqnhat,lambda);
    dldr = (sollambdahat(1) - sollambda(1))/dx;
    X(ii) = mu + sollambda(1)*L*xi;                    % calculate sample
    W(ii) = sollambda(1)+2*(xi'*xi)*dldr;              % calculate weights
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
print('target','-djpeg')

figure()
hist(X,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,7);
axis([lb ub 0 nbins])
title('proposal distribution')
xlabel('x')
ylabel('count')
print('proposal','-djpeg')

figure()
hist(x,n)
h = findobj(gca,'Type','patch');
h.FaceColor = Color(:,9);
axis([lb ub 0 nbins])
title('resampled ensemble')
xlabel('x')
ylabel('count')
print('resampled','-djpeg')
%%

toc()