tic()

colors                  % my own color palette for plotting
Ne = 100000;            % number of samples
n = 80;                 % number of bins in my histogram plot
ub = 4500;              % vertical scale for histogram
y = 2.5;                % observation from actual people
std = 0.5;              % guess for Gaussian proposal standard deviation
mu = 1.357;             % guess for Gaussian proposal mean
nexp = 1;               % number of experiments

rhoTemp = zeros(nexp,1);        % vector of rhos for each experiment
W = zeros(Ne,1);                % initialize empty weight vector
X = std.*randn(Ne,1) + mu;      % Ne draws from our proposal distribution

%% for plotting p(x)
Z = unifrnd(-4,8,Ne,1);         % uniform random variables for plotting target p(x)
Z = sort(Z);
funp = @(x)exp(-0.5.*((x-1).^2 + ((y-x.^3)./0.1).^2));  % p(x)
P = funp(Z);
size(P)
C = quadgk(funp,-Inf,Inf)       % scaling constant
P = P./C;
%%

for ll=1:nexp
    for ii=1:Ne
        lump = -0.5*((X(ii)-1)^2 + ((y-X(ii)^3)/0.1) -((X(ii)-mu)^2/(2*std^2)));
        W(ii) = sqrt(2*pi)*std*exp(lump)/C;  
    end
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
plot(Z,P,'Color',Color(:,19),'Linewidth',2)
%h = findobj(gca,'Type','patch');
%h.FaceColor = Color(:,8);
%h.EdgeColor = 'w';
axis([-4 8 0 25])
title('target distribution')
xlabel('x')
ylabel('p(x)')

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