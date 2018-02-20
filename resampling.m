function x = resampling(What,X,Ne,n)
%%
%  Resamping (for SMC)
% (see "A Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian
%       Bayesian Tracking," Arulampalam, Maskell, Gordon and Clapp
%       IEEE Trans. on Signal Processing, Vol. 50, No.2, 2002
%

% construct the cdf
% c = zeros(Ne,1);
% for jj=2:Ne+1
%     c(jj)=c(jj-1)+What(jj-1);
% end, clear jj

c = [0;cumsum(What(1:end-1))];

% c = cumsum(What);

% sample it and get the stronger particle more often
ii=1; % initialize
x=zeros(n,Ne);
u1=rand/Ne; % initialize
for jj=1:Ne-1
    u = u1+(jj-1)/Ne;
    while u>c(ii)
        ii=ii+1;
    end
    x(:,jj) = X(:,ii);
end
