function x = multn(X,mu,sigma)
%X should be a row vector, sigma a scalar, mu a scalar or row vector
n = size(X,2);
% Q = eye(n)./sigma;
% x = exp(-(X-mu)*Q*(X-mu)')./sqrt(2*pi*(sigma^n));
if(sigma==1)
   x = exp(-(X-mu)*(X-mu)')/sqrt(2*pi);
else
   x = exp(-(X-mu)*(X-mu)'/sigma)/sqrt(2*pi*(sigma^n)); 
end
end

