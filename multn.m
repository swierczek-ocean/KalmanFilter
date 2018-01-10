function x = multn(X,mu,sigma)
%X should be a row vector, sigma a scalar, mu a scalar or row vector
n = size(X,2);
Q = eye(n)./sigma;
x = exp(-(X-mu)*Q*(X-mu)')./sqrt(2*pi*(sigma^n));
end

