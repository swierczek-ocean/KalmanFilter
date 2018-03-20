function U = randrot(n)
x = ones(n,1)/sqrt(n);
Omega = 5*randn(n-1,n-1);
[~,~,V] = svd(Omega);
Omega2 = 5*randn(n,n);
Omega2(:,1)=x;
A = GramSchmidt(Omega2);
Lambda = zeros(n,n);
Lambda(2:end,2:end) = V';
Lambda(1,1) = 1;
U = A*Lambda*A';
end

