function X = lorenz63s(x,dt,sigma,beta,rho,Q)
[n,m] = size(x);
X = zeros(n,m);

X(1,:) = x(1,:) + dt*sigma.*(x(2,:)-x(1,:));
X(2,:) = x(2,:) + dt.*(x(1,:).*(rho-x(3,:))-x(2,:));
X(3,:) = x(3,:) + dt.*(x(1,:).*x(2,:)-beta.*x(3,:));

X = X + sqrt(dt)*sqrt(Q).*randn(n,m);
end

