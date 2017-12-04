function deriv = RK3deriv(x,L1,L2,F,n,dt)
f = @(x)(L1*x.*(L2*x)+(F-x));
df = @(x)(diag(L1*x)*L2 + diag(L2*x)*L1 - eye(n));

k1 = dt.*f(x);
k2 = dt.*df(x);
k3 = dt.*f(x+0.5.*k1);
k4 = dt.*df(x+0.5.*k1);
k5 = k4*(eye(n)+0.5.*k2);
k6 = x+2.*k3-k1;
k7 = dt.*df(k6);
k8 = eye(n) + 2.*k4*(eye(n)+0.5.*k2) - k2;

deriv = eye(n)+(1/6).*k2+(2/3).*k5+...
    +(1/6).*k7*k8;

end

