function Mlin = Mlin2(TimeSeries,L1,L2,F,n,dt,jump)

f = @(x)(L1*x.*(L2*x)+(F-x));
df = @(x)(diag(L1*x)*L2 + diag(L2*x)*L1 - eye(n));
deriv = @(x)(eye(n) + 0.5*dt.*df(x) + 0.5*dt.*df(x+dt.*f(x))*(eye(n) + dt.*df(x)));
Mlin = eye(n);

for i=1:jump
    Temp = deriv(TimeSeries(:,i));
    Mlin = (Temp)*Mlin;
end

end

