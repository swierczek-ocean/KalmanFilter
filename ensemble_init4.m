function ensemble = ensemble_init4(X_start,L1,L2,F,dt,ne,n)

num_iter = 3000/dt;
rando = randperm(num_iter,ne);
X = X_start;
f = @(x)(L1*x.*(L2*x)+(F-x));
ensemble = zeros(n,ne);
t=1;

for i=1:num_iter
    k1 = f(X);
    k2 = f(X+0.5*dt.*k1);
    k3 = f(X+2*dt.*k2 -dt.*k1);
    X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);
   if(sum(i==rando)==1)
       ensemble(:,t) = X;
       t=t+1;
   end
end

end
