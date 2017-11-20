function ensemble = ensemble_init4(X_start)
global L1, global L2, global F, global dt, global ne, global n

num_iter = 3000/dt;
rando = randperm(num_iter,ne);
X = X_start;
f = @(x)(L1*x.*(L2*x)+(F-x));
ensemble = zeros(n,ne);
t=1;

for i=1:num_iter
   k = f(X) + f(X+dt.*f(X)); 
   X = X + 0.5*dt.*k;
   if(sum(i==rando)==1)
       ensemble(:,t) = X;
       t=t+1;
   end
end

end
