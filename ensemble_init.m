function ensemble = ensemble_init(dt,ne,M,N,F,X_start)
num_iter = 1000/dt;
rando = randperm(num_iter,ne);
X = X_start;
ensemble = [];

for i=1:num_iter
   X = lorenz96s4(X,dt,M,N,F);
   if(sum(i==rando)==1)
       ensemble = [ensemble,X];
   end
end
end

