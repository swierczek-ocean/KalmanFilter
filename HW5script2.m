tic()
Ne = 3000;
n = 1000;
Nesz = size(Ne,2);

funp = @(x)(exp(-x^2/8)/sqrt(8*pi));
funf = @(x)(1*(x>=4));
funq = @(x)(exp(-(x-2)^2/2)/sqrt(2*pi));

% Efx = integral(funp,4,Inf);
Efx = 0.0227501319481792;
e = zeros(Nesz,1);

for ii=1:Nesz
    X = zeros(Ne(ii),1);
    evec = zeros(n,1);
    for jj=1:n
        x = normrnd(2,1,Ne(ii),1);
        for kk=1:Ne(ii)
           X(kk) = funf(x(kk))*funp(x(kk))/funq(x(kk));              
        end
        Ehat = sum(X)/Ne(ii);
        evec(jj) = abs(Efx-Ehat)/abs(Efx);
    end
    e(ii) = sum(evec)/n;
end

format long g

Output = [Ne',e*100]
toc()