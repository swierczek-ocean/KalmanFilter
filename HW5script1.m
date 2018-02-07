tic()

format long g
Ne = [1000000:1000000:15000000];
n = 1000;
Nesz = size(Ne,2);

% funp = @(x)(exp(-(x^2)/8)/sqrt(8*pi));
% funf = @(x)(1*(x>=4));
% funq = @(x)(exp(-(x^2)/2)/sqrt(2*pi));

% Efx = integral(funp,4,Inf);
Efx = 0.0227501319481792;
e = zeros(Nesz,1);
Output1 = [];

for ii=1:Nesz
    X = zeros(Ne(ii),1);
    evec = zeros(n,1);
    for jj=1:n
        x = randn(Ne(ii),1);
        parfor kk=1:Ne(ii)
            if(x(kk)<4)
                X(kk) = 0;
            else 
                X(kk) = divgauss1(0,2,0,1,x(kk));
            end
        end
        Ehat = sum(X)/Ne(ii);
        evec(jj) = abs(Efx-Ehat)/abs(Efx);
    end
    e(ii) = sum(evec)/n;
    Output1 = [Output1;Ne(ii),e(ii)*100]
    
end

save Output1

toc()


