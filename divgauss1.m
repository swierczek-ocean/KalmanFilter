function W = divgauss1(mup,stdp,muq,stdq,X)
% X must be a scalar or an array of scalars
[n,m] = size(X);
W = zeros(n,m);

for ii=1:n
    for jj=1:m
        lump = 0.5*(((X(ii,jj)-muq)^2)/(stdq^2) - ((X(ii,jj)-mup)^2)/(stdp^2));
        W(ii,jj) = (stdq/stdp)*exp(lump);
    end
end
end

