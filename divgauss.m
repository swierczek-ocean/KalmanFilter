function W = divgauss(mup,stdp,muq,stdq,X)
% X must be a row vector
% only works if covariances are diagonal

n = size(X,2);

if(stdp==1)
    lump = 0.5*((X-muq)*(X-muq)'/stdq - (X-mup)*(X-mup)');
    W = sqrt(stdq^n)*exp(lump);
elseif(stdq==1)
    lump = 0.5*((X-muq)*(X-muq)' - (X-mup)*(X-mup)'/stdp);
    W = sqrt(1/stdp^n)*exp(lump);
else
    lump = 0.5*((X-muq)*(X-muq)'/stdq - (X-mup)*(X-mup)'/stdp);
    W = sqrt((stdq/stdp)^n)*exp(lump);

end

