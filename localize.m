function L = localize(P,threshold,r)
[n,m] = size(P);
L = eye(n);
E = zeros(n);

for i=1:n-1
   nrm = sum(abs(diag(P,i)))/(n-i);
   if(nrm>threshold)
       L = L + diag(ones(n-i,1),i) + diag(ones(n-i,1),-i);
   end
end

for i=1:n
    for j=i+1:n
        E(i,j) = exp(-abs(i-j)/r) + exp(-abs(n-j+i)/r);
    end
end

mx = max(E(:));
E = E./mx;
E = eye(n) + E + E';

L = L.*E;

end

