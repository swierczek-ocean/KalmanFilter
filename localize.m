function L = localize(P,r)
[n,m] = size(P);
L = eye(n);

for i=1:n-1
   nrm = sum(abs(diag(P,i)))/(n-i);
   if(nrm>r)
       L = L + diag(ones(n-i,1),i) + diag(ones(n-i,1),-i);
   end
end

end

