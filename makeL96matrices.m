function [M,N,H] = makeL96matrices(n)

v1 = ones(n-2,1);
v2 = ones(n-1,1);
M = diag(v2,1) - diag(v1,-2);
M(n,1)=1; M(1,n-1)=-1; M(2,n)=-1;
N=diag(v2,-1);
N(1,n)=1;

H = zeros(n/2,n);

for i=1:n/2
   H(i,2*i-1)=1; 
end

end

