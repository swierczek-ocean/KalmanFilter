function [L1,L2,H] = prelim(n)

v1 = ones(n-2,1);
v2 = ones(n-1,1);
L1 = diag(v2,1) - diag(v1,-2);
L1(n,1)=1; L1(1,n-1)=-1; L1(2,n)=-1;
L2=diag(v2,-1);
L2(1,n)=1;
H = zeros(n/2,n);

for i=1:n/2
   H(i,2*i-1)=1; 
end

end

