function [R,U] = house_qr(A)
    % Householder reflections for QR decomposition.
    % [R,U] = house_qr(A) returns
    % R, the upper triangular factor, and
    % U, the reflector generators for use by house_apply.    
    H = @(u,x) x - u*(u'*x);
    [m,n] = size(A);
    U = zeros(m,n);
    R = A;
    for j = 1:min(m,n)
        u = house_gen(R(j:m,j));
        U(j:m,j) = u;
        R(j:m,j:n) = H(u,R(j:m,j:n));
        R(j+1:m,j) = 0;
    end
end

