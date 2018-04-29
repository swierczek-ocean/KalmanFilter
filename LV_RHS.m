function dX = LV_RHS(X,alpha,beta,gamma,delta)
dx = alpha*X(1)-beta*X(1)*X(2);
dy = -gamma*X(2)+delta*X(1)*X(2);
dX = [dx;dy];
% dX = [alpha,0;0,-gamma]*X + ([-beta,0;0,delta]*X).*flipud(X);
end

