function [l,x,w]=MyRandomMap(y,m,s,phi,xi)
func=@(l) rm_fun(l,m,s,phi,xi,y);
options = optimoptions('fsolve','Display','none');
l = fsolve(func,1,options);
x = m+l*s*xi;

rho = xi^2;

dx = 1e-3;
gradF = (F(x+dx,y)-F(x,y))/dx;
dldrho = 1/(gradF*xi);
w = l+2*rho*dldrho;

