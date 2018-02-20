function out = rm_fun(l,m,s,phi,xi,y)

x = m+l*s*xi;
out = F(x,y)-phi - .5*xi^2;
