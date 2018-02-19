%%
clear 
close clc


Ne = 1000000000


tic()

X = 2*randn(Ne,1)+5;

toc()

tic()

x = normrnd(5,2,Ne,1);

toc()


