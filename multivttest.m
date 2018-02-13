a=10;
b=0;

[X1,X2] = meshgrid(linspace(-5,15,200)',linspace(-5,15,200)');
X = [X1(:)-4.88 X2(:)-4.88];
C = [a,b;b,a];
df = 3;
p = mvtpdf(X,C,df);

figure;
surf(X1,X2,reshape(p,200,200))