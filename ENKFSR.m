function [X,mu_a,P_a] = ENKFSR(X,dt,jump,n,ne,H,R,L1,L2,F,r,alpha,y)

f = @(x)(L1*x.*(L2*x)+(F-x));
L = localize2(n,r);
m = size(y,1);
Rm = R*eye(m);

for i=1:jump
    k1 = f(X);
    k2 = f(X+0.5*dt.*k1);
    k3 = f(X+2*dt.*k2 -dt.*k1);
    X = X + dt.*((1/6).*k1+(2/3).*k2+(1/6).*k3);    % forecast ensemble
end

mu_f = mean(X,2);                                   % forecast mean
x_f = mu_f + sqrt(1+alpha).*(X-mu_f);               % ensemble inflation
X_f = (x_f - mu_f).*(1/sqrt(ne-1));                 % forecast perturbations
P_f = cov(X_f');                                    % forecast covariance
P_f = L.*P_f;                                       % localization
K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
mu_a = mu_f + K*(y-H*mu_f);                         % analysis mean
P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
V = (X_f')*(H');
A = (V/Rm)*(V');
A = 0.5.*(A+A');
[U,lambda] = eig(A);
Q = real(sqrtm(eye(ne)+lambda));
Z = U/Q;
%W = randomrotation(ne);
X_a = X_f*Z;                                        % analysis perturbations
X = mu_a + sqrt(ne-1).*X_a;                         % analysis ensemble

end

