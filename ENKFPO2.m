function [X,mu_a,P_a] = ENKFPO2(x,dt,n,Ne,H,R,r,alpha,y,M1,M2,M3,Q)

[nx,mx] = size(x);
L = localize2(n,r);
m = size(y,1);
Rm = R*eye(m);

X = lorenz63s4(x,dt,M1,M2,M3)+sqrt(dt)*sqrt(Q).*randn(nx,mx);

mu_f = mean(X,2);                                   % forecast mean
X = mu_f + sqrt(1+alpha).*(X-mu_f);                 % ensemble inflation
P = cov(X');                                        % forecast covariance
P = L.*P;                                           % localization
K = P*(H')/(H*P*H' + Rm);                           % Kalman Gain
Y_tilde = y+sqrt(R).*randn(m,Ne);                   % analysis perturbations
X = X + K*(Y_tilde-H*X);                            % analysis ensemble
mu_a = mu_f + K*(y-H*mu_f);                         % analysis mean
P_a = (eye(n)-K*H)*P;                               % analysis covariance

end

