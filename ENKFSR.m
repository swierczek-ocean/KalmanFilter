function [X,mu_a,spread] = ENKFSR(X,dt,jump,n,Ne,H,M,N,F,alpha,y,L,Rm,nobs)

for i=1:jump
    X = lorenz96s4(X,dt,M,N,F);    % forecast ensemble
end

P_f = (1+alpha)*L.*cov(X');
mu_f = mean(X,2);                                   % forecast mean
X_f = (X - mu_f);                                   % forecast perturbations
%X_t = mu_f + sqrt(1+alpha)*X_f;
K = P_f*H'*((H*P_f*H'+Rm)\eye(nobs));
mu_a = mu_f +K*(y-H*mu_f);

Z = (sqrt(1+alpha)/sqrt(Ne-1))*X_f;
tmp = Z'*H'*(Rm\(H*Z));
tmp = .5*(tmp+tmp');
[E,OM] = eig(tmp);
T = E*(sqrtm(eye(Ne)+OM)\E');
Za = Z*T;


X = mu_a + sqrt(Ne-1)*Za;
P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
spread = sqrt(trace(P_a)/n);

% mu_f = mean(X,2);                                   % forecast mean
% x_f = mu_f + sqrt(1+alpha).*(X-mu_f);               % ensemble inflation
% X_f = (x_f - mu_f).*(1/sqrt(ne-1));                 % forecast perturbations
% P_f = cov(X_f');                                    % forecast covariance
% P_f = L.*P_f;                                       % localization
% K = P_f*(H')/(H*P_f*H' + Rm);                       % Kalman Gain
% mu_a = mu_f + K*(y-H*mu_f);                         % analysis mean
% P_a = (eye(n)-K*H)*P_f;                             % analysis covariance
% V = (X_f')*(H');
% A = (V/Rm)*(V');
% A = 0.5.*(A+A');
% [U,lambda] = eig(A);
% Q = real(sqrtm(eye(ne)+lambda));
% Z = U/Q;
% X_a = X_f*Z;                                        % analysis perturbations
% X = sqrt(ne-1).*X_a + mu_a;                         % analysis ensemble

end

