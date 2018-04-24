function [X,mu_a,spread] = LETKF(X,dt,n,Ne,H,M1,M2,M3,alpha,y,L,Rm,nobs,sqrtQ)
[nx,mx] = size(X);
X = lorenz63s4(X,dt,M1,M2,M3)+sqrt(dt).*sqrtQ*randn(nx,mx);    % forecast ensemble

P_f = (1+alpha)*L.*cov(X');
mu_f = mean(X,2);                                   % forecast mean
X_f = (X - mu_f);                                   % forecast perturbations
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

end

