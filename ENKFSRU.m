function [X,mu_a,spread] = ENKFSRU(X,dt,jump,n,Ne,H,M,N,F,alpha,y,L,Rm,nobs)

for i=1:jump
    X = lorenz96s4(X,dt,M,N,F);    % forecast ensemble
end

P_f = (1+alpha)*L.*cov(X');
mu_f = mean(X,2);                                   % forecast mean
X_f = (X - mu_f);                                   % forecast perturbations
X_t = mu_f + sqrt(1+alpha)*X_f;
K = P_f*H'*((H*P_f*H'+Rm)\eye(nobs));
mu_a = mu_f +K*(y-H*mu_f);

Z = (sqrt(1+alpha)/sqrt(Ne-1))*X_f;
tmp = Z'*H'*(Rm\(H*Z));
tmp = .5*(tmp+tmp');
[E,OM] = eig(tmp);
T = E*(sqrtm(eye(Ne)+OM)\E');
Za = Z*T;
W = randrot(Ne);
X = mu_a + sqrt(Ne-1)*Za*W;
spread = sqrt(trace(cov(X'))/n);

end

