function [X,mu_a,P_a] = ENKFPO(X,dt,jump,n,ne,H,R,L1,L2,F,r,alpha,y)

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
X = mu_f + sqrt(1+alpha).*(X-mu_f);                 % ensemble inflation
P = cov(X');                                        % forecast covariance
P = L.*P;                                           % localization
K = P*(H')/(H*P*H' + Rm);                           % Kalman Gain
Y_tilde = y+normrnd(0,R,m,ne);                      % analysis perturbations
X = X + K*(Y_tilde-H*X);                            % analysis ensemble
mu_a = mu_f + K*(y-H*mu_f);                         % analysis mean
P_a = (eye(n)-K*H)*P;                               % analysis covariance

end

