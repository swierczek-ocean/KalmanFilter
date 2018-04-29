function X_obs = run_LV(theta,dt)
X_obs = [theta(5:6),zeros(2,10)];
obs = 2;
x = theta(5:6);
for ii=1:100
    x = LV_RK4(x,theta(1),theta(2),theta(3),theta(4),dt);
    if mod(ii,10)==0
        X_obs(:,obs) = x;
        obs = obs+1;
    end
end
end

