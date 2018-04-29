function X = LV_RK4(x,alpha,beta,gamma,delta,dt)
k1=LV_RHS(x,alpha,beta,gamma,delta);
k2=LV_RHS(x+0.5*dt.*k1,alpha,beta,gamma,delta);
k3=LV_RHS(x+0.5*dt.*k2,alpha,beta,gamma,delta);
k4=LV_RHS(x+dt.*k3,alpha,beta,gamma,delta);
X = x + (1/6)*dt.*(k1+2.*k2+2.*k3+k4);

%X = x + dt.*LV_RHS(x,alpha,beta,gamma,delta);
end

