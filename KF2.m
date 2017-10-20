function [TrueState,SynthData,mu] = KF2(dt,nsteps,xi,omega,R,x)
tic();

time = 0:dt:((nsteps-1)*dt);
A = [0,1;-omega*omega,-2*xi*omega];
M = eye(2) + dt.*A + 0.5*dt*dt.*A*A + ...
    + (1/6)*dt*dt*dt.*A*A*A + (1/24)*dt*dt*dt*dt.*A*A*A*A
H = [1,0];
eta = normrnd(0,R,1,nsteps);
KalmanGain = zeros(2,nsteps);
TrueState = [x,zeros(2,nsteps-1)];
P=eye(2);
mu = zeros(2,nsteps);
MSE = zeros(1,nsteps);
trP = [1,zeros(1,nsteps-1)];

for i=1:nsteps-1
    TrueState(:,i+1) = M*TrueState(:,i);
end

SynthData = TrueState(1,:) + eta;

figure
plot(time,TrueState(1,:),'Linewidth',3,'Color','blue')
hold on
plot(time,TrueState(2,:),'Linewidth',3,'Color','magenta')
title('True State')
xlabel('time')
legend('position','velocity')
print('TrueState','-djpeg')

figure
plot(time,TrueState(1,:),'Linewidth',3,'Color','blue')
title('True Position')
xlabel('time')
ylabel('position')
print('TrueStatePosition','-djpeg')

figure
plot(time,TrueState(2,:),'Linewidth',3,'Color','magenta')
title('True Velocity')
xlabel('time')
ylabel('velocity')
print('TrueStateVelocity','-djpeg')

figure
plot(time,TrueState(1,:),'Linewidth',3,'Color','blue')
hold on
plot(time,SynthData,'*','MarkerSize',5,'Color','red')
title('True State vs. Perturbed')
xlabel('time')
ylabel('position')
legend('true position','perturbed position')
print('TrueStatePerturbedPosition','-djpeg')

for j=1:nsteps-1
    muf=M*mu(:,j);
    P = M*P*transpose(M);
    Temp = H*P*transpose(H)+R;
    K = P*transpose(H)/Temp;
    KalmanGain(:,j) = K;
    mu(:,j+1)=muf+K*(SynthData(j+1)-H*muf);
    P = (eye(2)-K*H)*P;    
    trP(j+1) = 0.5*(P(1,1)+P(2,2));
end

error = TrueState-mu;

for i=1:nsteps
   MSE(i) = 0.5.*transpose(error(:,i))*error(:,i);
end

figure
plot(time,TrueState(1,:),'Linewidth',3,'Color','blue')
hold on
plot(time,SynthData(1,:),'*','MarkerSize',5,'Color','red')
plot(time,mu(1,:),'o','MarkerSize',5,'Color','green')
title('True State vs. Perturbed vs. Recovered')
xlabel('time')
ylabel('position')
legend('true position','perturbed position','recovered position')
print('TrueStatePerturbedKalmanPosition','-djpeg')

figure
plot(time,TrueState(2,:),'Linewidth',3,'Color','magenta')
hold on
plot(time,mu(2,:),'o','MarkerSize',5,'Color','green')
title('True State vs. Recovered')
xlabel('time')
ylabel('velocity')
legend('true velocity','recovered velocity')
print('TrueStateKalmanVelocity','-djpeg')

figure
plot(time,KalmanGain(1,:),'*','MarkerSize',5,'Color','red')
hold on
plot(time,KalmanGain(2,:),'o','MarkerSize',5,'Color','green')
title('Kalman Gain')
xlabel('time')
legend('position gain','velocity gain')
print('KalmanGain','-djpeg')

figure
plot(time,MSE,'*','MarkerSize',5,'Color','red')
hold on
plot(time,trP,'o','MarkerSize',5,'Color','green')
title('Errors')
xlabel('time')
legend('mean square error','trace of analysis covariance')
print('Errors','-djpeg')

toc()
end

