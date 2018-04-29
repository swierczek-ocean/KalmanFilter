tic()
%% preliminaries
colors
load Data
HareLynx = Data(:,2:3)';
dim = 6;
Ne = 5e4;
theta = zeros(dim,Ne);
theta(:,1) = [0.5861; 0.2345; 0.7780; 0.1768; 2.5786; 3.8248];
num_sims = 100;
dt = 0.1;
time_steps = 10/dt;
D = [0.1,0.1,0.1,0.1,0.2,0.2];
PropCov = diag(D);
sigma = 0.215;
%% 

%% Markov Chain
accept = 0;
D_current = run_LV(theta(:,1),dt);
size(D_current)
F_current = sum(sum((D_current-HareLynx).^2));

for jj=2:Ne
    theta_prop = theta(:,jj-1) + sigma.*PropCov*randn(dim,1);
    
    if max(theta_prop)>10 || min(theta_prop)<0
        alpha = 0;
    else
        D_prop = run_LV(theta_prop,dt);
        F_prop = sum(sum((D_prop-HareLynx).^2));
        alpha = exp(-0.5*(F_prop-F_current));
    end
    
    if rand<alpha
        theta(:,jj) = theta_prop;
        D_current = D_prop;
        F_current = F_prop;
        accept = accept + 1;
        fprintf('Sample %g / %g\r',jj,Ne)
    else
        theta(:,jj) = theta(:,jj-1);
    end
end
accept_rate = accept/Ne


%% plots
color1 = 31;
color2 = 29;
color3 = 11;
color4 = 9;
Index = randperm(Ne-100,num_sims)+100;
Temp = run_LV(theta(:,Index(1)),dt);
figure()
plot(0:10,Temp(1,:),'LineWidth',1,'Color',Color(:,color1))
hold on
plot(0:10,Temp(2,:),'LineWidth',1,'Color',Color(:,color2))
xlabel('time in years','FontSize',20)
ylabel('populations','FontSize',20)
title('RWM MCMC simulations')
for ll=2:num_sims
    Temp = run_LV(theta(:,Index(ll)),dt);
    plot(0:10,Temp(1,:),'LineWidth',1,'Color',Color(:,color1))
    plot(0:10,Temp(2,:),'LineWidth',1,'Color',Color(:,color2))
end
h1 = plot(0:10,HareLynx(1,:),'LineWidth',3,'Color',Color(:,color3));
h2 = plot(0:10,HareLynx(2,:),'LineWidth',3,'Color',Color(:,color4));
legend([h1(1),h2(1)],'Hare','Lynx')
hold off
print('LV_sim','-djpeg')
%%

figure()
TrianglePlot(theta(:,3000:end),1)

[~,~,~,tau_avg,~,~] = UWerr_fft(theta')


toc()