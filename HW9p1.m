tic()
clc
clear
close all
set(0,'DefaultFigureVisible','off')

%% preliminaries
colors
Ne = 1000;
size_k = 100;
R = 1;
w = zeros(1,Ne);
X_init = 10.*randn(1,Ne);
name = 'SPF_HW9';
Ens = [X_init;zeros(size_k,Ne)];
True_X = [10.*randn,zeros(1,size_k)];
Obs = zeros(1,size_k+1);
SPF = [mean(X_init,2),zeros(1,size_k)];
%% 

%% standard particle filter
f = @(x,k)(0.5.*x + 25.*x./(1+x.^2)+8*cos(1.2*k)+10*randn(1,size(x,2)));   % stochastic model
h = @(x)(0.05.*x.^2 + randn(1,size(x,2)));                                 % observation function

for k=1:size_k
    True_X(k+1) = f(True_X(k),k);
    Obs(k+1) = h(True_X(k+1));
    Ens(k+1,:) = f(Ens(k,:),k);
    Fake_Obs = h(Ens(k+1,:));
    for ii=1:Ne
       w(ii) = 0.5.*((Obs(k+1)-Fake_Obs(ii))^2)/R; 
    end 
    W = normalizeweights(w);
    Ens(k+1,:) = resamplingmmo(W,Ens(k+1,:),Ne,1);
    SPF(:,k+1) = mean(Ens(k+1,:),2);
    
%     figure
%     histogram(Ens(k+1,:),40,'Normalization','pdf');
%     axis([-36 36 0 0.6])
%     xlabel('x')
%     ylabel('p(x)')
%     title(['ensemble distribution at k = ',num2str(k)])
%     print(['p(x)_k_',num2str(k)],'-djpeg')
end

%%

%% plots

figure
h1 = plot(0:size_k,True_X,'LineWidth',2,'Color',Color(:,35));
hold on
h2 = plot(0:size_k,SPF,'LineWidth',2,'Color',Color(:,12));
title('true vs. SPF')
xlabel('step')
ylabel('x')
legend([h1(1),h2(1)],'true','SPF')
print('SPF','-djpeg')
hold off

Error = True_X - SPF;
Error = abs(Error);
average_RMSE = mean(Error,2)

figure
plot(0:size_k,Error,'LineWidth',2,'Color',Color(:,16))
title('RMSE')
xlabel('step')
ylabel('error')
print('RMSE','-djpeg')


toc()








