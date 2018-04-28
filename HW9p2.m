tic()
clc
clf
clear
close all
munlock('UWerr_v6')
munlock('UWerr_fft')
clc
clf
clear
close all
% set(0,'DefaultFigureVisible','on')

%% preliminaries
% n = [10, 50 100, 200, 500];
n = [10, 20];
n_sz = length(n);
accept_rate = zeros(n_sz,1);
tau_avg = zeros(n_sz,1);
Ne = 10^5;
%%


%% Markov Chains

for ii=1:n_sz
   x = zeros(n(ii),Ne);
   x(:,1) = randn(n(ii),1);
   sigma = 1/sqrt(n(ii));
   accept = 0;
   for jj=2:Ne
       x_prop = x(:,jj-1) + sigma.*randn(n(ii),1);
       alpha = min(1,exp(-0.5*(x_prop'*x_prop-x(:,jj-1)'*x(:,jj-1))));
       if rand<alpha
           x(:,jj) = x_prop;
           accept = accept + 1;
       else
           x(:,jj) = x(:,jj-1);
       end
   end
   accept_rate(ii) = accept/Ne;
   [~,~,~,tau_avg(ii),~,~] = UWerr_fft(x');
end
%%

%% plots
colors

figure(11)
plot(n,accept_rate,'LineWidth',3,'Color',Color(:,10))
hold on
plot(n,accept_rate,'.','MarkerSize',13,'Color',Color(:,10))
xlabel('dimension','FontSize',20)
ylabel('acceptance rate','FontSize',20)
hold off

figure(12)
plot(n,tau_avg,'LineWidth',3,'Color',Color(:,18))
hold on
plot(n,tau_avg,'.','MarkerSize',13,'Color',Color(:,18))
xlabel('dimension','FontSize',20)
ylabel('tau','FontSize',20)
hold off

%%

accept_rate
tau_avg

toc()