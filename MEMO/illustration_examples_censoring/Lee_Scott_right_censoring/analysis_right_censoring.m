clear all;
close all;
clc;

%% Load analysis parameters and set options
load dataset_righ_censoring

visualize = true;
% visualize = false;

n_replicates = 30;
thinning = 3;

mu = 1;
sigma = sigma(1:thinning:end);
cen_time = cen_time(1:thinning:end);
n_replicates = round(n_replicates/thinning);

%% Initialization
counter = 0;

N_correct_AIC_censored = zeros(length(sigma),length(cen_time));
N_correct_BIC_censored = zeros(length(sigma),length(cen_time));

N_correct_AIC_truncated = zeros(length(sigma),length(cen_time));
N_correct_BIC_truncated = zeros(length(sigma),length(cen_time));

%% Loop
for i1 = 1:length(sigma) % sigmas
    for i2 = 1:n_replicates % relization
        for i3 = 1:length(cen_time) % censoring time

% Counter update
counter = counter + 1;

% Progress report
disp([num2str(100*counter/(length(sigma)*n_replicates*length(cen_time))) '%']);

%% Generation of data (without censoring and truncation)
paramt.pp = 1;
paramt.mu = mu;
paramt.C = sigma(i1)^2;

y = util_gmm(paramt.pp, paramt.mu, paramt.C, sample_size);


%% Loop: Number of modes
N = 5;
AIC_truncated = nan(1,N);
BIC_truncated = nan(1,N);

for n = 1:5
    
%% Data processing - censoring
% Scenario: Censoring
c_lb = -5;
c_ub = log(cen_time(i3));
t_lb = -50;
t_ub = 50;

% % Scenario: Truncation
% c_lb = -5;
% c_ub = 50;
% t_lb = -50;
% t_ub = 1;


idxtu = y > t_ub;    % truncation
idxtl = y < t_lb;    % truncation
idxt = any(idxtu | idxtl,2);

idxcu = y > c_ub;    % censoring
idxcl = y < c_lb;    % censoring
idxc = any(idxcu | idxcl,2);

x = y;
x(idxcu) = c_ub;
x(idxcl) = c_lb;
x = x(~idxt,:);

pattern = em_censor_pattern(x,[t_lb;c_ub]);
xc = x(pattern.complete,:);

%% Fitting of mixture model - censoring
% Initialize by k-means 
%   run k-means five times and initialize each algorithm with the
%   parameters that achieves the highest likelihood
par_k = init_kmeans(x, n, 5);

% truncated and censored data EM
par_k_c_ub_ll = zeros(5,1);
for m = 1:5
    [temp,par_k_c_ub_ll(m)] = em_tc1_post(x, par_k{m}.pp, par_k{m}.mu, par_k{m}.C, pattern, [t_lb;t_ub]);
end
[temp,par_k_c_ub_idx] = max(par_k_c_ub_ll);
param_init = par_k{par_k_c_ub_idx};
param = em_tc1(x,n,[t_lb;t_ub],[c_lb;c_ub],param_init);

[temp,idx] = sort(param.mu(:,1), 'ascend');
param.pp = param.pp(idx);
param.mu = param.mu(idx,:);
param.C = param.C(:,:,idx);

AIC_censored(n) = param.AIC;
BIC_censored(n) = param.BIC;

%% Visualization
if visualize
    figure(2*n-1); hold off;

    % grid
    [xg] = paramt.mu + sqrt(paramt.C)*linspace(-5,5,1000);

    % data points
    plot(x(:,1),0,'xr','markersize', 5); hold on

    plot(paramt.mu(:,1),0, '+k', ...
         param.mu(:,1),0, 'ok', 'linewidth', 2)

    pdf_sum = zeros(size(xg));
    for k = 1:size(paramt.mu,1)
        pdf_sum = pdf_sum + paramt.pp(k)*pdf('normal',xg,paramt.mu(k,:),sqrt(paramt.C(:,:,k)));
        plot(xg,paramt.pp(k)*pdf('normal',xg,paramt.mu(k,:),sqrt(paramt.C(:,:,k))), '--b', 'linewidth', 2);
    end
    plot(xg,pdf_sum, '-b', 'linewidth', 2);

    pdf_sum = zeros(size(xg));
    for k=1:size(param.mu,1)
        pdf_sum = pdf_sum + param.pp(k)*pdf('normal',xg,param.mu(k,:),sqrt(param.C(:,:,k)));
        plot(xg,param.pp(k)*pdf('normal',xg,param.mu(k,:),sqrt(param.C(:,:,k))), '--k', 'linewidth', 2);
    end
    plot(xg,pdf_sum, '-k', 'linewidth', 2);

    title('censored EM', 'FontSize', 15)
    xlim(xg([1,end])), axis square
    
    drawnow
end

%% Data processing - truncation
% Scenario: Truncation
c_lb = -5;
c_ub = 50;
t_lb = -50;
t_ub = log(cen_time(i3));

idxtu = y > t_ub;    % truncation
idxtl = y < t_lb;    % truncation
idxt = any(idxtu | idxtl,2);

idxcu = y > c_ub;    % censoring
idxcl = y < c_lb;    % censoring
idxc = any(idxcu | idxcl,2);

x = y;
x(idxcu) = c_ub;
x(idxcl) = c_lb;
x = x(~idxt,:);

pattern = em_censor_pattern(x,[t_lb;c_ub]);
xc = x(pattern.complete,:);

%% Fitting of mixture model - truncation
% Initialize by k-means 
%   run k-means five times and initialize each algorithm with the
%   parameters that achieves the highest likelihood
par_k = init_kmeans(x, n, 5);

% truncated and censored data EM
par_k_c_ub_ll = zeros(5,1);
for m = 1:5
    [temp,par_k_c_ub_ll(m)] = em_tc1_post(x, par_k{m}.pp, par_k{m}.mu, par_k{m}.C, pattern, [t_lb;t_ub]);
end
[temp,par_k_c_ub_idx] = max(par_k_c_ub_ll);
param_init = par_k{par_k_c_ub_idx};
param = em_tc1(x,n,[t_lb;t_ub],[c_lb;c_ub],param_init);

[temp,idx] = sort(param.mu(:,1), 'ascend');
param.pp = param.pp(idx);
param.mu = param.mu(idx,:);
param.C = param.C(:,:,idx);

AIC_truncated(n) = param.AIC;
BIC_truncated(n) = param.BIC;

%% Visualization - truncation
if visualize
    figure(2*n); hold off;

    % grid
    [xg] = paramt.mu + sqrt(paramt.C)*linspace(-5,5,1000);

    % data points
    plot(x(:,1),0,'xr','markersize', 5); hold on

    plot(paramt.mu(:,1),0, '+k', ...
         param.mu(:,1),0, 'ok', 'linewidth', 2)

    pdf_sum = zeros(size(xg));
    for k = 1:size(paramt.mu,1)
        pdf_sum = pdf_sum + paramt.pp(k)*pdf('normal',xg,paramt.mu(k,:),sqrt(paramt.C(:,:,k)));
        plot(xg,paramt.pp(k)*pdf('normal',xg,paramt.mu(k,:),sqrt(paramt.C(:,:,k))), '--b', 'linewidth', 2);
    end
    plot(xg,pdf_sum, '-b', 'linewidth', 2);

    pdf_sum = zeros(size(xg));
    for k=1:size(param.mu,1)
        pdf_sum = pdf_sum + param.pp(k)*pdf('normal',xg,param.mu(k,:),sqrt(param.C(:,:,k)));
        plot(xg,param.pp(k)*pdf('normal',xg,param.mu(k,:),sqrt(param.C(:,:,k))), '--k', 'linewidth', 2);
    end
    plot(xg,pdf_sum, '-k', 'linewidth', 2);

    title('truncated EM', 'FontSize', 15)
    xlim(xg([1,end])), axis square
    
    drawnow
end

%% end of loop over number of mixture components
end

%% Selected model index
N_correct_AIC_censored(i1,i3) = N_correct_AIC_censored(i1,i3) + (AIC_censored(1) <= min(AIC_censored(2:end)));
N_correct_BIC_censored(i1,i3) = N_correct_BIC_censored(i1,i3) + (BIC_censored(1) <= min(BIC_censored(2:end)));

N_correct_AIC_truncated(i1,i3) = N_correct_AIC_truncated(i1,i3) + (AIC_truncated(1) <= min(AIC_truncated(2:end)));
N_correct_BIC_truncated(i1,i3) = N_correct_BIC_truncated(i1,i3) + (BIC_truncated(1) <= min(BIC_truncated(2:end)));

        end
    end
end

%% Normalization
N_correct_AIC_censored = 100*N_correct_AIC_censored/n_replicates;
N_correct_BIC_censored = 100*N_correct_BIC_censored/n_replicates;

N_correct_AIC_truncated = 100*N_correct_AIC_truncated/n_replicates;
N_correct_BIC_truncated = 100*N_correct_BIC_truncated/n_replicates;

%% Visualization
figure;

subplot(2,2,1);
colormap('gray')
imagesc(N_correct_AIC_censored(end:-1:1,:),[0,100]);
set(gca,'xtick',1:length(cen_time),'xticklabel',round(10*cen_time)/10,...
        'ytick',1:length(sigma),'yticklabel',sigma(end:-1:1))
title('AIC - censored');
xlabel('censoring time');
ylabel('sigma');
h = colorbar;
ylabel(h,'selection accuracy')

subplot(2,2,2);
colormap('gray')
imagesc(N_correct_BIC_censored(end:-1:1,:),[0,100]);
set(gca,'xtick',1:length(cen_time),'xticklabel',round(10*cen_time)/10,...
        'ytick',1:length(sigma),'yticklabel',sigma(end:-1:1))
title('BIC - censored');
xlabel('censoring time');
ylabel('sigma');
h = colorbar;
ylabel(h,'selection accuracy')

subplot(2,2,3);
colormap('gray')
imagesc(N_correct_AIC_truncated(end:-1:1,:),[0,100]);
set(gca,'xtick',1:length(cen_time),'xticklabel',round(10*cen_time)/10,...
        'ytick',1:length(sigma),'yticklabel',sigma(end:-1:1))
title('AIC - truncated');
xlabel('censoring time');
ylabel('sigma');
h = colorbar;
ylabel(h,'selection accuracy')

subplot(2,2,4);
colormap('gray')
imagesc(N_correct_BIC_truncated(end:-1:1,:),[0,100]);
set(gca,'xtick',1:length(cen_time),'xticklabel',round(10*cen_time)/10,...
        'ytick',1:length(sigma),'yticklabel',sigma(end:-1:1))
title('BIC - truncated');
xlabel('censoring time');
ylabel('sigma');
h = colorbar;
ylabel(h,'selection accuracy')

print('preliminary_result','-depsc2')

save all
