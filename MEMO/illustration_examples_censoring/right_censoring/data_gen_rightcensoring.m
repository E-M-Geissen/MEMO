
clc;
clear all;
close all;


rng('shuffle')


% % Options for dataset
dx = [0]; % continous data
N = 200; % number of data points drawn from log-normal distribution
Q = 20;  % number of realization = number of times data drawn from log-normal distribution


% parameter \mu and \sigma of log-normal distribution
mus=1;  % true mu for data generation
sigma(1)=0.2; % first sigma in grid
mean=exp(mus+((sigma(1)^2)/2)); % calculate distribution mean for initial parameters


%% Calculate more log-normal sigmas calculated defined by coresponding to means defined by given censoring values
multiples=linspace(1,5,10); % multiples of mean used for censoring value
cen_time=multiples*mean;


% calculate sigmas that correspond to the means determined by censoring
% value = multiples of mean(mus, sigma(1))
for i=1:length(cen_time)
    sigma_hilf(i)= sqrt(2*(log(cen_time(i))-mus));
end
% use only every second  sigma determined bay mean=cen_time and calculate corresponding log-normal stds
b=0;
for i=1:2:length(cen_time)
    b=b+1;
    sigma(i)= sigma_hilf(b);
    std(i)=sqrt((exp(sigma(i)^2)-1)*exp(2*mus+sigma(i)^2));
end

% for equal spacing on y-axis use sigmas corresponding to log-normal stds in between the ones determined by mean=cen_time
for i=2:2:(length(cen_time)-1)
    sigma(i)=sqrt(log(0.5+sqrt( ((std(i-1)+0.5*(std(i+1)-std(i-1)))^2)/exp(2)+0.25)));
end

% %% Options for model
% o_max = 5; % maximal number of models (= maximal number of subpopulations)
%
% % lower and upper limits for estimations
% mu_min = -1;
% mu_max = +3;
%
% sigma_min = 1e-3;
% sigma_max = 1e+3;
%
%
% % Options for multi-start optimization
% options_MS.n_starts = 30;
% options_MS.plot = plot_opt;
% options_MS.proposal = 'latin hypercube';
% options_MS.fmincon = optimset('algorithm','interior-point',...%'active-set',...
%     'display','off',...
%     'GradObj','on',...
%     'MaxIter',1000,...
%     'MaxFunEvals',3000);

%% DATA GENERATION
sample_fun = @(N,sigma) exp(sigma*randn(N,1)+mus); % mu=1; sigma from grid



D=zeros(length(sigma), Q,length(cen_time),2,N);
% Loop: sigmas

for p = 1:length(sigma)
    
    % Loop: realizations
    for i = 1:Q
        
        % Generation of artificial data, data will be evaluated for every censoring time
        % and model complexity
        
        x = sample_fun(N,sigma(p));
        
        % Interval censoring?
        if dx > 0
            xm_raw = dx*ceil(x/dx);
        else
            xm_raw = x;
        end
        
        % Loop: Censoring values
        for j=1:length(cen_time)
            
            %% 3 different approaches
            %% V1: censored data is refused (=truncation)
            xm=xm_raw(xm_raw<=cen_time(j));
            
            D(p,i,j,1,:)=xm;
            %% V2: censored data is set to censoring time
            xm=xm_raw;
            xm(xm>cen_time(j))=cen_time(j);
            
             D(p,i,j,2,:)=xm;
            
        end
        
    end
    
end

% %% calculate percent fraction of realizations correctly identfying the single source distribution correctly
%
% [~,selection_AIC_v1]  = min(AIC_v1,[],1);    selection_AIC_v1   = squeeze(selection_AIC_v1);
% [~,selection_AICc_v1] = min(AICc_v1,[],1);   selection_AICc_v1  = squeeze(selection_AICc_v1);
% [~,selection_BIC_v1]  = min(BIC_v1,[],1);    selection_BIC_v1   = squeeze(selection_BIC_v1);
%
% for p=1:length(cen_time)
%     selection_AIC_tr_v1(:,:,p)   = transpose(selection_AIC_v1(:,:,p));
%     selection_AICc_tr_v1(:,:,p)  = transpose(selection_AICc_v1(:,:,p));
%     selection_BIC_tr_v1(:,:,p)   = transpose(selection_BIC_v1(:,:,p));
% end
%
% for p=1:length(cen_time)
%     for j = 1:length(sigma)
%         for o = 1:o_max
%             perc_selection_AIC_v1(o,j,p)   = (sum(selection_AIC_tr_v1(j,:,p)   == o))/Q;
%             perc_selection_AICc_v1(o,j,p)  = (sum(selection_AICc_tr_v1(j,:,p)  == o))/Q;
%             perc_selection_BIC_v1(o,j,p)   = (sum(selection_BIC_tr_v1(j,:,p)   == o))/Q;
%         end
%     end
% end
%
% %%
% [~,selection_AIC_v2]  = min(AIC_v2,[],1);    selection_AIC_v2   = squeeze(selection_AIC_v2);
% [~,selection_AICc_v2] = min(AICc_v2,[],1);   selection_AICc_v2  = squeeze(selection_AICc_v2);
% [~,selection_BIC_v2]  = min(BIC_v2,[],1);    selection_BIC_v2   = squeeze(selection_BIC_v2);
%
% for p=1:length(cen_time)
%     selection_AIC_tr_v2(:,:,p)   = transpose(selection_AIC_v2(:,:,p));
%     selection_AICc_tr_v2(:,:,p)  = transpose(selection_AICc_v2(:,:,p));
%     selection_BIC_tr_v2(:,:,p)   = transpose(selection_BIC_v2(:,:,p));
% end
%
% for p=1:length(cen_time)
%     for j = 1:length(sigma)
%         for o = 1:o_max
%             perc_selection_AIC_v2(o,j,p)   = (sum(selection_AIC_tr_v2(j,:,p)   == o))/Q;
%             perc_selection_AICc_v2(o,j,p)  = (sum(selection_AICc_tr_v2(j,:,p)  == o))/Q;
%             perc_selection_BIC_v2(o,j,p)   = (sum(selection_BIC_tr_v2(j,:,p)   == o))/Q;
%         end
%     end
% end
%
% %%
% [~,selection_AIC_v3]  = min(AIC_v3,[],1);    selection_AIC_v3   = squeeze(selection_AIC_v3);
% [~,selection_AICc_v3] = min(AICc_v3,[],1);   selection_AICc_v3  = squeeze(selection_AICc_v3);
% [~,selection_BIC_v3]  = min(BIC_v3,[],1);    selection_BIC_v3   = squeeze(selection_BIC_v3);
%
% for p=1:length(cen_time)
%     selection_AIC_tr_v3(:,:,p)   = transpose(selection_AIC_v3(:,:,p));
%     selection_AICc_tr_v3(:,:,p)  = transpose(selection_AICc_v3(:,:,p));
%     selection_BIC_tr_v3(:,:,p)   = transpose(selection_BIC_v3(:,:,p));
% end
%
% for p=1:length(cen_time)
%     for j = 1:length(sigma)
%         for o = 1:o_max
%             perc_selection_AIC_v3(o,j,p)   = (sum(selection_AIC_tr_v3(j,:,p)   == o))/Q;
%             perc_selection_AICc_v3(o,j,p)  = (sum(selection_AICc_tr_v3(j,:,p)  == o))/Q;
%             perc_selection_BIC_v3(o,j,p)   = (sum(selection_BIC_tr_v3(j,:,p)   == o))/Q;
%         end
%     end
% end
%
%
% %% SAVE RESULTS
%
% save('right_cen_data_BIC','perc_selection_BIC_v1','perc_selection_BIC_v2','perc_selection_BIC_v3', 'mus','sigma', 'dx','o_max' ,'cen_time');
% save('right_cen_data_AIC','perc_selection_AIC_v1','perc_selection_AIC_v2','perc_selection_AIC_v3', 'mus','sigma', 'dx','o_max' ,'cen_time');
% save('right_cen_data_AICc','perc_selection_AICc_v1','perc_selection_AICc_v2','perc_selection_AICc_v3', 'mus','sigma', 'dx','o_max' ,'cen_time');


% %% VISUALIZATION OF RESULTS
%
% % set figure size
% width=5;
% height=5;
%
% TextSizes.DefaultAxesFontSize = 6;
% TextSizes.DefaultTextFontSize =6;
% set(0,TextSizes);
%
%
% % Legend constcurtion
% for o = 1:o_max
%     if o == 1
%         ls{o} = [num2str(o,'%d') ' component'];
%     else
%         ls{o} = [num2str(o,'%d') ' components'];
%     end
% end
%
% % BIC:
%
% for i = 1:3
%
%     switch i
%         case 1, perc = perc_selection_BIC_v1;  ts = 'censored data refused';name='v1'
%         case 2, perc = perc_selection_BIC_v2; ts = 'censored data = censoring time'; name='v2';
%         case 3, perc = perc_selection_BIC_v3; ts = 'censoring time considered';  name='v3';
%     end
%
%     hFig(i)=figure;
%     set(hFig(i), 'Units','centimeters','Position', [10 10 width height])
%
%     t=1;
%     for p=1:length(sigma)
%         for j=1:length(cen_time)
%             xx(t)=cen_time(j);
%             y(t)=sqrt((exp(sigma(p)^2)-1)*exp(2*mus+sigma(p)^2));
%             z(t)=perc(1,p,j)*100;
%             t=t+1;
%         end
%     end
%
%
%     s=zeros(size(xx));
%     s(:,:)=20;
%
%     colormap gray
%
%     scatter(xx, y,s,z,'fill');
%     hold on
%
%     h=scatter(xx, y,s,[0,0,0],'LineWidth',1);
%     hold all
%     plot(cen_time(1:1:5),sqrt((exp(sigma(1:2:9).^2)-1).*exp(2*mus+sigma(1:2:9).^2)),'-','Color',[0 0.749019622802734 0.749019622802734],'LineWidth',1);
%     set(gca, 'LooseInset', [0,0,0,0]);
%     hold all
%     ylim([-0.5 21.5])
%     xlim([2.2 14.5])
%     axis square
%
%     ylabel('standard deviation');
%     xlabel('censoring time');
%     caxis([0 100]);
%
%
%     hc = colorbar;
%     ylabel(hc,'% identified correctly')
%
%     title(ts);
%     set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
%
%
%     print('-depsc2','-r1000',['./figs/right_cen_',name]);
%
%
% end
%
%
% %% illustrate effact on mean
%
% % set figure size
% width=5
% height=5
%
% TextSizes.DefaultAxesFontSize = 6;
% TextSizes.DefaultTextFontSize = 6;
% set(0,TextSizes);
%
% load mycmap
% clear x
% clear mean
%
% for i=1:3
%     switch i
%
%         case 1
%
%             means_vx=means_v1;ts = 'censored data refused';name='v1';
%
%         case 2
%             means_vx=means_v2;  ts='censored data = censoring time';name='v2';
%
%         case 3
%             means_vx=means_v3;  ts = 'censoring considered';name='v3';
%     end
%
%
%
%
%
%     for j=1:length(cen_time)
%         for p=1:length(sigma)
%
%             mean_m_vx(p,j)= mean(means_vx(:,p,j));
%             mean_nom=exp(1+(sigma(p).^2)./2)
%             rel(p,j)= (mean_m_vx(p,j)-mean_nom)/mean_nom*100;
%         end
%     end
%       hFig(i)=figure;
%     set(hFig(i), 'Units','centimeters','Position', [10 10 width height])
%
%
%
%
%
%
%     t=1;
%     for p=1:length(sigma)
%         for j=1:length(cen_time)
%             xvalue(t)=cen_time(j);
%             yvalue(t)=sqrt((exp(sigma(p)^2)-1)*exp(2*mus+sigma(p)^2));
%             z(t)=rel(p,j);
%             t=t+1;
%         end
%     end
%
%
%
%     s=zeros(size(xvalue));
%     s(:,:)=20;
%
%     set(gcf,'Colormap',mycmap)
%
%     scatter(xvalue, yvalue,s,z,'fill');
%     hold on
%     scatter(xvalue, yvalue,s,z,'fill');
%     hold on
%     h=scatter(xvalue, yvalue ,s,[0,0,0],'LineWidth',1);
%     hold all
%      ylim([2 21])
%      xlim([2.3 14])
%     axis square
%     ylabel('standard deviation');
%     xlabel('censoring time');
%     caxis([-100 100]);
%     h = colorbar;
%     ylabel(h,'% mean deviation from true mean')
%
%     title(ts)
%     set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
%
%
%     print('-depsc2','-r1000',['./right_cen_mean',name]);
% end
%
