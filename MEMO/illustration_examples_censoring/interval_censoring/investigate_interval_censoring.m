clc;
clear all;
close all;
rng('shuffle')


%% OPTIONS
plot_opt ='false';%'true' ;%'false'; % MS plot is shown
plot_opt_2 ='false'; % 'true';% 'false'; % Fit is shown


% Options for dataset
% Large analysis
dx = [0,10.^linspace(-3,0,10)];
N = 100; % number of data points
Q = 10;  % number of realization

ratio=0.2;
mu=1;
sigma = sqrt(log(1/2 + sqrt(1/4 + (ratio.*dx).^2./exp(2*mu))));

% Options for model
o_max = 5; % maximal number of models (distributions)
mu_min = -1;
mu_max = +3;
sigma_min = 1e-6;
sigma_max = 1e+3;



% Options for multi-start optimization
options_MS.n_starts = 30;
options_MS.plot = plot_opt;
options_MS.proposal = 'latin hypercube';
options_MS.fmincon = optimset('algorithm','interior-point',...%'active-set',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',1000,...
    'MaxFunEvals',3000);

%% DATA GENERATION

sample_fun = @(N,sigma) exp(sigma*randn(N,1)+mu); % mu=1; sigma=0.2


%% MODELS
for o = 1:o_max
    % Generation of symbolic variables
    if o >= 2
        str_w = '[';
        str_mu = '[';
        str_log10_sigma = '[';
        for j = 1:o
            if j >= 2
                str_w = [str_w ';'];
                str_mu = [str_mu ';'];
                str_log10_sigma = [str_log10_sigma ';'];
            end
            str_w = [str_w 'w_' num2str(j,'%d')];
            str_mu = [str_mu 'mu_' num2str(j,'%d')];
            str_log10_sigma = [str_log10_sigma 'log10_sigma_' num2str(j,'%d')];
        end
        str_w = [str_w ']'];
        str_mu = [str_mu ']'];
        str_log10_sigma = [str_log10_sigma ']'];
        w = sym(str_w);
        w(1) = 1-sum(w(2:end));
        mu = sym(str_mu);
        log10_sigma = sym(str_log10_sigma);
    else
        w = 1;
        mu = sym('mu');
        log10_sigma = sym('log10_sigma');
    end
    % Model generation
    M{o}.mixture.type = data_opt;
    %    M{o}.mixture.type = 'log-normal';
    M{o}.label.x = 'time [min]';
    M{o}.label.y = 'probability density';
    %
    M{o}.experiment(1).name  = ['order ' num2str(o,'%d')];
    for j = 1:o
        M{o}.experiment(1).w{j}     = w(j);
        M{o}.experiment(1).mu{j}    = mu(j);
        M{o}.experiment(1).sigma{j} = 10^log10_sigma(j);
    end
    M{o}.experiment(1).size  = o;
    % Construct parameter struct
    for j = 1:o
        if j == 1
            parameters{o}.sym = [mu(1);log10_sigma(1)];
            parameters{o}.name = {['\mu_{' num2str(j,'%d') '}'];...
                ['log_{10}(\sigma_{' num2str(j,'%d') '})']};
            parameters{o}.min = [mu_min;log10(sigma_min)];
            parameters{o}.max = [mu_max;log10(sigma_max)];
        else
            parameters{o}.sym = [parameters{o}.sym;mu(j);log10_sigma(j);w(j)];
            parameters{o}.name(end+1:end+3) = {['\mu_{' num2str(j,'%d') '}'];...
                ['log_{10}(\sigma_{' num2str(j,'%d') '})'];...
                ['w_{' num2str(j,'%d') '}']};
            parameters{o}.min = [parameters{o}.min;mu_min;log10(sigma_min);0];
            parameters{o}.max = [parameters{o}.max;mu_max;log10(sigma_max);1];
        end
    end
    parameters{o}.number = length(parameters{o}.sym);
    % Compile model
    [M{o},parameters{o}.constraints] = getMixtureModel(M{o},parameters{o}.sym);
end


%% EVALUATION OF FITTING STATISTICS
AIC = zeros(o_max,Q,length(dx),length(sigma));
AICc = zeros(o_max,Q,length(dx),length(sigma));
BIC = zeros(o_max,Q,length(dx),length(sigma));
cAIC = zeros(o_max,Q,length(dx),length(sigma));
cAICc = zeros(o_max,Q,length(dx),length(sigma));
cBIC = zeros(o_max,Q,length(dx),length(sigma));

% Loop: realizations

for p=2:length(sigma)
    for i = 1:Q
        
        % Generation of artificial data (equal data for every model and every dx)
        x = sample_fun(N,sigma(p));
        
        for j = 1:length(dx)
            if dx(j) > 0
                
                offset=unifrnd(0,dx(j));
                xm = dx(j)*ceil((x-offset)/dx(j))+offset;
            else
                xm = x;
                length(unique(xm))
            end
            
            D{1}.name = 'test data';
            D{1}.description = [];
            D{1}.data.uncensored = xm;
            D{1}.data.censored = [];
            
            
            
            % Fitting
            for o = 1:o_max
                
                
                if o==1
                    parameters{o}.guess=[1;log10(sigma(p))];
                    
                    
                elseif o==2
                    
                    parameters{o}.guess = [1;log10(sigma(p));1;log10(sigma(p));0.6];
                    
                elseif o==3
                    
                    parameters{o}.guess = [1;log10(sigma(p));1;log10(sigma(p));0.3;1;log10(sigma(p));0.2];
                    
                elseif o==4
                    
                    parameters{o}.guess = [1;log10(sigma(p));1;log10(sigma(p));0.3;1;log10(sigma(p));0.2;1;log10(sigma(p));0.1];
                    
                elseif o==5
                    
                    parameters{o}.guess = [1;log10(sigma(p));1;log10(sigma(p));0.35;1;log10(sigma(p));0.3;1;log10(sigma(p));0.2;1;log10(sigma(p));0.1];
                    
                    
                    
                end
                % Process data for faster computations
                D{1}.observation_interval = dx(j);
                D = processData(D);
                k = parameters{o}.number;
                
                % Run estimation
                [parameters{o},M{o}.fh.fit] = optimizeMultiStart(parameters{o},@(theta,opt) logLikelihoodMM(theta,M{o},D,opt),options_MS);
                % Evaluation of model selection criterion
                cparameter_est{o,i,p,j} = parameters{o};
                
                if o==1
                    mui=parameters{o}.MS.MAP.par(1,1);
                    sigi=10^parameters{o}.MS.MAP.par(2,1);
                    cmeans(j,i,p)=exp(mui+((sigi^2)/2));
                    cstds(j,i,p)=sqrt((exp(sigi^2)-1)*exp(2*mui+sigi^2));
                end
                
                cAIC(o,i,j,p)  = -2*parameters{o}.MS.MAP.logPost + 2*k;
                cAICc(o,i,j,p) = -2*parameters{o}.MS.MAP.logPost + 2*N*(k+1)/(N-k-1);
                cBIC(o,i,j,p)  = -2*parameters{o}.MS.MAP.logPost + k*log(N);
                
                % Plot
                if strcmp(plot_opt_2,'true')
                    
                    plotMixtureModel4(parameters{o},M{o},[],D);
                    title(['censoring considered \sigma=',num2str(sigma(p)),' dx=',num2str(dx(j))]);
                    options.plot_type = 'pdf';
                    plotMixtureModel4(parameters{o},M{o},[],D,options);
                    title(['censoring considered \sigma=',num2str(sigma(p)),' dx=',num2str(dx(j))]);
                end
                
                
                % Process data for faster computations
                D{1}.observation_interval = 0;
                D = processData(D);
                k = parameters{o}.number;
                
                % Run estimation
                [parameters{o},M{o}.fh.fit] = optimizeMultiStart(parameters{o},@(theta,opt) logLikelihoodMM(theta,M{o},D,opt),options_MS);
                % Evaluation of model selection criterion
                parameter_est{o,i,p,j} = parameters{o};
                
                if o==1
                    mui=parameters{o}.MS.MAP.par(1,1);
                    sigi=10^parameters{o}.MS.MAP.par(2,1);
                    means(j,i,p)=exp(mui+((sigi^2)/2));
                    stds(j,i,p)=sqrt((exp(sigi^2)-1)*exp(2*mui+sigi^2));
                end
                
                AIC(o,i,j,p)  = -2*parameters{o}.MS.MAP.logPost + 2*k;
                AICc(o,i,j,p) = -2*parameters{o}.MS.MAP.logPost + 2*N*(k+1)/(N-k-1);
                BIC(o,i,j,p)  = -2*parameters{o}.MS.MAP.logPost + k*log(N);
                % Plot
                if strcmp(plot_opt_2,'true')
                    plotMixtureModel4(parameters{o},M{o},[],D);
                    options.plot_type = 'pdf';
                    title(['censoring considered \sigma=',num2str(sigma(p)),' dx=',num2str(dx(j))]);
                    plotMixtureModel4(parameters{o},M{o},[],D,options);
                    title(['censoring considered \sigma=',num2str(sigma(p)),' dx=',num2str(dx(j))]);
                end
            end
            
        end
        save([date '_sigmaNr_' num2str(p,'%d') '_real_=' num2str(i,'%d') ]);
    end
    save([date '_sigmaNr_' num2str(p,'%d')]);
end



%% SAVE RESULTS
save([date '_N=' num2str(N,'%d') '_Q=' num2str(Q,'%d') '_O=' num2str(o_max,'%d')]);



%% EVALUATION OF STATISTICS

[~,selection_AIC]  = min(AIC,[],1);
[~,selection_cAIC]  = min(cAIC,[],1);
[~,selection_AICc] = min(AICc,[],1);
[~,selection_cAICc] = min(cAICc,[],1);
[~,selection_BIC]  = min(BIC,[],1);
[~,selection_cBIC]  = min(cBIC,[],1);

%%
selection_cAIC  = squeeze(selection_cAIC);
selection_cAICc = squeeze(selection_cAICc);
selection_cBIC  = squeeze(selection_cBIC);
selection_AIC   = squeeze(selection_AIC);
selection_AICc  = squeeze(selection_AICc);
selection_BIC   = squeeze(selection_BIC);



%%
for p=1:length(sigma)
    selection_cAIC_tr(:,:,p)  = transpose(selection_cAIC(:,:,p));
    selection_cAICc_tr(:,:,p) = transpose(selection_cAICc(:,:,p));
    selection_cBIC_tr(:,:,p)  = transpose(selection_cBIC(:,:,p));
    selection_AIC_tr(:,:,p)   = transpose(selection_AIC(:,:,p));
    selection_AICc_tr(:,:,p)  = transpose(selection_AICc(:,:,p));
    selection_BIC_tr(:,:,p)   = transpose(selection_BIC(:,:,p));
end

%%
for p=1:length(sigma)
    for j = 1:length(dx)
        for o = 1:o_max
            perc_selection_cAIC(o,j,p)  = sum(selection_cAIC_tr(j,:,p)  == o)/Q;
            perc_selection_cAICc(o,j,p) = sum(selection_cAICc_tr(j,:,p) == o)/Q;
            perc_selection_cBIC(o,j,p)  = sum(selection_cBIC_tr(j,:,p)  == o)/Q;
            perc_selection_AIC(o,j,p)   = sum(selection_AIC_tr(j,:,p)   == o)/Q;
            perc_selection_AICc(o,j,p)  = sum(selection_AICc_tr(j,:,p)  == o)/Q;
            perc_selection_BIC(o,j,p)   = sum(selection_BIC_tr(j,:,p)   == o)/Q;
        end
    end
end

%%
save('grid_data_BIC','perc_selection_BIC','perc_selection_cBIC', 'sigma', 'dx','o_max' );
save('grid_data_AIC','perc_selection_AIC','perc_selection_cAIC', 'sigma', 'dx','o_max' );
save('grid_data_AICc','perc_selection_AICc','perc_selection_cAICc', 'sigma', 'dx','o_max' );
%%
%% set figure size
width=6;
height=5;

%%
TextSizes.DefaultAxesFontSize = 6;
TextSizes.DefaultTextFontSize =6;
set(0,TextSizes);


%%
% Legend construction
for o = 1:o_max
    if o == 1
        ls{o} = [num2str(o,'%d') ' component'];
    else
        ls{o} = [num2str(o,'%d') ' components'];
    end
end

% BIC:
for i = 1:2

    switch i
        case 1, perc = perc_selection_BIC;  ts = 'censoring disregarded';
        case 2, perc = perc_selection_cBIC; ts = 'censoring considered';
    end
    
    hFig(i)=figure;
    set(hFig(i), 'Units','centimeters','Position', [10 10 width height])
    
    t=1;
    for p=1:length(sigma)
        for j=1:length(dx)
            x(t)=log10(dx(j));
            y(t)=log10(sqrt((exp(sigma(p)^2)-1)*exp(2*1+sigma(p)^2)));
            %y(t)=log10(sigma(p));
            z(t)=perc(1,j,p)*100;
            t=t+1;
        end
    end
    
    x(x==-inf)=-4;
    
    s=zeros(size(x));
    s(:,:)=20;
    
    
    colormap gray
    
    scatter(x, y,s,z,'fill');
    set(gca, 'LooseInset', [0,0,0,0]);
    hold on
    scatter(x, y,s,[0,0,0],'LineWidth',1);
    hold all
    line([-3.75, -3.75],[-4.5, -0.4],'LineStyle','--','Color','k')
    line([-3.25, -3.25],[-4.5, -0.4],'LineStyle','--','Color','k')
    line([-3.6:0.1:-3.4],zeros(length([-3.6:0.1:-3.4]))-2.5,'LineStyle','none','Marker','.','Color','k')
    
    TextSizes.DefaultAxesFontSize = 6;
    TextSizes.DefaultTextFontSize = 6;
    set(0,TextSizes);
    
    ylabel('standard deviation');
    xlabel('censoring interval \Deltat');
    caxis([0 100]);
    h = colorbar;
    ylabel(h,'% identified correctly')
    
    axis equal
    ylim([-4 -0.4])
    xlim([-4.2 0])
    set(gca,'XTick',[-4 -3 -2 -1 0 ])
    tck = get(gca,'XTick')';
    tck=10.^(tck);
    tck(1)=0;
    set(gca,'XTickLabel',tck)
    
    set(gca,'YTick',[-4 -3 -2 -1 0 ])
    tck = get(gca,'YTick')';
    tck=10.^(tck);
    
    set(gca,'YTickLabel',tck)
    
    
    title(ts);
    
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
    if i==1
        print('-depsc2','-r1000',['./figs/grid_BIC_cen_dis_BW']);
        print('-dpdf','-r1000',['./figs/grid_BIC_cen_dis_BW']);
    else
        print('-depsc2','-r1000',['./figs/grid_BIC_cen_BW']);
    end
end


% pro sigma ein plot
%% VISUALIZATION OF RESULTS
col_bar = [1 0 0; repmat(linspace(0.5,0.95,o_max-1)',1,3)];
%% set figure size
width2=10;
height2=10;
for p=2:length(sigma)
    % Legend constcurtion
    for o = 1:o_max
        if o == 1
            ls{o} = [num2str(o,'%d') ' component'];
        else
            ls{o} = [num2str(o,'%d') ' components'];
        end
    end
    
    % BIC:
    hFig2=figure('name','BIC');
    set(hFig2, 'Units','centimeters','Position', [10 10 width2 height2])
    for i = 1:2
        % Open subplot
        subplot(2,1,i)
        % Case
        switch i
            case 1, perc = perc_selection_BIC(:,:,p);  ts =['s.d.= ' num2str(sqrt((exp(sigma(p)^2)-1)*exp(2*1+sigma(p)^2)),'%6.3g') ' censoring disregarded'];
            case 2, perc =perc_selection_cBIC(:,:,p); ts =['s.d.= ' num2str(sqrt((exp(sigma(p)^2)-1)*exp(2*1+sigma(p)^2)),'%6.3g') ' censoring considered'];
                sigma(p)
        end
        % different levels of censoring
        for o = 1:o_max
            dx(1)=1e-4;
            if o == 1
                lh(o) = fill([log10(dx(1:end)),log10(dx(end:-1:1))],...
                    [zeros(1,length(dx)),perc(1,end:-1:1)],...
                    col_bar(o,:)); hold on;
            else
                lh(o) = fill([log10(dx(1:end)),log10(dx(end:-1:1))],...
                    [sum(perc(1:(o-1),1:end),1),sum(perc(1:o,end:-1:1),1)],...
                    col_bar(o,:)); hold on;
            end
            xlabel('censoring interval \Deltat');
            ylabel('percentage');
            ylim([0 1])
            set(gca,'XTick',[-4 -3 -2 -1 0 ])
            tck = get(gca,'XTick')';
            tck=10.^(tck);
            tck(1)=0;
            set(gca,'XTickLabel',tck)
            set(gca,'YTick',[0 0.5 1 ])
            tck = get(gca,'YTick')';
            tck=100*(tck);
            
            set(gca,'YTickLabel',tck)
        end
        title(ts);
        if p == 2 && i==1
            legend(lh,ls);
        end
    end
    
    
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
    print('-depsc2','-r1000',['./figs/',['detail_' num2str(p,'%d')],'_BIC']);
end