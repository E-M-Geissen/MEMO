% plotMixtureModel plots a mixture model and the corresponding data.
%   It allows for the consideration of multiple experiments
%
% USAGE:
% ======
% [fh] = plotMixtureModel(parameters,M,D)
% [fh] = plotMixtureModel(parameters,M,D,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing at least parameters.MS.MAP.par.
% M ... model structure
% D ... data structure
% options ... options for plottion:
%   .fh ... figure handle. If no figure handle is provided, a new figure is
%           created.
%
% Outputs:
% ========
% fh ... figure handle
%
% 2012/05/16 Jan Hasenauer
% 2012/07/11 Jan Hasenauer

% function [fh] = plotMixtureModel(parameters,M,D,options)
function [fh] = plotMixtureModel3(varargin)


TextSizes.DefaultAxesFontSize = 8;
TextSizes.DefaultTextFontSize =8;
set(0,TextSizes);


%% CHECK AND ASSIGN INPUTS
if nargin >= 3
    parameters = varargin{1};
    M = varargin{2};
    Mc = varargin{3};
    D = varargin{4};
else
    error('Not enought inputs.')
end

options.fh = [];
options.plot_type = 'cdf';
options.plot_data.style = 'discrete';
options.plot_data.scaling = 0.3*1e-1;
options.plot_model.style = 'continuous';
options.plot_model.type = 'cumulative';
options.plot_model.subtype = 'MAP';
options.plot_model.subpop = 'true';
options.plot_mixture='true';
%options.plot_model.subtype = 'MCMC';
options.unc_ana.N = 1e3;
options.width = 18;
options.height = 12;

for j = 1:length(D)
    options.xmin(j) = min([D{j}.data.uncensored(:);D{j}.data.censored(:)]);
    options.xmax(j) = max([D{j}.data.uncensored(:);D{j}.data.censored(:)]);
end
if nargin == 5
    options = setdefault(varargin{5},options);
end
width=options.width;
height=options.height;
%% PREPARE FIGURE
if isempty(options.fh)
    fh = figure;
else
    fh = figure(options.fh);
end

set(fh, 'Units','centimeters','Position', [5 5 width height])

s = round(sqrt(length(D))*[1,1]);
if prod(s) < length(D)
    s(2) = s(2) + 1;
end

% Loop: data
for j = 1:length(D)
    
    lh = [];
    ll = {};
    
    %% DATA
    X  = D{j}.data.uncensored(:);
    Xc = D{j}.data.censored(:);
    dx = D{j}.observation_interval;
    ncells = length(X) + length(Xc);
    
    %% LOWER AND UPPER BOUNDS
    
    if min(X)<0
        xmin=min(X)-abs(min(X)/5);
    else
        xmin=0;
    end
    
    xmax = -inf;
    if ~isempty(X)
        xmax = max(xmax,max(X));
    end
    if ~isempty(Xc)
        xmax = max(xmax,max(Xc));
    end
    xmax = 10^(ceil(log10(xmax)*10+1)/10);
    
    %%%%%%%%%%%
    %%%%%%%%%%%
    %%%%%%%%%%%
    if strcmp(options.plot_type,'cdf') && (dx > 0)
        %% GRID
        x = xmin:dx:xmax+dx;
        
        %% EVALUATION OF DENSITIES AND CUMULATIVE DENSITIES
        % Data
        Fec_D = cumsum(histc(X,x))/ncells;
        Fce_D = cumsum(histc(Xc,x)/ncells);
        
        % Model
        % Construction of model for experiment j
        M_j.mixture = M.mixture;
        M_j.experiment(1) = M.experiment(j);
        if ~isempty(Mc)
            Mc_j.mixture = Mc.mixture;
            Mc_j.experiment(1) = Mc.experiment(j);
        else
            Mc_j = [];
        end
        D_j{1} = D{j};
        % Calculation of densities
        switch options.plot_model.subtype
            case 'MAP'
                Dens = getMixtureDensity(x,parameters.MS.MAP.par,M_j,Mc_j);
            case 'MCMC'
                [percMeas,percModelPred] = computeCondifenceIntervals(parameters,M_j,Mc_j,D_j,options.unc_ana);
        end
        
        %% PLOT INDIVIDUAL EXPERIMENTS
        % Open subplot
        subplot(s(1),s(2),j);
        col = 1*[linspace(0,0.7,M.experiment(j).size+1);
            linspace(0,0.7,M.experiment(j).size+1);
            zeros(1,M.experiment(j).size+1) ]';
        
        %% PLOT EVENT DENSITIES
        % Model
        switch options.plot_model.subtype
            case 'MAP'
                % Uncensored data
                lh(end+1) = stairs(x,Dens.experiment(1).Fec,'-','color',[0,0,0],'linewidth',1); hold on;
                ll{end+1} = ['Model: counts'];
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = stairs(x,Dens.experiment(1).Fce,'-','color',[169/255,169/255,169/255],'linewidth',1); hold on;
                    ll{end+1} = ['Model: censoring times'];
                    
                    % Subpopulations
                    if M.experiment(j).size>= 2 && strcmp(options.plot_model.subpop ,'true')
                        for k = 1:M.experiment(j).size
                            lh(end+1) = stairs(x,Dens.experiment(1).Fekc{k},'-','color',col(k+1,:),'linewidth',2); hold on;
                            ll{end+1} = ['M subpopulation ' num2str(k,'%d')];
                        end
                    end
                    
                else
                    
                    % Subpopulations
                    if M.experiment(j).size>= 2 && strcmp(options.plot_model.subpop ,'true')
                        for k = 1:M.experiment(j).size
                            lh(end+1) = stairs(x,Dens.experiment(1).Fek{k},'--','color',col(k+1,:),'linewidth',1); hold on;
                            ll{end+1} = ['M subpopulation ' num2str(k,'%d') ];
                        end
                    end
                end
            case 'MCMC'
                % Determine index of median
                im = find(percMeas{1}.percentiles == 50);
                % Uncensored data
                for k = 1:(size(percMeas{1}.Fec,2)-1)/2
                    Xb = [x(1:end-1),x(end:-1:2);x(2:end),x(end-1:-1:1)];
                    Yb = [percMeas{1}.Fec(1:end-1,k)',percMeas{1}.Fec(end-1:-1:1,end+1-k)';...
                        percMeas{1}.Fec(1:end-1,k)',percMeas{1}.Fec(end-1:-1:1,end+1-k)'];
                    fill(Xb(:),Yb(:),[87/255,87/255,87/255],'facealpha',.6,'EdgeColor','none');hold on;
                end
                
                ll{end+1} = ['Model: prometaphase lengths'];
                % Censored data
                if ~isempty(Mc)
                    for k = 1:(size(percMeas{1}.Fce,2)-1)/2
                        Xb = [x(1:end-1),x(end:-1:2);x(2:end),x(end-1:-1:1)];
                        Yb = [percMeas{1}.Fce(1:end-1,k)',percMeas{1}.Fce(end-1:-1:1,end+1-k)';...
                            percMeas{1}.Fce(1:end-1,k)',percMeas{1}.Fce(end-1:-1:1,end+1-k)'];
                        % plot(Xb,Yb,'c-');
                        fill(Xb(:),Yb(:),[210/255,210/255,210/255],'facealpha',1,'EdgeColor','none');
                        lh(end+1) = stairs(x,percMeas{1}.Fec(:,im),'-','color',[0,0,0],'linewidth',1); hold on;
                        lh(end+1) = stairs(x,percMeas{1}.Fce(:,im),'-','color',[169/255,169/255,169/255],'linewidth',1); hold on;
                    end
                    ll{end+1} = ['Model: censoring times'];
                end
        end
        
        % Data
        lh(end+1) = stairs(x,Fec_D,'-','color',[0,0,0],'linewidth',2);
        ll{end+1} = ['Data: counts'];
        if ~isempty(Mc)
            lh(end+1) = stairs(x,Fce_D,'-','color',[169/255,169/255,169/255],'linewidth',2);
            ll{end+1} = ['Data: censoring times'];
        end
         xlim([0 x(end)]);                                                                                                                                     
        % ylim([0,1.1]);
        %%%%%%%%%%%
        %%%%%%%%%%%
        %%%%%%%%%%%
    elseif strcmp(options.plot_type,'cdf') && (dx == 0)
        %% GRID
        x = linspace(xmin,xmax,200);
        
        %% EVALUATION OF DENSITIES AND CUMULATIVE DENSITIES
        % Data
        Fec_D = cumsum(histc(X,x))/ncells;
        Fce_D = cumsum(histc(Xc,x)/ncells);
        
        % Model
        % Construction of model for experiment j
        M_j.mixture = M.mixture;
        M_j.experiment(1) = M.experiment(j);
        if ~isempty(Mc)
            Mc_j.mixture = Mc.mixture;
            Mc_j.experiment(1) = Mc.experiment(j);
        else
            Mc_j = [];
        end
        D_j{1} = D{j};
        % Calculation of densities
        switch options.plot_model.subtype
            case 'MAP'
                Dens = getMixtureDensity(x,parameters.MS.MAP.par,M_j,Mc_j);
            case 'MCMC'
                [percMeas,percModelPred] = computeCondifenceIntervals(parameters,M_j,Mc_j,D_j,options.unc_ana);
        end
        
        %% PLOT INDIVIDUAL EXPERIMENTS
        % Open subplot
        subplot(s(1),s(2),j);
        col = 1*[linspace(0,0.7,M.experiment(j).size+1);
            linspace(0,0.7,M.experiment(j).size+1);
            ones(1,M.experiment(j).size+1) ]';
        
        %% PLOT EVENT DENSITIES
        % Model
        switch options.plot_model.subtype
            case 'MAP'
                % Uncensored data
                lh(end+1) = plot(x,Dens.experiment(1).Fec,'-','color',[0,0.5,0],'linewidth',2); hold on;
                ll{end+1} = ['Model'];
                %                 % Censored data
                %                 if ~isempty(Mc)
                %                     lh(end+1) = plot(x,Dens.experiment(1).Fce,'-','color',[1.0,0,0],'linewidth',2); hold on;
                %                     ll{end+1} = ['M: censored'];
                %                 end
                %                 % Subpopulations
                %                 if M.experiment(j).size>= 2
                %                     for k = 1:M.experiment(j).size
                %                         lh(end+1) = plot(x,Dens.experiment(1).Fekc{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                %                         ll{end+1} = ['M sub ' num2str(k,'%d') ': uncensored'];
                %                     end
                %                 end
                
                
                
                
                
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = plot(x,Dens.experiment(1).Fce,'-','color',[1.0,0,0],'linewidth',2); hold on;
                    ll{end+1} = ['M: censored'];
                    
                    % Subpopulations
                    if M.experiment(j).size>= 2
                        for k = 1:M.experiment(j).size
                            lh(end+1) = plot(x,Dens.experiment(1).Fekc{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                            ll{end+1} = ['Model subpop. ' num2str(k,'%d') ];
                        end
                    end
                    
                else
                    
                    % Subpopulations
                    if M.experiment(j).size>= 2
                        for k = 1:M.experiment(j).size
                            lh(end+1) = plot(x,Dens.experiment(1).Fek{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                            ll{end+1} = ['Model subpop. ' num2str(k,'%d')];
                        end
                    end
                end
                
                
                
                
                
                
            case 'MCMC'
                % Determine index of median
                im = find(percMeas{1}.percentiles == 50);
                % Uncensored data
                lh(end+1) = plot(x,percMeas{1}.Fec(:,im),'-','color',[0,0.5,0],'linewidth',2); hold on;
                for k = 1:(size(percMeas{1}.Fec,2)-1)/2
                    Xb = [x(1:end),x(end:-1:1)];
                    Yb = [percMeas{1}.Fec(1:end,k)',percMeas{1}.Fec(end:-1:1,end+1-k)'];
                    fill(Xb(:),Yb(:),[0,0.5,0],'facealpha',.1,'EdgeColor','none');
                end
                ll{end+1} = ['Model'];
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = plot(x,percMeas{1}.Fce(:,im),'-','color',[1.0,0,0],'linewidth',2); hold on;
                    for k = 1:(size(percMeas{1}.Fce,2)-1)/2
                        Xb = [x(1:end),x(end:-1:2)];
                        Yb = [percMeas{1}.Fce(1:end,k)',percMeas{1}.Fce(end:-1:1,end+1-k)'];
                        % plot(Xb,Yb,'c-');
                        fill(Xb(:),Yb(:),[1.0,0,0],'facealpha',.1,'EdgeColor','none');
                    end
                    ll{end+1} = ['M: censored'];
                end
        end
        
        % Data
        lh(end+1) = stairs(x,Fec_D,'-','color',[0,0.5,0]);
        ll{end+1} = ['Data'];
        if ~isempty(Mc)
            lh(end+1) = stairs(x,Fce_D,'-','color',[1,0,0]);
            ll{end+1} = ['D: censored'];
        end
    elseif strcmp(options.plot_type,'pdf')
        
        %% GRID
        x = linspace(xmin,xmax,400);
        
        %% EVALUATION OF DENSITIES AND CUMULATIVE DENSITIES
        
        % Model
        % Construction of model for experiment j
        M_j.mixture = M.mixture;
        M_j.experiment(1) = M.experiment(j);
        if ~isempty(Mc)
            Mc_j.mixture = Mc.mixture;
            Mc_j.experiment(1) = Mc.experiment(j);
        else
            Mc_j = [];
        end
        D_j{1} = D{j};
        % Calculation of densities
        switch options.plot_model.subtype
            case 'MAP'
                Dens = getMixtureDensity(x,parameters.MS.MAP.par,M_j,Mc_j);
            case 'MCMC'
                [percMeas,percModelPred] = computeCondifenceIntervals(parameters,M_j,Mc_j,D_j,options.unc_ana);
        end
        
        %% PLOT INDIVIDUAL EXPERIMENTS
        % Open subplot
        subplot(s(1),s(2),j);
        col = 1*[linspace(0,0.7,M.experiment(j).size+1);
            linspace(0,0.7,M.experiment(j).size+1);
            zeros(1,M.experiment(j).size+1) ]';
        %% plot data
        % histogramms
        if strcmp(options.plot_data.type,'hist')
            
            
            
            % Plot uncensored data
            
            
            if ~isempty(X)
                
                X_min = inf;
                X_max = -inf;
                
                X_min = min(X_min, min(X));
                X_max = max(X_max, max(X));
                X_hist = linspace(X_min,X_max, options.plot_data.bins+1);
                h=hist(X, X_hist);
                lh(end+1) = bar(X_hist,h/abs(length(X)*(X_hist(2)-X_hist(1))), 'Edgecolor',[0.5 0.5 0.5], 'Facecolor', [0.5 0.5 0.5],'Barwidth',1); hold on
                ll{end+1} = ['data'];
            end
            % Plot censored data
            
            
            if ~isempty(Xc)
                Xc_min = inf;
                Xc_max = -inf;
                
                Xc_min = min(Xc_min, min(Xc));
                Xc_max = max(Xc_max, max(Xc));
                
                
                h=hist(Xc, X_hist);
                bar(Xc_hist,h/abs(length(Xc)*(X_hist(2)-X_hist(1))),'Edgecolor',[1 0.5 0.5], 'Facecolor', [1 0.5 0.5],'Barwidth',1); hold on
                ll{end+1} = ['rightcensored data'];
            end
            
            % plot single points of raw data
        elseif strcmp(options.plot_datatype,'raw')
            
            % Plot uncensored data
            if ~isempty(X)
                X_red = unique(X);
                for i = 1:length(X_red)
                    num_rep = sum(X_red(i) == X);
                    if i == 1
                        lh(end+1) = plot(X_red(i)*ones(num_rep,1),0.06*ymax*[0:1:num_rep-1],'o','color',[0,0.5,0]); hold on;
                        ll{end+1} = ['uncensored data'];
                    else
                        plot(X_red(i)*ones(num_rep,1),0.03*ymax*[0:1:num_rep-1],'o','color',[0,0.5,0]); hold on;
                    end
                end
            end
            
            % Plot censored data
            if ~isempty(Xc)
                Xc_red = unique(Xc);
                for i = 1:length(Xc_red)
                    num_rep = sum(Xc_red(i) == Xc);
                    if i == 1
                        lh(end+1) = plot(Xc_red(i)*ones(num_rep,1),0.06*ymax*[0:1:num_rep-1],'r^'); hold on;
                        ll{end+1} = ['rightcensored data'];
                    else
                        plot(Xc_red(i)*ones(num_rep,1),0.03*ymax*[0:1:num_rep-1],'r^'); hold on;
                    end
                end
            end
        end
        %% PLOT EVENT DENSITIES
        % Model
        switch options.plot_model.subtype
            case 'MAP'
                % Subpopulations
                if strcmp(options.plot_subpop,'true')
                    if M.experiment(j).size>= 2
                        if ~isempty(Mc)
                            for k = 1:M.experiment(j).size
                                lh(end+1) = plot(x,Dens.experiment(1).fekc{k},'color',col(k+1,:),'linewidth',2); hold on;
                                ll{end+1} = ['subpopulation ' num2str(k,'%d') ];
                            end
                        else
                            for k = 1:M.experiment(j).size
                                lh(end+1) = plot(x,Dens.experiment(1).fek{k},'color',col(k+1,:),'linewidth',2); hold on;
                                ll{end+1} = ['subpopulation ' num2str(k,'%d') ];
                            end
                        end
                    end
                end
                if strcmp(options.plot_mixture,'true')
                    % Uncensored data
                    lh(end+1) = plot(x,Dens.experiment(1).fec,'-','color',[0,0.5,0],'linewidth',2); hold on;
                    ll{end+1} = ['overall'];
                    % Censored data
                    if ~isempty(Mc)
                        lh(end+1) = plot(x,Dens.experiment(1).fce,'-','color',[1.0,0,0],'linewidth',2); hold on;
                        ll{end+1} = ['M: censored'];
                    end
                end
            case 'MCMC'
                % Determine index of median
                im = find(percMeas{1}.percentiles == 50);
                % Uncensored data
                lh(end+1) = plot(x,percMeas{1}.fec(:,im),'-','color',[0,0.5,0],'linewidth',2); hold on;
                for k = 1:(size(percMeas{1}.fec,2)-1)/2
                    Xb = [x(1:end),x(end:-1:1)];
                    Yb = [percMeas{1}.fec(1:end,k)',percMeas{1}.fec(end:-1:1,end+1-k)'];
                    fill(Xb(:),Yb(:),[0,0.5,0],'facealpha',.1,'EdgeColor','none');
                end
                ll{end+1} = ['M: uncensored'];
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = plot(x,percMeas{1}.fce(:,im),'-','color',[1.0,0,0],'linewidth',2); hold on;
                    for k = 1:(size(percMeas{1}.fce,2)-1)/2
                        Xb = [x(1:end),x(end:-1:2)];
                        Yb = [percMeas{1}.fce(1:end,k)',percMeas{1}.fce(end:-1:1,end+1-k)'];
                        % plot(Xb,Yb,'c-');
                        fill(Xb(:),Yb(:),[1.0,0,0],'facealpha',.1,'EdgeColor','none');
                    end
                    ll{end+1} = ['M: censored'];
                end
        end
        
        ymax = max(max(Dens.experiment(1).fec),max(Dens.experiment(1).fce));
    end
    
    %% LABEL AND BOUNDS
    % Labels
    xlabel(M.label.x);
    ylabel('cumulative probability');
    title(M.experiment(j).name);
   % if j==length(D)
        legend(lh,ll,'Location','East');
    %end
     %  set(gca,'xtick',10.^[-2:1],'xticklabel',[0.01,0.1,1,10]);
   %set(gca,'xtick',10.^[-1:1],'xticklabel',[0.1,1,10]);
    %     % Bounds
        %xlim([-1,2]);
         %set(gca,'xtick',[-1:2],'xticklabel',[0.1,1,10,100]);
    %     ylim([0,1.1]);
    
end

drawnow;