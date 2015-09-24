% plotMixtureModel plots a mixture model and the corresponding data.
%   It allows for the consideration of multiple experiments
%
% USAGE:
% ======
% [fh] = plotMixtureModel(parameters,M,Mc,D)
% [fh] = plotMixtureModel(parameters,M,Mc,D,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing at least parameters.MS.MAP.par.
% M ... model structure
% D ... data structure
% options ... options for plot:
%   .fh ... figure handle. If no figure handle is provided, a new figure is
%           created.
%   .plot_type plot cdf or pdf of distribution (default = 'cdf')
%   .plot_model.subtype plot maximum a posteriori estimate ('MAP') or MCMC ('MCMC')    (default  = 'MAP')
%                       MCMC only valid for .plot_type plot = 'pdf':
%   .plot_subpop plot distributions of subpopulations(default  = 'false')
%                only available in combination with and .plot_model.subtype = 'MAP'.
%                 an dplot_model.subtype = 'MAP'.
%   .plot_mixture plot mixture distribution (default = 'true')
%   .plot_overall plot overall distribution (default = 'false'), option
%                 only available in combination with .plot_type = 'pdf'
%                 and plot_model.subtype = 'MAP'.
%   .plot_data.type plot data as histogram ('hist') or as single points ('single') (default ='single')
%   .plot_data.bins number of bins for histograms (default =30)
%   .xlim limits for x axis of subplots [min max], either for all subplots at once [min max] or for each subplot extra [min1 max1; min2 max2; ... ](default = [])
%
% Outputs:
% ========
% fh ... figure handle
%
% 2012/05/16 Jan Hasenauer
% 2012/07/11 Jan Hasenauer

% function [fh] = plotMixtureModel(parameters,M,D,options)
function [fh] = plotMixtureModel(varargin)

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
options.plot_type = 'cdf';  %'pdf'
options.plot_model.subtype = 'MAP'; % 'MCMC';
options.plot_subpop = 'false';
options.plot_mixture = 'true';
options.plot_overall = 'false';
options.plot_data.type ='single'; %'hist'
options.plot_data.bins=30;
options.xlim = [];
options.unc_ana.N = 1e3;

for j = 1:length(D)
    options.xmin(j) = min([D{j}.data.uncensored(:);D{j}.data.censored(:)]);
    options.xmax(j) = max([D{j}.data.uncensored(:);D{j}.data.censored(:)]);
end
if nargin == 5
    options = setdefault(varargin{5},options);
end

%% PREPARE FIGURE
if isempty(options.fh)
    fh = figure;
else
    fh = figure(options.fh);
end

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
            ones(1,M.experiment(j).size+1) ]';
        
        %% PLOT EVENT DENSITIES
        % Model
        switch options.plot_model.subtype
            case 'MAP'
                % Uncensored data
                lh(end+1) = stairs(x,Dens.experiment(1).Fec,'-','color',[0,0.5,0],'linewidth',2); hold on;
                ll{end+1} = ['M: uncensored'];
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = stairs(x,Dens.experiment(1).Fce,'-','color',[1.0,0,0],'linewidth',2); hold on;
                    ll{end+1} = ['M: censored'];
                end
                
                % Subpopulations
                if strcmp(options.plot_subpop,'true')
                    if M.experiment(j).size>= 2
                        for k = 1:M.experiment(j).size
                            lh(end+1) = stairs(x,Dens.experiment(1).Fekc{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                            ll{end+1} = ['M sub ' num2str(k,'%d') ': uncensored'];
                        end
                    end
                end
                
                %                 % plot overall event distribution
                %                 if strcmp(options.plot_overall,'true')
                %                     lh(end+1) = stairs(x,Dens.experiment(1).Fe,'--','color',[0,0.5,0],'linewidth',2); hold on;
                %                     ll{end+1} = ['M: overall'];
                %                 end
                
            case 'MCMC'
                % Determine index of median
                im = find(percMeas{1}.percentiles == 50);
                % Uncensored data
                lh(end+1) = stairs(x,percMeas{1}.Fec(:,im),'-','color',[0,0.5,0],'linewidth',2); hold on;
                for k = 1:(size(percMeas{1}.Fec,2)-1)/2
                    Xb = [x(1:end-1),x(end:-1:2);x(2:end),x(end-1:-1:1)];
                    Yb = [percMeas{1}.Fec(1:end-1,k)',percMeas{1}.Fec(end-1:-1:1,end+1-k)';...
                        percMeas{1}.Fec(1:end-1,k)',percMeas{1}.Fec(end-1:-1:1,end+1-k)'];
                    fill(Xb(:),Yb(:),[0,0.5,0],'facealpha',.1,'EdgeColor','none');
                end
                ll{end+1} = ['M: uncensored'];
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = stairs(x,percMeas{1}.Fce(:,im),'-','color',[1.0,0,0],'linewidth',2); hold on;
                    for k = 1:(size(percMeas{1}.Fce,2)-1)/2
                        Xb = [x(1:end-1),x(end:-1:2);x(2:end),x(end-1:-1:1)];
                        Yb = [percMeas{1}.Fce(1:end-1,k)',percMeas{1}.Fce(end-1:-1:1,end+1-k)';...
                            percMeas{1}.Fce(1:end-1,k)',percMeas{1}.Fce(end-1:-1:1,end+1-k)'];
                        fill(Xb(:),Yb(:),[1.0,0,0],'facealpha',.1,'EdgeColor','none');
                    end
                    ll{end+1} = ['M: censored'];
                end
        end
        
        % Data
        lh(end+1) = stairs(x,Fec_D,'-','color',[0,0.5,0]);
        ll{end+1} = ['D: uncensored'];
        if ~isempty(Mc)
            lh(end+1) = stairs(x,Fce_D,'-','color',[1,0,0]);
            ll{end+1} = ['D: censored'];
        end
        % Bounds
        xlim(x([1,end]));
        ylim([0,1.1]);
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
                ll{end+1} = ['M: uncensored'];
                % Censored data
                if ~isempty(Mc)
                    lh(end+1) = plot(x,Dens.experiment(1).Fce,'-','color',[1.0,0,0],'linewidth',2); hold on;
                    ll{end+1} = ['M: censored'];
                end
                % Subpopulations
                if strcmp(options.plot_subpop,'true')
                    if M.experiment(j).size>= 2
                        if ~isempty(Mc)
                            for k = 1:M.experiment(j).size
                                lh(end+1) = plot(x,Dens.experiment(1).Fekc{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                                ll{end+1} = ['M sub ' num2str(k,'%d') ': uncensored'];
                            end
                        else
                            for k = 1:M.experiment(j).size
                                lh(end+1) = plot(x,Dens.experiment(1).Fek{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                                ll{end+1} = ['M sub ' num2str(k,'%d') ': uncensored'];
                            end
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
                ll{end+1} = ['M: uncensored'];
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
        ll{end+1} = ['D: uncensored'];
        if ~isempty(Mc)
            lh(end+1) = stairs(x,Fce_D,'-','color',[1,0,0]);
            ll{end+1} = ['D: censored'];
        end
        
        % Bounds
        xlim(x([1,end]));
        ylim([0,1.1]);
        
        %%%%%%%%%%%
        %%%%%%%%%%%
        %%%%%%%%%%%
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
        col = 1*[zeros(1,M.experiment(j).size+1);
            linspace(0.6,0.9,M.experiment(j).size+1);
            zeros(1,M.experiment(j).size+1) ]';
        
        %% PLOT DATA
        
        if strcmp(options.plot_data.type,'hist')  %  as histogramms
            
            
            % Plot uncensored data
            if ~isempty(X)
                
                X_min = inf;
                X_max = -inf;
                
                X_min = min(X_min, min(X));
                X_max = max(X_max, max(X));
                X_hist = linspace(X_min,X_max, options.plot_data.bins+1);
                [n, h]=hist(X, X_hist);
                lh(end+1) = bar(h,n/abs(length(X)*(X_hist(2)-X_hist(1))), 'Edgecolor',[0.1 0.7 0], 'Facecolor', [0.1 0.7 0],'Barwidth',1); hold on
                ll{end+1} = ['D: uncensored'];
            end
            
            
            % Plot censored data
            if ~isempty(Xc)
                Xc_min = inf;
                Xc_max = -inf;
                
                Xc_min = min(Xc_min, min(Xc));
                Xc_max = max(Xc_max, max(Xc));
                Xc_hist = linspace(X_min,X_max, options.plot_data.bins+1);
                [n,h]=hist(Xc, Xc_hist);
                lh(end+1) = bar(h,n/abs(length(Xc)*(X_hist(2)-X_hist(1))), 'Edgecolor',[1 0.3 0.1], 'Facecolor', [1 0.3 0.1],'Barwidth',1); hold on
                ll{end+1} = ['D: censored'];
            end
            
            
        elseif strcmp(options.plot_data.type,'single') % as single points
            
            if ~strcmp(options.plot_model.subtype ,'MCMC')
                ymax = max(max(Dens.experiment(1).fec),max(Dens.experiment(1).fce));
                % Plot uncensored data
                if ~isempty(X)
                    X_red = unique(X);
                    for i = 1:length(X_red)
                        num_rep = sum(X_red(i) == X);
                        if i == 1
                            lh(end+1) = plot(X_red(i)*ones(num_rep,1),0.06*ymax*[0:1:num_rep-1],'o','color',[0,0.5,0]); hold on;
                            ll{end+1} = ['D: uncensored'];
                        else
                            plot(X_red(i)*ones(num_rep,1),0.03*ymax*[0:1:num_rep-1],'o','color',[0,0.5,0]); hold on;
                        end
                    end
                end
            else
                error('The option ''MCMC'' is only valid for ''cdfs'' !');
            end
            % Plot censored data
            if ~isempty(Xc)
                Xc_red = unique(Xc);
                for i = 1:length(Xc_red)
                    num_rep = sum(Xc_red(i) == Xc);
                    if i == 1
                        lh(end+1) = plot(Xc_red(i)*ones(num_rep,1),0.06*ymax*[0:1:num_rep-1],'r^'); hold on;
                        ll{end+1} = ['D: censored'];
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
                                lh(end+1) = plot(x,Dens.experiment(1).fekc{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                                ll{end+1} = ['M sub ' num2str(k,'%d') ': uncensored'];
                            end
                        else
                            for k = 1:M.experiment(j).size
                                lh(end+1) = plot(x,Dens.experiment(1).fek{k},'--','color',col(k+1,:),'linewidth',2); hold on;
                                ll{end+1} = ['M sub ' num2str(k,'%d') ': uncensored'];
                            end
                        end
                    end
                end
                if strcmp(options.plot_mixture,'true')
                    % Uncensored data
                    lh(end+1) = plot(x,Dens.experiment(1).fec,'-','color',[0,0.5,0],'linewidth',2); hold on;
                    ll{end+1} = ['M: uncensored'];
                    % Censored data
                    if ~isempty(Mc)
                        lh(end+1) = plot(x,Dens.experiment(1).fce,'-','color',[1.0,0,0],'linewidth',2); hold on;
                        ll{end+1} = ['M: censored'];
                    end
                end
                
                if strcmp(options.plot_overall,'true')
                    % Uncensored data
                    lh(end+1) = plot(x,Dens.experiment(1).fe,'--','color',[0,0.5,0],'linewidth',2); hold on;
                    ll{end+1} = ['M'];
                    % Censored data
                    if ~isempty(Mc)
                        lh(end+1) = plot(x,Dens.experiment(1).fc,'--','color',[1.0,0,0],'linewidth',2); hold on;
                        ll{end+1} = ['Mc'];
                    end
                end
                
            case 'MCMC'
                error('The option ''MCMC'' is only valid for ''cdfs'' !');
                
        end
        
        xlim([xmin,xmax]);
        ymax = max(max(Dens.experiment(1).fec),max(Dens.experiment(1).fce));
    end
    
    %% LABEL AND BOUNDS
    % Labels
    xlabel(M.label.x);
    ylabel(M.label.y);
    if ~isempty(options.xlim)
        if size(options.xlim,1) >= length(D)
            xlim(options.xlim(j,:));
        else
            xlim(options.xlim(1,:));
        end
    end
    title(M.experiment(j).name);
    %if j == 1
    legend(lh,ll,'Location','NorthWest');
    %end
    
end

drawnow;