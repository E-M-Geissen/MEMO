% plotP plots the profiles stored in parameters.
%
% USAGE:
% ======
% fh = plotP(parameters)
% fh = plotP(parameters,fh)
% fh = plotP(parameters,fh,I)
% fh = plotP(parameters,fh,I,options)
%
% INPUTS:
% =======
% parameters ... parameter struct containing information about parameters
%       and profiles.
% fh ... handle of figure in which profiles is plotted. If no
%       figure handle is provided, a new figure is opened.
% I ... index of subplot (parameter) which is updated. If no index is
%       provided the profiles for all parameters are updated.
% options ... options of plotting
%       .mark_constraint ... if 'true', points on the profile for which one
%           ore more parameters are close to the constraints are marked
%           with a cross (default = 'false').
%       .scale ... scale of y-axis (default = 'lin').
%
% Outputs:
% ========
% fh .. figure handle
%
% 2012/05/31 Jan Hasenauer

% function fh = plotP(parameters,fh,I,options)
function fh = plotP(varargin)

TextSizes.DefaultAxesFontSize = 8;
TextSizes.DefaultTextFontSize =8;
set(0,TextSizes);


%% CHECK AND ASSIGN INPUTS
% Assign parameters
if nargin >= 1
    parameters = varargin{1};
    parameters_mcmc = varargin{2};
else
    error('plotPL requires a parameter object as input.');
end

% Open figure
if nargin >= 3
    if ~isempty(varargin{3})
        fh = figure(varargin{3});
    else
        fh = figure;
    end
else
    fh = figure;
end


% Index of subplot which is updated
if nargin >= 4
    I = varargin{4};
else
    I = 1:parameters.number;
end

% Options
options.mark_constraint = 'false';
options.alpha_level = 0.99;
options.type = 'R';
options.interval = 'dynamic';%'static';
options.hold_on = 'false';
options.linewidth = 1.5;
options.parameter_number = [];
options.plot_local = 'local';
options.bins = 50;
if nargin == 5
    options = setdefault(varargin{5},options);
end

width=options.width;
height=options.height;
set(fh, 'Units','centimeters','Position', [5 5 width height])

if ~isempty(options.parameter_number)
    npar_alpha = options.parameter_number;
else
    npar_alpha = parameters.number;
end


%% CONSTRUCT INTERIOR BOUNDS
p = 1e-4;
xmin = (1-p)*parameters.min +    p *parameters.max;
xmax =    p *parameters.min + (1-p)*parameters.max;

%% PLOT PROFILE LIKELIHOODS
% Compute number of subfigure
s = round(sqrt(parameters.number)*[1,1]);
if prod(s) < parameters.number
    s(2) = s(2) + 1;
end

% Loop: Parameter
for i = I
        lh = [];
    ll = {};
    % Open subplot
    subplot(s(1),s(2),i);
    if strcmp(options.hold_on,'true');
        hold on;
    else
        hold off;
    end
    
    % Plot profile likelihood
    if ~isempty(parameters.P(i).par)
        switch options.type
            case 'R'
                % Histograms
                [N,X] = hist(parameters_mcmc.MCMC.sample.par(i,:),options.bins);hold on;
                    lh(end+1)=bar(X,N/max(N),1);hold on;
                    ll{end+1} = ['MCMC sample'];
                    lh(end+1)=plot(parameters.P(i).par(i,:),exp(parameters.P(i).logPost - parameters.MS.MAP.logPost),'r-','linewidth',options.linewidth); hold on;
                    ll{end+1} = ['$\mathrm{PL}(\theta_i)/P(D|\theta^\mathrm{ML})$'];
            case 'PL'
                plot(parameters.P(i).par(i,:),parameters.P(i).logPost,'r-','linewidth',options.linewidth); hold on;
        end
    end
    
    
    % Determine index of points which are in the interior
    if strcmp(options.mark_constraint,'true')
        ind = find(sum(bsxfun(@gt,xmin,parameters.P(i).par)+bsxfun(@gt,parameters.P(i).par,xmax),1));
        if ~isempty(ind)
            switch options.type
                case 'R'
                    plot(parameters.P(i).par(i,ind),exp(parameters.P(i).logPost(ind) - parameters.MS.MAP.logPost),'rx','linewidth',options.linewidth); hold on;
                case 'PL'
                    plot(parameters.P(i).par(i,ind),parameters.P(i).logPost(ind),'rx','linewidth',options.linewidth); hold on;
            end
        end
    end
    
    % Plot optimum
    switch options.type
        case 'R'
           lh(end+1)= plot(parameters.MS.MAP.par(i),1,'ro','linewidth',options.linewidth); hold on;
           ll{end+1} = ['MLE'];
        case 'PL'
            plot(parameters.MS.MAP.par(i),parameters.MS.MAP.logPost,'ro','linewidth',options.linewidth); hold on;
            %
            ind = find(parameters.MS.MAP_list.logPost >= parameters.MS.MAP.logPost-chi2inv(options.alpha_level,npar_alpha)/2);
            plot(parameters.MS.MAP_list.par(i,ind),parameters.MS.MAP_list.logPost(ind),'ro','linewidth',options.linewidth); hold on;
    end
    
    % Plot confidence levels
    
    % Limits
    switch options.interval
        case 'static'
            xl = [parameters.min(i),parameters.max(i)];
        case 'dynamic'
            xl = [min([parameters.P(i).par(i,:),min(parameters.MS.MAP_list.par(i,:))]),...
                max([parameters.P(i).par(i,:),max(parameters.MS.MAP_list.par(i,:))])];
        otherwise
            error('This option is not available.');
    end
    xlim(xl);
    switch options.type
        case 'R'
            ylim([0,1.1]);
        case 'PL'
            ylim([parameters.MS.MAP.logPost-4,parameters.MS.MAP.logPost]);
    end
    
    % Confidence levels
    switch options.type
        case 'R'
            plot(xl,[1,1]*exp(-chi2inv(0.99,1)/2),'k--');
            plot(xl,[1,1]*exp(-chi2inv(0.99,parameters.number)/2),'k:');
        case 'PL'
            plot(xl,[1,1]*(parameters.MS.MAP.logPost-chi2inv(options.alpha_level,1)/2),'k--');
            plot(xl,[1,1]*(parameters.MS.MAP.logPost-chi2inv(options.alpha_level,npar_alpha)/2),'k:');
    end
    
    % Plot local approximation
    if strcmp(options.plot_local,'true')
        H = parameters.MS.MAP.hessian(i,i);
        x = parameters.P(i).par(i,:);
        x_opt = parameters.MS.MAP.par(i);
        switch options.type
            case 'R'
                plot(x,exp(-0.5*H*(x-x_opt).^2),'b-','linewidth',options.linewidth); hold on;
            case 'PL'
                plot(x,parameters.MS.MAP.logPost-0.5*H*(x-x_opt).^2,'b-','linewidth',options.linewidth); hold on;
        end
    end
    
    % Labels
    xlabel(parameters.name(i));
    if mod(i,s(2)) == 1
        switch options.type
            case 'R'
                %ylabel('ration, R');
            case 'PL'
                ylabel('$\mathrm{PL}(\theta_i)$','interpreter','latex');
        end
    else
        set(gca,'Ytick',[]);
    end
    if i==max(I);
    legend(lh,ll,'Location','EastOutside','interpreter','latex');
    end
end

drawnow;
