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
%
% Outputs:
% ========
% fh .. figure handle
%
% 2012/05/31 Jan Hasenauer

% function fh = plotP(parameters,fh,I,options)
function fh = plotP(varargin)

%% CHECK AND ASSIGN INPUTS
% Assign parameters
if nargin >= 1
    parameters = varargin{1};
else
    error('plotPL requires a parameter object as input.');
end

% Open figure
if nargin >= 2
    if ~isempty(varargin{2})
        fh = figure(varargin{2});
    else
        fh = figure;
    end
else
    fh = figure;
end

% Index of subplot which is updated 
if nargin >= 3
    I = varargin{3};
else
    I = 1:parameters.number;
end

% Options
options.mark_constraint = 'false';
options.interval = 'static';
options.hold_on = 'false';
options.linewidth = 1.5;
if nargin == 4
    options = setdefault(varargin{4},options);
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

    % Open subplot
    subplot(s(1),s(2),i);
    if strcmp(options.hold_on,'true');
        hold on;
    else
        hold off;
    end
    
    % Plot profile likelihood
    if ~isempty(parameters.P(i).par)
        plot(parameters.P(i).par(i,:),exp(parameters.P(i).logPost - parameters.MS.MAP.logPost),'r-','linewidth',options.linewidth); hold on;
    end
    
    % Determine index of points which are in the interior
    if strcmp(options.mark_constraint,'true')
        ind = find(sum(bsxfun(@gt,xmin,parameters.P(i).par)+bsxfun(@gt,parameters.P(i).par,xmax),1));
        if ~isempty(ind)
            plot(parameters.P(i).par(i,ind),exp(parameters.P(i).logPost(ind) - parameters.MS.MAP.logPos),'rx','linewidth',options.linewidth); hold on;    
        end
    end
    
    % Plot optimum
    plot(parameters.MS.MAP.par(i),1,'ro','linewidth',options.linewidth); hold on;
    % Limits
    switch options.interval
        case 'static'
            xlim([parameters.min(i),parameters.max(i)]);
        case 'dynamic'
            xlim([min(parameters.P(i).par(i,:)),max(parameters.P(i).par(i,:))]);
        case 'none'
            %
        otherwise
            error('This option is not available.');
    end
    ylim([0,1.1]);
    % Labels
    xlabel(parameters.name(i));
    if mod(i,s(2)) == 1
        ylabel('ration, R');
    else
        set(gca,'Ytick',[]);
    end

end

drawnow;
