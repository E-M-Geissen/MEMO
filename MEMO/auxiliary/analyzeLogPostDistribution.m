% analyzeLogPostDistribution estimates the expected distribution of
%   likelihood function values, assuming that the model M with parameters
%   theta describes the processes.
%   To estimate this distribution, artificial measurement data are 
%   generated from the mixture model, and the log-posterior function 
%   is evaluated.
%
% USAGE:
% ======
% [logPost,fh] = analyzeLogPostDistribution(theta,M,D)
% [logPost,fh] = analyzeLogPostDistribution(theta,M,D,options)
%
% INPUTS:
% =======
% parameters ... parameter object
% M ... mixture model for given data
% D ... measurement data
% options ... options of this algorithm
%   .N ... size of sample drawn from distribution (default = 1000)
%   .bins ... number of bins of the plotted histogram (default = 25)
%   .plot ... plot the results during the computation (default = 'true').
%   .fh ... figure handle. If no figure handle is provided, a new figure
%           is created.
%
% Outputs:
% ========
% logL ... sample of log-likelihood values from the model
% fh ... figure handle
%
% 2012/05/18 Jan Hasenauer

% function [logPost,fh] = analyzeLogPostDistribution(theta,M,D,options)
function [logPost,fh] = analyzeLogPostDistribution(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 4
    parameters = varargin{1};
    M = varargin{2};
    Mc = varargin{3};
    D = varargin{4};
else
    error('analyzeLogLDistribution requires at least four inputs.')
end

% Set default and assign option
options.N = 1000;
options.bins = 25;
options.plot = 'true';
options.reoptimize = 'true';
options.fh = [];
options.optim_opt.n_starts = 0;
options.optim_opt.plot = 'false';
options.optim_opt.mode = 'silent';
options.optim_opt.fmincon = optimset('algorithm','active-set',...
                                     'display','off',...
                                     'GradObj','on',...
                                     'MaxIter',1000,...
                                     'MaxFunEvals',1000*parameters.number);
if nargin == 5
    options = setdefault(varargin{5},options);
end

%% INITIALIZATION
logPost = [];
par.number = parameters.number;
par.name = parameters.name;
par.guess = parameters.MS.MAP.par;
par.min = parameters.min;
par.max = parameters.max;
% Figure handle
if strcmp(options.plot,'true')
    if isempty(options.fh)
        fh = figure;
    else
        fh = figure(options.fh);
    end
else
    fh = [];
end


%% GENERATE ARTIFICIAL SAMPLE OF LOG-LIKELIHOOD AND REOPTIMIZE PARAMETERS
% Loop: Realization of artificial data
for i = 1:options.N
    % Generate artificial data
    aD = generateArtificialData(parameters.MS.MAP.par,M,Mc,D);
    aD = processData(aD);
    switch options.reoptimize
        case 'false'
            % Evaluate likelihood function
            if isempty(Mc)
                logPost = [logPost;logLikelihoodMM(parameters.MS.MAP.par,M,aD)];
            else
                logPost = [logPost;logLikelihoodMMc(parameters.MS.MAP.par,M,Mc,aD)];
            end
        case 'true'
            % Reoptimization
            if isempty(Mc)
                logLikelihood = @(theta,opt) logLikelihoodMM(theta,M,aD,opt);
            else
                logLikelihood = @(theta,opt) logLikelihoodMMc(theta,M,Mc,aD,opt);
            end
            pari = optimizeMultiStart(par,logLikelihood,options.optim_opt);
            logPost = [logPost;pari.MS.MAP.logPost];
            % Output
            disp([num2str(i) '/' num2str(options.N)]);
        otherwise
            error('This option is not available.')
    end
    
    %% PLOT STATISTIC OF LOG-LIKELIHOOD
    if strcmp(options.plot,'true')
        % Open figure
        figure(fh); hold off;
        % Plot histogram
        [N,X] = hist(logPost,options.bins);
        bar(X,N/max(N),1,'facecolor',0.7*[1,1,1]); hold on;
        lh = plot(parameters.MS.MAP.logPost,0,'ro','markersize',8,'linewidth',2);
        % Label, legend, and ticks
        xlabel('log-likelihood');
        ylabel('relative frequency');
        legend(lh,'estimated logL');
        drawnow;
    end    
end
disp(' ');
disp(['Sampling of log-likelihood distribution. -> done']);


    