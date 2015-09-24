
% function [logPost,fh] = computeCondifenceIntervals(parameters,M,Mc,D,options)
function [percMeas,percModelPred] = computeCondifenceIntervals(varargin)

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
options.N = 100;
options.plot = 'true';
options.fh = [];
options.conf_level = [5,50,95];
if nargin == 5
    options = setdefault(varargin{5},options);
end

%% INITIALIZATION
%
for j = 1:length(D)
    % Data
    X  = D{j}.data.uncensored(:);
    Xc = D{j}.data.censored(:);
    dx = D{j}.observation_interval;
    ncells(j) = length(X) + length(Xc);
    % Lower and upper bounds
    xmin = 0;
    xmax = -inf;
    if ~isempty(X)
        xmax = max(xmax,max(X));
    end
    if ~isempty(Xc)
        xmax = max(xmax,max(Xc));
    end
    xmax = 10^(ceil(log10(xmax)*10+1)/10);
    % Grid
    x{j} = xmin:dx:xmax+dx;
    % Fine grid
    r = ceil(400/length(x{j}));
    xf{j} = xmin:dx/r:xmax+dx;
    
    % Sample of measured cumulative distribution
    measCDF{j}.Fec = zeros(length(x{j}),options.N);
    measCDF{j}.Fce = zeros(length(x{j}),options.N);
    % Sample of model pdfs
    modelPred{j}.fec = zeros(length(xf{j}),options.N);
    modelPred{j}.fce = zeros(length(xf{j}),options.N);
    modelPred{j}.Fec = zeros(length(xf{j}),options.N);
    modelPred{j}.Fce = zeros(length(xf{j}),options.N);
end

%% EVALUATION OF MEASUREMENT CONFIDENCE INTERVALS
% Loop: Realization of artificial data and parameter value
for i = 1:options.N
    % Generate artificial data
    theta = parameters.MCMC.sample.par(:,randi(size(parameters.MCMC.sample.par,2)));
    aD = generateArtificialData(theta,M,Mc,D);
    % Evaluate artificial date
    for j = 1:length(D)
        % Percentile intervals of measured cumulative distribution
        H  = histc(aD{j}.data.uncensored,x{j});
        Hc = histc(aD{j}.data.censored  ,x{j});
        measCDF{j}.Fec(:,i) = cumsum(H )'/ncells(j);
        measCDF{j}.Fce(:,i) = cumsum(Hc)'/ncells(j);
    end
end
% Evaluation of percentile intervals
for j = 1:length(D)
    percMeas{j}.Fec = prctile(measCDF{j}.Fec,options.conf_level,2);
    percMeas{j}.Fce = prctile(measCDF{j}.Fce,options.conf_level,2);
    percMeas{j}.percentiles = options.conf_level;
end

%% EVALUATION OF MODEL CONFIDENCE INTERVALS
% Loop: Realization of artificial data and parameter value
for i = 1:options.N
    for j = 1:length(D)
        % Evaluation of model densities for sampled parameter
        theta = parameters.MCMC.sample.par(:,randi(size(parameters.MCMC.sample.par,2)));
        Dens = getMixtureDensity(xf{j},theta,M,Mc);
        % Percentile intervals of measured cumulative distribution
        modelPred{j}.fec(:,i) = Dens.experiment(j).fec;
        modelPred{j}.fce(:,i) = Dens.experiment(j).fce;
        modelPred{j}.Fec(:,i) = Dens.experiment(j).Fec;
        modelPred{j}.Fce(:,i) = Dens.experiment(j).Fce;
    end
end
% Evaluation of percentile intervals
for j = 1:length(D)
    percModelPred{j}.fec = prctile(modelPred{j}.fec,options.conf_level,2);
    percModelPred{j}.fce = prctile(modelPred{j}.fce,options.conf_level,2);
    percModelPred{j}.Fec = prctile(modelPred{j}.Fec,options.conf_level,2);
    percModelPred{j}.Fce = prctile(modelPred{j}.Fce,options.conf_level,2);
end
    