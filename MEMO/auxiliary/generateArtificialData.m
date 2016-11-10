% generateArtificialData generates an artificial dataset using M and theta,
%   which has the same properties as D, concerning its size and the 
%   censoring time distribution.
%
% USAGE:
% ======
% aD = generateArtificialData(theta,M,Mc,D)
%
% INPUTS:
% =======
% theta ... parameter value for model
% M ... mixture model for given data
% Mc ... mixture model for censoring times
% D ... measurement data
%
% Outputs:
% ========
% aD{i} ... artificial dataset for i-th experiment
%
% 2012/05/18 Jan Hasenauer

% function aD = generateArtificialData(theta,M,Mc,D)
function aD = generateArtificialData(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin == 4
    theta = varargin{1};
    M = varargin{2};
    Mc = varargin{3};
    D = varargin{4};
else
    error('generateArtificialData requires three inputs.')
end

%% INITIALIZE DATASET
for j = 1:length(D)
    aD{j}.name = D{j}.name;
    aD{j}.description = D{j}.description;
    aD{j}.observation_interval = D{j}.observation_interval;
end

%% CENSORING DISTRIBUTION
if isempty(Mc)
    % DETERMINE EMPIRICAL CENSORING DISTRIBUTION
    % -> This methods is actually wrong, as the observed distribution of
    %    censoring times is not the true distribution of censoring times,
    %    as it competes against the ever event.
    Xc_data = [];
    % Loop: Experiments
    for j = 1:length(D)
        % Collect censored data
        if ~isempty(D{j}.data.censored)
            Xc_data = [Xc_data; D{j}.data.censored];
        end
    end
    % Function handle for censoring times
    if ~isempty(Xc_data)
        Xc_fun = @(n) Xc_data(randi(length(Xc_data),[n,1]));
    else
        Xc_fun = @(n) inf(n,1);
    end
end

%% NUMBER OF DATA POINTS IN EXPERIMENTS
n = zeros(length(D),1);
for j = 1:length(D)
    if ~isempty(D{j}.data.uncensored)
        n(j) = n(j) + length(D{j}.data.uncensored);
    end
    if ~isempty(D{j}.data.censored)
        n(j) = n(j) + length(D{j}.data.censored);
    end
end

%% SAMPLE FROM EVENT AND CENSORING MODEL
options_sampling.min = D{1}.observation_interval;
% Event
Se = getMixtureSample(max(n),theta,M,options_sampling);
% Censoring
if ~isempty(Mc)
    Sc = getMixtureSample(max(n),theta,Mc,options_sampling);
else
    for j = 1:length(D)
        Sc{j} = Xc_fun(n(j));
    end
end
    
%% GENERATE ARTIFICAL DATA
% (The size of the individual dataset will be conserved.)
% Loop: Experiments
for j = 1:length(D)
    % Assign times
    X  = Se{j}(1:n(j));
    Xc = Sc{j}(1:n(j));
    % Construct censored and uncensored dataset
    ind_X  = X  <= Xc;
    ind_Xc = Xc <= X;
    % Assign artificial data
    % (it has to be considered that events are only observed at discrete
    %  points in time, spaced according to D{2}.observation_interval.)
    aD{j}.data.uncensored = D{j}.observation_interval * ceil(X(ind_X)/D{j}.observation_interval);
    aD{j}.data.censored   = Xc(ind_Xc);
end
