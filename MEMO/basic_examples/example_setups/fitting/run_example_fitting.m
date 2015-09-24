% Loads a model, generates data for given parameters, estimates parameters
% from generated data, plots model data fit
%
% MODEL 
% =======
% three experimental conditions, linked by a common censoring distriution
% condition 1: 1 subpopulation
% condition 2: 2 subpopulations
% condition 3: 3 subpopulations
%
% DATA
% =======
% Data is generated from given true paramters


clear all

%% load model
example_model_fitting

printModel(M);
printModel(Mc);

%% set true parameters 
theta= [3; log(1); ...
        2; log(1); 6; log(1); 0.3;...
        2; log(1); 6; log(1); 0.3; 6.5; log(0.1); 0.5;...
        7; log(0.1)];

%% GENERATE ARTIFICAL DATA
% NUMBER OF DATA POINTS IN EXPERIMENT
n(1)=1000;
n(2)=1000;
n(3)=1000;

% SAMPLE FROM EVENT AND CENSORING MODEL
options_sampling.min = D{1}.observation_interval;;
% Event
Se = getMixtureSample(n,theta,M,options_sampling);
% Censoring
Sc = getMixtureSample(n,theta,Mc,options_sampling);

    
% (The size of the individual dataset will be conserved.)
% Loop: Experiments
for j = 1:length(D)
    % Assign times
    X  = Se{j}(1:n);
    Xc = Sc{j}(1:n);
    % Construct censored and uncensored dataset
    ind_X  = X  <= Xc;
    ind_Xc = Xc <= X;
    % Assign artificial data
    % (it has to be considered that events are only observed at discrete
    %  points in time, spaced according to D{2}.observation_interval.)
    D{j}.data.uncensored = D{j}.observation_interval * ceil(X(ind_X)/D{j}.observation_interval);
    D{j}.data.censored   = D{j}.observation_interval * floor(Xc(ind_Xc)/D{j}.observation_interval);
end


% Process data for faster computations
D = processData(D);

%% OPTIMIZATION
% Options
options_fit.n_starts = 15;
options_fit.plot = false;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...%'active-set',...
                           'display','off',...
                           'GradObj','on',...
                           'MaxIter',4000,...
                           'MaxFunEvals',4000*parameters.number);
% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMMc(theta,M,Mc,D,opt),options_fit);

% Print result of estimation
printModel(M,parameters);
printModel(Mc,parameters);

%% plot model data fit

options.plot_subpop = 'true';
options.plot_type = 'cdf';
[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D,options);

