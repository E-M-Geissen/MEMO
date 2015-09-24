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
fprintf('\n')
disp('Subpopulation structure of ''true'' model:');
example_model_selection_true_model
printModel(M);

%% set true parameters 
theta= [3; log(0.2); ...
        2; log(0.2); 4; log(0.2); 0.3;...
        2; log(0.2); 4; log(0.2); 0.1];

%% GENERATE ARTIFICAL DATA
fprintf('\n')
disp('Generating data from true model...');
% NUMBER OF DATA POINTS IN EXPERIMENT
n(1)=100;
n(2)=100;
n(3)=100;

% SAMPLE FROM EVENT AND CENSORING MODEL
options_sampling.min = D{1}.observation_interval;
% Event
Se = getMixtureSample(n,theta,M,options_sampling);


    
% (The size of the individual dataset will be conserved.)
% Loop: Experiments
for j = 1:length(D)
    % Assign times
    X  = Se{j}(1:n);
    
    % Construct censored and uncensored dataset
   
    % Assign artificial data
    % (it has to be considered that events are only observed at discrete
    %  points in time, spaced according to D{2}.observation_interval.)
    D{j}.data.uncensored = D{j}.observation_interval * ceil(X/D{j}.observation_interval);
    D{j}.data.censored =[];
end


% Process data for faster computations
D = processData(D);


%% load model with maximum number of subpopulations for backward model selection
clear M
clear parameters
fprintf('\n')
disp('Loading model with maximum number of subpopulations for successive reduction ...');
example_model_selection
printModel(M);
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
fprintf('\n')                       
disp('Find MLE for full model...');

% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMMc(theta,M,Mc,D,opt),options_fit);

% calculate information criterions
parameters= eval_performance(D,parameters);

% Print result of estimation
printModel(M,parameters);


%% plot model data fit

options.plot_subpop = 'true';
options.plot_type = 'cdf';
[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D,options);

%% MODEL SELECTION


% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');

fprintf('\n')
disp('Perform model selection ...');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);

% Print 
fprintf('\n')
fprintf('\n')
fprintf('\n')
disp('Full model:');
printModel(M,parameters)
fprintf('\n')
disp('Best reduced model:');
printModel(M_red,parameters_red)

[M.fh.model_data] = plotMixtureModel(parameters_red,M_red,Mc_red,D,options);
