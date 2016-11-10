% file used for the generation of data used in file run_example_model_selection.m
% caution data is not written to files

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
