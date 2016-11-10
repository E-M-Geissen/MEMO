clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'false';%'true'; % plots are shown

%% Model
model_Mad3_delta; % only Mad3 delta, discrete interval: 5min
Mc=[];

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

%% OPTIMIZATION
% Options
options_fit.n_starts = 100;
options_fit.plot = plot_opt;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...%'active-set',...
                           'display','off',...
                           'GradObj','on',...
                           'MaxIter',4000,...
                           'MaxFunEvals',4000*parameters.number);
% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_fit);

% calculate information criterions for MLE
parameters= eval_performance(D,parameters);


% Print result of estimation
printModel(M,parameters);


% save optimization results
save('optimization','parameters','M','D','-v7.3');

%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);

% Print full and best reduced model
printModel(M,parameters);
printModel(M_red,parameters_red)

% save result of model selection
save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');


