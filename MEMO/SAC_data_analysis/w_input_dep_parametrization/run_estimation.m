clc;
clear all;
close all;

%% OPTIONS


%% Model
% all datasets with reduced Mad abundance, w parameterized hill function of
% input
model_w_dependency_on_input

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

%% OPTIMIZATION
% Options
options_fit.n_starts = 1000;
options_fit.plot = plot_opt;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...%'active-set',...
                               'display','off',...
                               'GradObj','on',...
                               'MaxIter',4000,...
                               'MaxFunEvals',4000*parameters.number);
% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMMc(theta,M,Mc,D,opt),options_fit);

% calculate information criterions for MLE
parameters= eval_performance(D,parameters);

% Print result of estimation
printModel(M,parameters);
printModel(Mc,parameters);

% save results of optimization
save('optimization','parameters','M','Mc','D','-v7.3');

