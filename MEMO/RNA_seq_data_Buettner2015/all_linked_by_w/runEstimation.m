%clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'true'; % plots are shown


%% Model
model_all_un_data_cen; % data of all four genes linked by one weighting parameter
Mc=[];

%%
% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

%% OPTIMIZATION
% Options
options_fit.n_starts = 200;
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


% save
save('optimization','parameters','M','D','-v7.3');


[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D);

%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');

% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);

% Print optimal full model and best reduced model to compare BICs
printModel(M,parameters);
printModel(M_red,parameters_red)


save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');