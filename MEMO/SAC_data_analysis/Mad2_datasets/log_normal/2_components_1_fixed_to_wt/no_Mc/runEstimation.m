clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'true'; % plots are shown


%% Model
model_Mad2_all_fv_noMc; % only Mad2 delta, discrete interval: 5min

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
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_fit);

% calculate information criterions for MLE
parameters= eval_performance(D,parameters);

% Print result of estimation
printModel(M,parameters);


% save
save('optimization','parameters','M','D','-v7.3');




%% MODEL SELECTION
Mc=[];
% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);
% Print and plot optimal model
printModel(M,parameters);
printModel(M_red,parameters_red)

%if strcmp(plot_opt,'true')
    [M_red.fh.model_data] = plotMixtureModel3(parameters_red,M_red,Mc_red,D);
%end


save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');

     