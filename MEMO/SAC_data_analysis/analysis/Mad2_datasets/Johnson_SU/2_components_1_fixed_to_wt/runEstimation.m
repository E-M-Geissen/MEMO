clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'false'; % plots are shown

%% Model
model_Mad2_all_fv_JohnSU; % only Mad2 delta, discrete interval: 5min

% Print model on screen
printModel(M);
printModel(Mc);

% Process data for faster computations
D = processData(D);

%% OPTIMIZATION
% Options
options_fit.n_starts = 5000;
options_fit.plot = plot_opt;
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

% save
save('optimization','parameters','M','Mc','D','-v7.3');



%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);
% Print and plot optimal model
printModel(S{1,1}.model, S{1,1}.parameters)
printModel(M_red,parameters_red)


save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');


