clc;
clear all;
close all;

%% OPTIONS
plot_opt =  'false'; % plots are not shown
options.compute_P = 'true'; % computation of profiles

%% Model
model_Mad2_all_fv; % all Mad2 data, discrete interval: 5min

% Print model on screen
printModel(M);
printModel(Mc);
% Process data for faster computations
D = processData(D);
 
%% OPTIMIZATION
 % Options
options_fit.n_starts = 500;
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

[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D);

% save
save('optimization','parameters','M','Mc','D','-v7.3');




%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);


save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');

% Print and plot optimal model
printModel(M,parameters);
printModel(Mc_red,parameters_red);

[M_red.fh.model_data] = plotMixtureModel(parameters_red,M_red,Mc_red,D);

%% ANALYZE DISTRIBUTION OF LOG-LIKELIHOODS (REDUCED MODEL)
% Options
options_fit_ana.N = 10000;
options_fit_ana.bins = 25;
options_fit_ana.plot = plot_opt;
options_fit_ana.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Analyze distribution of log-likelihoods
[parameters_red.fit_analysis.logL_sample,M_red.fh.fit_ana] = analyzeLogPostDistribution(parameters_red,M_red,Mc_red,D,options_fit_ana);


save('Ldis_red','parameters_red','-v7.3');

%% COMPUTE PROFILES (REDUCED MODEL)

options_PL.plot = 'true';
if strcmp(options.compute_P,'true')
    [parameters_red,M_red.fh.PL] = computeProfiles(parameters_red,@(theta,opt) logLikelihoodMMc(theta,M_red,Mc_red,D,opt),options_PL);
end

save('prof_red','parameters_red','-v7.3');

%% PERFORM MCMC-SAMPLING (REDUCED MODEL) 
 % Options
 options_MCMC.nsimu_warmup = 1e5;
 options_MCMC.nsimu_run    =3e5;

% Run model selection
[parameters_red,M_red.fh.MCMC.logL_trace,M_red.fh.MCMC.par_trace,M_red.fh.MCMC.par_dis] = ...
    computeMCMCsample_DRAM(parameters_red,@(theta) logLikelihoodMMc(theta,M_red,Mc_red,D),options_MCMC);

save('mcmc_red','parameters_red','-v7.3');



