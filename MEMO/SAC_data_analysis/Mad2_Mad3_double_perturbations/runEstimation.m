% estimate parameters of double perturbation datasets given fixed values
% for wild type and censoring distributions
clc;
clear all;
close all;

%% OPTIONS


%% LOAD PARAMETERS OF WT AND CENSORING DISTRIBUTION
% parameters of wild type distribution and censoring distribution derived from simultaneous analysis 
% of all single perturbation datasets with when WT fraction is a hill 
% function of the input 
load parameters_WT.mat

mu_WT = parameters_WT.MS.MAP.par(1);
esigma_WT=parameters_WT.MS.MAP.par(2);

gamma_cen = parameters_WT.MS.MAP.par(31);
esigma_WTen  =   parameters_WT.MS.MAP.par(32);
elambda_cen   =  parameters_WT.MS.MAP.par(33);
exi_cen=  parameters_WT.MS.MAP.par(34);

%% MODEL
model_Mad2_Mad3_double;

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
% Print result of estimation
printModel(M,parameters);
printModel(Mc,parameters);

% save
save('optimization','parameters','M','Mc','D','-v7.3');


%% PERFORM MCMC-SAMPLING (FULL MODEL)
% needed for uncertainty in wild type fractions
% Options
options_MCMC.nsimu_warmup = 3e3;
options_MCMC.nsimu_run    =1e5;

% Run MCMC-sampling
[parameters,M.fh.MCMC.logL_trace,M.fh.MCMC.par_trace,M.fh.MCMC.par_dis] = ...
    computeMCMCsample_DRAM(parameters,@(theta) logLikelihoodMMc(theta,M,Mc,D),options_MCMC);

% calculate information criterions for MLE
parameters= eval_performance(D,parameters);

% save results
save('mcmc_full','parameters','-v7.3');

%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);
% Print and plot optimal model
printModel(M,parameters);
printModel(M_red,parameters_red)


% save results
save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');


