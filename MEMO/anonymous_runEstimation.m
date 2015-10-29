% exemplary run_estimations file, that uses most of the functionality
% provided by MEMO
% just replace "your_model_file" with the correct file name

clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'true'; % plots are shown
options.compute_P = 'true'; % computation of profiles enabled

%% Model
your_model_file; % load your personal model

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

% %% OPTIMIZATION
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


%% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)
if strcmp(plot_opt,'true')
    [M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D);
end


%% ANALYZE DISTRIBUTION OF LOG-LIKELIHOODS (FULL MODEL)
% Options
options_fit_ana.N = 10000;
options_fit_ana.bins = 25;
options_fit_ana.plot = plot_opt;
options_fit_ana.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');

% Analyze distribution of log-likelihoods
[parameters.fit_analysis.logL_sample,M.fh.fit_ana] = analyzeLogPostDistribution(parameters,M,Mc,D,options_fit_ana);

%% COMPUTE PROFILES (FULL MODEL)
if strcmp(options.compute_P,'true')
    % Options
    options_PL.plot = plot_opt;
    % Compute profile likelihoods
    [parameters,M.fh.PL] = computeProfiles(parameters,@(theta,opt) logLikelihoodMMc(theta,M,Mc,D,opt),options_PL);
end

%% PERFORM MCMC-SAMPLING (FULL MODEL) 
% Options
options_MCMC.nsimu_warmup = 1e5;
options_MCMC.nsimu_run    =3e5;

% Run MCMC-sampling
[parameters,M.fh.MCMC.logL_trace,M.fh.MCMC.par_trace,M.fh.MCMC.par_dis] = ...
    computeMCMCsample_DRAM(parameters,@(theta) logLikelihoodMMc(theta,M,Mc,D),options_MCMC);



% Comparison of mixture model and data (considering the uncertainties)
if strcmp(plot_opt,'true')
    options_DMplot.plot_model.subtype = 'MCMC';
    [M.fh.model_data_unc] = plotMixtureModel(parameters,M,Mc,D,options_DMplot);
end



%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');

% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);

% Print and plot optimal model
printModel(M_red,parameters_red)

if strcmp(plot_opt,'true')
    [M_red.fh.model_data] = plotMixtureModel(parameters_red,M_red,Mc_red,D);
end

%% ANALYZE DISTRIBUTION OF LOG-LIKELIHOODS (REDUCED MODEL)
% Options
options_fit_ana.N = 10000;
options_fit_ana.bins = 25;
options_fit_ana.plot = plot_opt;
options_fit_ana.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');

% Analyze distribution of log-likelihoods
[parameters_red.fit_analysis.logL_sample,M_red.fh.fit_ana] = analyzeLogPostDistribution(parameters_red,M_red,Mc_red,D,options_fit_ana);


%% COMPUTE PROFILES (REDUCED MODEL)

options_PL.plot = 'true';
if strcmp(options.compute_P,'true')
    [parameters_red,M_red.fh.PL] = computeProfiles(parameters_red,@(theta,opt) logLikelihoodMMc(theta,M_red,Mc_red,D,opt),options_PL);
end


%% PERFORM MCMC-SAMPLING (REDUCED MODEL) 
% Options
options_MCMC.nsimu_warmup = 1e5;
options_MCMC.nsimu_run    =3e5;

% Run model selection
[parameters_red,M_red.fh.MCMC.logL_trace,M_red.fh.MCMC.par_trace,M_red.fh.MCMC.par_dis] = ...
    computeMCMCsample_DRAM(parameters_red,@(theta) logLikelihoodMMc(theta,M_red,Mc_red,D),options_MCMC);


% Comparison of mixture model and data (considering the uncertainties)
if strcmp(plot_opt,'true')
    options_DMplot.plot_model.subtype = 'MCMC';
    [M_red.fh.model_data_unc] = plotMixtureModel(parameters_red,M_red,Mc_red,D,options_DMplot);
end




    
 