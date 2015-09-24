clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'true'; % plots are shown
options.compute_P = 'true'; % computation of profiles

%% Model
%model_kinetics_fun_mu; %
model_kinetics_fun_mu_new_param_log10_cen02; %

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

Mc=[];
%% OPTIMIZATION
% Options
options_fit.n_starts = 300;
options_fit.plot = plot_opt;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...%'active-set',...
                           'display','off',...
                           'GradObj','on',...
                           'MaxIter',4000,...
                           'MaxFunEvals',4000*parameters.number);
%

parameters.guess= [-1.31224933892433;-0.800170025950256;0.308686895511533;-1.00664990398034;-0.793472192038934;-0.958362694943307;-0.566965993149609;-0.880727397693913;-0.782505437418861;-0.844955577642641;-0.659188864221808;-0.893688851987427;-0.609383901991252;-0.805528955667771;-0.336775651990680;-0.904264595011374;-0.402351929567410;-0.849032552318242;-0.533867263234241;-0.924174324473709;-0.677188465471327;-1.01589979044468;-0.735719735078803;0.000396284244769350;2.03477266321198;-6.17071496092370;3.13002830013107;-4.52353825267214;-1.64103414327923;2.89605292538153];
%parameters.guess=[1+zeros(30,1)];
% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_fit);
% Print result of estimation
printModel(M,parameters);

%%
% save
save('optimization','parameters','M','Mc','D','-v7.3');
parameters_plot=parameters;
parameters_plot.MS.MAP.par=parameters.guess;
 [M.fh.model_data] = plotMixtureModel3(parameters_plot,M,Mc,D);
 printModel(M,parameters_plot);
% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)
%if strcmp(plot_opt,'true')
options.plot_type='cdf';
    [M.fh.model_data] = plotMixtureModel3(parameters,M,Mc,D);
%end
%print('-depsc2','-r1000',['./figs/full_model_datafit1']);

%% ANALYZE DISTRIBUTION OF LOG-LIKELIHOODS (FULL MODEL)
% % Options
% options_fit_ana.N = 10000;
% options_fit_ana.bins = 25;
% options_fit_ana.plot = plot_opt;
% options_fit_ana.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% % Analyze distribution of log-likelihoods
% [parameters.fit_analysis.logL_sample,M.fh.fit_ana] = analyzeLogPostDistribution(parameters,M,Mc,D,options_fit_ana);
% 
% 
% save('Ldis_full','parameters','-v7.3');
% print('-depsc2','-r1000',['./figs/full_loglike']);
% COMPUTE PROFILES (FULL MODEL)
if strcmp(options.compute_P,'true')
    % Options
    options_PL.plot = 'true';
options_PL.P_next_step.min = 1e-1;
    options_PL.parameter_index=[24:30];
    % Compute profile likelihoods
    [parameters,M.fh.PL] = computeProfiles(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_PL);
end
save('prof_full_2','parameters','-v7.3');

opt.interval = 'static';

plotP(parameters,[],[24:30],opt);
% print('-depsc2','-r1000',['./figs/profile_full']);
% %% PERFORM MCMC-SAMPLING (FULL MODEL) 
% % Options
% options_MCMC.nsimu_warmup = 3e3;
% options_MCMC.nsimu_run    =1e5;
% 
% % Run MCMC-sampling
% [parameters,M.fh.MCMC.logL_trace,M.fh.MCMC.par_trace,M.fh.MCMC.par_dis] = ...
%     computeMCMCsample_DRAM(parameters,@(theta) logLikelihoodMMc(theta,M,Mc,D),options_MCMC);
% 
% save('mcmc_full','parameters','-v7.3');
% 
% % Comparison of mixture model and data (considering the uncertainties)
% if strcmp(plot_opt,'true')
%     options_DMplot.plot_model.subtype = 'MCMC';
%     [M.fh.model_data_unc] = plotMixtureModel3(parameters,M,Mc,D,options_DMplot);
% end



%% MODEL SELECTION

% Options
options_modsel.criterion = 'BIC';
options_modsel.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% Run model selection
[M_red,Mc_red,parameters_red,S,R] = performModelSelection(parameters,M,Mc,D,options_modsel);
% Print and plot optimal model
printModel(S{1,1}.model, S{1,1}.parameters)
printModel(M_red,parameters_red)

%if strcmp(plot_opt,'true')
    [M_red.fh.model_data] = plotMixtureModel3(parameters_red,M_red,Mc_red,D);
%end


save('model_sel','M_red','Mc_red','parameters_red','S','R','-v7.3');
%print('-depsc2','-r1000',['./figs/red_model_datafit']);

%%
printModel(S{1,1}.model, S{1,1}.parameters)
printModel(M_red,parameters_red)
printModel(Mc_red, parameters_red)
% %% ANALYZE DISTRIBUTION OF LOG-LIKELIHOODS (REDUCED MODEL)
% % Options
% options_fit_ana.N = 10000;
% options_fit_ana.bins = 25;
% options_fit_ana.plot = plot_opt;
% options_fit_ana.optim_opt.fmincon = optimset(options_fit.fmincon,'algorithm','active-set');
% % Analyze distribution of log-likelihoods
% [parameters_red.fit_analysis.logL_sample,M_red.fh.fit_ana] = analyzeLogPostDistribution(parameters_red,M_red,Mc_red,D,options_fit_ana);
% 
% 
% save('Ldis_red','parameters_red','-v7.3');
% %print('-depsc2','-r1000',['./figs/red_loglike']);
% %% COMPUTE PROFILES (REDUCED MODEL)
% 
% options_PL.plot = 'true';
% if strcmp(options.compute_P,'true')
%     [parameters_red,M_red.fh.PL] = computeProfiles(parameters_red,@(theta,opt) logLikelihoodMMc(theta,M_red,Mc_red,D,opt),options_PL);
% end
% 
% save('prof_red','parameters_red','-v7.3');
% print('-depsc2','-r1000',['./figs/profile_red']);
% %% PERFORM MCMC-SAMPLING (REDUCED MODEL) 
%  % Options
%  options_MCMC.nsimu_warmup = 5e3;
%  options_MCMC.nsimu_run    =3e5;
% 
% % Run model selection
% [parameters_red,M_red.fh.MCMC.logL_trace,M_red.fh.MCMC.par_trace,M_red.fh.MCMC.par_dis] = ...
%     computeMCMCsample_DRAM(parameters_red,@(theta) logLikelihoodMMc(theta,M_red,Mc_red,D),options_MCMC);
% 
% save('mcmc_red','parameters_red','-v7.3');
% 
% print(M_red.fh.MCMC.par_dis,'-depsc2','-r1000',['./figs/MCMC_red_pardis']);
% print(M_red.fh.MCMC.par_trace,'-depsc2','-r1000',['./figs/MCMC_red_partrace']);
% print(M_red.fh.MCMC.logL_trace,'-depsc2','-r1000',['./figs/MCMC_red_logLtrace']);
% 
% % Comparison of mixture model and data (considering the uncertainties)
% if strcmp(plot_opt,'true')
%     options_DMplot.plot_model.subtype = 'MCMC';
%     [M_red.fh.model_data_unc] = plotMixtureModel33(parameters_red,M_red,Mc_red,D,options_DMplot);
% end
% %print('-depsc2','-r1000',['./figs/red_model_data_unc']);
%  set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 20 16])
% print('-depsc',['./figs/red_model_data_unc']);
% %%
% printModel(S{1,1}.model, S{1,1}.parameters)
% printModel(M_red,parameters_red)
% 
%   % Parameter histogram
%   options.interval = 'static';
%     par_dis=plotMCMC(parameters_red)
%     print(par_dis,'-depsc2','-r1000',['./figs/MCMC_red_pardis']);
%     