%% LOAD OPTIMIZATION RESULTS AND SAMPLING RESULTS
load('optimization.mat');
load('model_sel.mat');
load('mcmc_red.mat')
load('prof_red')
%% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)
options.plot_subpop = 'true';
[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D,options);
saveas(M.fh.model_data,'./figs/full_model_datafit');


%% COMPARISION OF MIXTURE MODEL AND DATA (REDUCED MODEL)

[M.fh.model_data] = plotMixtureModel(parameters_red,M_red,Mc_red,D);
saveas(M.fh.model_data,'./figs/red_model_datafit');

%% COMPARISION OF MIXTURE MODEL AND DATA WITH UNCERTAINTY (REDUCED MODEL)
        options.plot_model.subtype = 'MCMC';
    [M_red.fh.model_data_unc] = plotMixtureModel(parameters_red,M_red,Mc_red,D,options);
    
%% PLOT PROFILE LIKELIHOODs FOR REDUCED MODEL
options.interval = 'static';% 'dynamic';
options.mark_constraint='false';
options.type= 'R';

plotP(parameters_red,[], [1:parameters_red.number], options);
