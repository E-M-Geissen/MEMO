%% LOAD OPTIMIZATION RESULTS AND SAMPLING RESULTS
load('optimization.mat')
load('mcmc_full.mat')
load('model_sel.mat')
%% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)

[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D);
savefig('./figs/full_model_datafit');


% Comparison of mixture model and data (considering the uncertainties)
options_DMplot.plot_model.subtype = 'MCMC';
[M.fh.model_data_unc] = plotMixtureModel(parameters,M,Mc,D,options_DMplot);
savefig('./figs/full_model_datafit_unc');

%% COMPARISION OF MIXTURE MODEL AND DATA (REDUCED MODEL)

[M.fh.model_data] = plotMixtureModel(parameters_red,M_red,Mc_red,D_red);
savefig('./figs/full_model_datafit');