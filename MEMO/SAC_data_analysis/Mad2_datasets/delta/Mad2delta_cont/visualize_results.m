%% LOAD OPTIMIZATION RESULTS AND SAMPLING RESULTS
load('optimization.mat');
load('model_sel.mat');
Mc=[];
%% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)
options.plot_subpop = 'true';
[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D,options);
saveas(M.fh.model_data,'./figs/full_model_datafit');


%% COMPARISION OF MIXTURE MODEL AND DATA (REDUCED MODEL)

[M.fh.model_data] = plotMixtureModel(parameters_red,M_red,Mc_red,D);
saveas(M.fh.model_data,'./figs/red_model_datafit');