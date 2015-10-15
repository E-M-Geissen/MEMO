%% LOAD OPTIMIZATION RESULTS AND SAMPLING RESULTS
load('optimization.mat');

%% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)
options.plot_subpop = 'true';
[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D,options);
saveas(M.fh.model_data,'./figs/_model_datafit');

