%% LOAD OPTIMIZATION RESULTS
load('optimization.mat')

%% COMPARISION OF MIXTURE MODEL AND DATA (FULL MODEL)
[M.fh.model_data] = plotMixtureModel(parameters,M,Mc,D);

savefig('./figs/model_data_fit.fig');


