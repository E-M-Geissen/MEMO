current_path= pwd;
addpath(genpath([current_path, filesep, 'auxiliary']));
addpath(genpath([current_path,filesep,'SAC_data_analysis',filesep,'project']));
addpath(genpath([current_path,filesep,'NGF-Erk_data_analysis',filesep,'project']));
addpath(genpath([current_path,filesep,'basic_examples',filesep,'project']));

% for mcmc sampling the matlab toolbox DRAM (Haario2006) has to be in the
% matlab path
 % addpath(genpath('D:\modelling\matlab\Toolboxen\mcmcstat'));


% General properties
TextSizes.DefaultAxesFontSize = 12;
TextSizes.DefaultTextFontSize = 14;
set(0,TextSizes);