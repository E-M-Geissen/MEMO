%% model illustration MEMO model setup for a model with log-normally distributed events of interest and Johnson SU distributed censoring events. 
%%
%% SPECIFY MODEL PARAMETERS
syms    mu_exp1_sub1 mu_exp1_sub2 log_sigma_exp1_sub1 log_sigma_exp1_sub2 w_exp1_sub1  ... 
        gamma_cen  esigma_cen  elambda_cen  exi_cen;

parameters.sym  = [mu_exp1_sub1; mu_exp1_sub2; log_sigma_exp1_sub1; log_sigma_exp1_sub2; w_exp1_sub1;  ... 
										gamma_cen; esigma_cen; elambda_cen; exi_cen];

% parameter names for plotting and other vizualization
parameters.name = {'\mu_{exp1,sub1}'; 'log(\sigma_{exp1,sub1})'; '\mu_{exp1,sub2}'; 'log(\sigma_{exp1,sub2})'; 'w_{exp1,sub1}'; ...
									 '\gamma_{cen}'; 'log(\sigma_{cen})'; 'log(\lambda_{cen})'; 'log(\xi_{cen})' };

parameters.number = length(parameters.sym);

parameters.min  = [ log(5); log(1e-1); 0; ...
										-20   ; log(1e-4); log(5); log(5) ];

parameters.max  = [  log(2e3); log(1e1); 1; ...
										20; log(1e4); log(2e3); log(2e3)];

%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal'; % specification of distribution type 
M.label.x = 'time [min]'; % x-label for model-data-comparison plots, in units of measurement data 
M.label.y = 'probability density'; % y-label for model-data-comparison plots

%% Experimental condition No. 1
% File were data for Exp.1 is stored

i = 1 ;

data_file_exp_1;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;  % number of mixture components=subpopulations
M.experiment(i).w     = {1-w_exp1_sub1,w_exp1_sub1}; % weights
M.experiment(i).mu    = {mu_exp1_sub2,mu_exp1_sub1};  %  one §\mcommentfont $\mu$§ per subpopulations
M.experiment(i).sigma = {exp(log_sigma_exp1_sub2),exp(log_sigma_exp1_sub1)}; % one §\mcommentfont $\sigma$§ per subpopulations, exp(sigma) for parmeter estimation in logarithmic space


D{i}.name = tit;  % Titel defined in data file data_file_exp_1.m
D{i}.description = []; % optional
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5; % in this example data is interval censored and inter observation time is 5 min

%% Experimental condition No. 2 
%i = i+1 ;

% and so on


% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
%Mc = getMixtureModel(Mc,parameters.sym);
