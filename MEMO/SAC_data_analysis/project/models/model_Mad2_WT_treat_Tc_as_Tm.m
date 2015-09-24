% Data:
% Mad2 wild type, right censored measurements treated as prometaphase lengths 
%
% Model:
%  log-normal 
%  two components variable across experiments
% no censoring model because no right censored data included in dataset

%% SPECIFY MODEL PARAMETERS 
syms mu_M3_WT_1 log_sigma_M3_WT_1 mu_M3_WT_2 log_sigma_M3_WT_2 w_M3_WT ;

parameters.sym  = [ mu_M3_WT_1 ; log_sigma_M3_WT_1 ;mu_M3_WT_2 ; log_sigma_M3_WT_2 ; w_M3_WT];

parameters.name = {'mu_1_{M3,WT}' ; 'log(\sigma_1_{M3,WT})' ;'mu_1_{M3,WT}' ; 'log(\sigma_1_{M3,WT})' ;'w_{M3,WT}'  };

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-2);log(5); log(1e-2);0 ];  


parameters.max  = [ log(2e3); log(1e1);log(2e3); log(1e1);1];


%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';



% EXPERIMENT
i = 1;

% load data from file
data_WT_Mad2_treat_Tc_as_Tm;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_WT,w_M3_WT};
M.experiment(i).mu    = {mu_M3_WT_1,mu_M3_WT_2};
M.experiment(i).sigma = {exp(log_sigma_M3_WT_1),exp(log_sigma_M3_WT_2)};




% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
