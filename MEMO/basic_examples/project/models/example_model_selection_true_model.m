%% SPECIFY MODEL PARAMETERS 
syms mu_1_exp1     log_sigma_1_exp1          ...
     mu_1_exp2     log_sigma_1_exp2    mu_2_exp2     log_sigma_2_exp2    w_exp2      ...
     mu_1_exp3     log_sigma_1_exp3    mu_2_exp3     log_sigma_2_exp3    w_exp3;
 
parameters.sym  = [mu_1_exp1;     log_sigma_1_exp1;         ...
                   mu_1_exp2;     log_sigma_1_exp2;    mu_2_exp2;     log_sigma_2_exp2;    w_exp2;      ...
                   mu_1_exp3;     log_sigma_1_exp3;    mu_2_exp3;     log_sigma_2_exp3;    w_exp3];
               
parameters.name = {'\mu_1_{exp1}';     'log(\sigma_1_{exp1})';     ...
                   '\mu_1_{exp2}';     'log(\sigma_1_{exp2})';    'mu_2_{exp2}';     'log(\sigma_2_{exp2})';    'w_{exp2}';  ...
                   '\mu_1_{exp3}';     'log(\sigma_1_{exp3})';    'mu_2_{exp3}';     'log(\sigma_2_{exp3})';    'w_{exp3}'}; 
               
parameters.number = length(parameters.sym);

parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-1); ...  
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); log(5); log(1e-1); 0 ];       
                     
parameters.max  = [ log(2e3); log(1e1);  ...
                    log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); log(2e3); log(1e1); 1 ];
                
parameters.guess = parameters.max - parameters.min;

%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';


% EXPERIMENT
i = 1;
% exp_cond_1; % experimental condition without subpopulations 
% 
D{i}.name = 'condition 1';
% D{i}.description = [];
% D{i}.data.uncensored = Tm;
% D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = 'condition 1';
M.experiment(i).size   = 1;
M.experiment(i).w      = {1};
M.experiment(i).mu     = {mu_1_exp1};
M.experiment(i).sigma  = {exp(log_sigma_1_exp1)};



% % EXPERIMENT
 i = i + 1;
% exp_cond_2; % experimental condition with 2 subpopulations 

D{i}.name = 'condition 2';
D{i}.description = [];
% D{i}.data.uncensored = Tm;
% D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = 'condition 2';
M.experiment(i).size   = 2;
M.experiment(i).w      = {w_exp2, 1-w_exp2};
M.experiment(i).mu     = {mu_1_exp2, mu_2_exp2};
M.experiment(i).sigma  = {exp(log_sigma_1_exp2),exp(log_sigma_2_exp2)};



% EXPERIMENT
i = i + 1;
% exp_cond_3; % experimental condition with 3 subpopulations 
% 
D{i}.name = 'condition 3';
D{i}.description = [];
% D{i}.data.uncensored = Tm;
% D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = 'condition 3';
M.experiment(i).size   = 2;
M.experiment(i).w      = {w_exp3, 1-w_exp3};
M.experiment(i).mu     = {mu_1_exp3, mu_2_exp3};
M.experiment(i).sigma  = {exp(log_sigma_1_exp3),exp(log_sigma_2_exp3)};


% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M ,parameters.sym);
Mc = [];
