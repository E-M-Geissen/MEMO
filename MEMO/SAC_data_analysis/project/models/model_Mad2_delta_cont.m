% Data:
% Mad2 delta 

% Model:
%  log-normal  
%  two components variable across experiments
%  observation interval for model = 0 min
 
%% SPECIFY MODEL PARAMETERS 
syms mu_M2_0_1 log_M2_0_1 mu_M2_0_2 log_M2_0_2 w_M2_0 ;

parameters.sym  = [ mu_M2_0_1 ; log_M2_0_1 ;mu_M2_0_2 ; log_M2_0_2 ; w_M2_0];

parameters.name = {'\mu_1' ; 'log(\sigma_1)' ;'\mu_2' ; 'log(\sigma_2)' ;'w'};

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-5);log(5); log(1e-5);0];      

parameters.max  = [log(2e3); log(1e1);log(2e3); log(1e1);1];



%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'Johnson SU';
Mc.label.x = 'time [min]';
Mc.label.y = 'probability density';

% EXPERIMENT
i = 1;
data_Mad2_delta;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 0;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_0,w_M2_0};
M.experiment(i).mu    = {mu_M2_0_1,mu_M2_0_2};
M.experiment(i).sigma = {exp(log_M2_0_1),exp(log_M2_0_2)};




% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);

