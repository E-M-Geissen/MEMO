% Data:
% Mad3 delta 

% Model:
%   log-normal
%   data treatead as continuous (D{i}.observation_interval = 0;)
%   two components variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_M3_0_1 esigma_M3_0_1 mu_M3_0_2 esigma_M3_0_2 w_M3_0 ;

parameters.sym  = [ mu_M3_0_1 ; esigma_M3_0_1 ;mu_M3_0_2 ; esigma_M3_0_2 ; w_M3_0];

parameters.name = {'mu_1_{M3,0}' ; 'esigma_1_{M3,0}' ;'mu_2_{M3,0}' ; 'esigma_2_{M3,0}' ;'w_{M3,0}'};

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-2);log(5); log(1e-2);0];

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
% load data from file
data_Mad3_delta;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 0;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_0,w_M3_0};
M.experiment(i).mu    = {mu_M3_0_1,mu_M3_0_2};
M.experiment(i).sigma = {exp(esigma_M3_0_1),exp(esigma_M3_0_2)};

% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
