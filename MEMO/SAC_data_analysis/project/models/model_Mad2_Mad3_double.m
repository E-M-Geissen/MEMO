% Data:
%   Mad2 Mad3 double perturbations
%
% Model:
%   model can only be executed if additional parameter are loaded from
%   parameters_WT.mat
%   log-normal distribution for event data
%   Johnson SU distribution for censoring times
%   one component fixed across experiments = Wild type distribution
%   one component variable across experiments
 
%% SPECIFY MODEL PARAMETERS
syms mu_M2_65_P188_M3_30   esigma_M2_65_P188_M3_30     w_M2_65_P188_M3_30    ...
     mu_M2_65_P188_M3_60   esigma_M2_65_P188_M3_60     w_M2_65_P188_M3_60    ...
     mu_M2_65_P188_M3_120  esigma_M2_65_P188_M3_120    w_M2_65_P188_M3_120 ;

parameters.sym  = [mu_M2_65_P188_M3_30;   esigma_M2_65_P188_M3_30;     w_M2_65_P188_M3_30;    ...
                   mu_M2_65_P188_M3_60;   esigma_M2_65_P188_M3_60;     w_M2_65_P188_M3_60;    ...
                   mu_M2_65_P188_M3_120;  esigma_M2_65_P188_M3_120;    w_M2_65_P188_M3_120];


parameters.name = {'mu_{M2,65_P188-M3,30}'   ; 'esigma_{M2,65_P188-M3,30}'   ; 'w_{M2,M2,65_P188-M3,30}'   ; ...
                   'mu_{M2,65_P188-M3,60}'   ; 'esigma_{M2,65_P188-M3,60}'   ; 'w_{M2,M2,65_P188-M3,60}'   ; ...
                   'mu_{M2,65_P188-M3,120}'   ; 'esigma_{M2,65_P188-M3,120}'   ; 'w_{M2,M2,65_P188-M3,120}'  } ;

parameters.min  = [ log(5); log(1e-2); 0; ...
                    log(5); log(1e-2); 0; ...
                    log(5); log(1e-2); 0];

parameters.max  = [ log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1];

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);


%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'Johnson SU';
Mc.label.x = 'time [min]';
Mc.label.y = 'probability density';


% EXPERIMENT
i = 1;
data_65_Mad2_P188_30_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P188_M3_30,w_M2_65_P188_M3_30};
M.experiment(i).mu    = {mu_WT,mu_M2_65_P188_M3_30};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_65_P188_M3_30)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_WTen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i+1;
data_65_Mad2_P188_60_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P188_M3_60,w_M2_65_P188_M3_60};
M.experiment(i).mu    = {mu_WT,mu_M2_65_P188_M3_60};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_65_P188_M3_60)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_WTen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


% EXPERIMENT
i = i+1;
data_65_Mad2_P188_120_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P188_M3_120,w_M2_65_P188_M3_120};
M.experiment(i).mu    = {mu_WT,mu_M2_65_P188_M3_120};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_65_P188_M3_120)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_WTen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
Mc = getMixtureModel(Mc,parameters.sym);