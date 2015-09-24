% Data:
%   control
%   all Ma32 data (no double mutants)
% Model:
%   Johnson SU
%   one component fixed across experiments = wild type
%   one component variable across experiments
 
%% SPECIFY MODEL PARAMETERS 
syms  gamma_wt       log_sigma_wt         log_lambda_wt         log_xi_wt              ...
      gamma_M3_120   log_sigma_M3_120     log_lambda_M3_120     log_xi_M3_120    w_M3_120  ...
      gamma_M3_60    log_sigma_M3_60      log_lambda_M3_60      log_xi_M3_60     w_M3_60  ...
      gamma_M3_30    log_sigma_M3_30      log_lambda_M3_30      log_xi_M3_30     w_M3_30  ...
      gamma_M3_0     log_sigma_M3_0       log_lambda_M3_0       log_xi_M3_0      w_M3_0  ...
      gamma_cen      log_sigma_cen        log_lambda_cen        log_xi_cen;
 
 
parameters.sym  = [gamma_wt       log_sigma_wt         log_lambda_wt         log_xi_wt              ...
                   gamma_M3_120   log_sigma_M3_120     log_lambda_M3_120     log_xi_M3_120    w_M3_120  ...
                   gamma_M3_60    log_sigma_M3_60      log_lambda_M3_60      log_xi_M3_60     w_M3_60  ...
                   gamma_M3_30    log_sigma_M3_30      log_lambda_M3_30      log_xi_M3_30     w_M3_30  ...
                   gamma_M3_0     log_sigma_M3_0       log_lambda_M3_0       log_xi_M3_0      w_M3_0  ...
                   gamma_cen      log_sigma_cen        log_lambda_cen        log_xi_cen;];
 
 
parameters.name = {'\gamma_{wt}' ;     'log(\sigma_{wt})' ;      'log(\lambda_{wt})' ;     'log(\xi_{wt})'    ;                     ...
                   '\gamma_{M3,120}';  'log(\sigma_{M3,120})';   'log(\lambda_{M3,120})';  'log(\xi_{M3,120})'; 'w_{M3,120}'; ...
                   '\gamma_{M3,60}';   'log(\sigma_{M3,60})';    'log(\lambda_{M3,60})';   'log(\xi_{M3,60})';  'w_{M3,60}'; ...
                   '\gamma_{M3,30}';   'log(\sigma_{M3,30})';    'log(\lambda_{M3,30})';   'log(\xi_{M3,30})';  'w_{M3,30}'; ...
                   '\gamma_{M3,0}';    'log(\sigma_{M3,0})';     'log(\lambda_{M3,0})';    'log(\xi_{M3,0})';   'w_{M3,0}'; ...
                   '\gamma_{cen}'  ;   'log(\sigma_{cen})' ;     'log(\lambda_{cen})' ;    'log(\xi_{cen)})' };
               
               
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  -20   ; log(1e-4); log(5); log(5);  ...
                      -20   ; log(1e-4); log(5); log(5);0; ...       
                     -20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5);0; ...
                     -20   ; log(1e-4); log(5); log(5)];  
                 
parameters.max  = [ 20; log(1e4); log(2e3); log(2e3);    ...
                    20; log(1e4); log(2e3); log(2e3);1; ...
                    20; log(1e4); log(2e3); log(2e3);1; ...
                    20; log(1e4); log(2e3); log(2e3); 1;...
                    20; log(1e4); log(2e3); log(2e3);1; ...
                    20; log(1e4); log(2e3); log(2e3)];
                      


%% SPECIFY MODEL AND DATA
M.mixture.type = 'Johnson SU';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'Johnson SU';
Mc.label.x = 'time [min]';
Mc.label.y = 'probability density';





% EXPERIMENT
i = 1;
data_120_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_120,w_M3_120};
M.experiment(i).gamma  = { gamma_wt, gamma_M3_120 };
M.experiment(i).sigma  = {  exp(log_sigma_wt),exp(log_sigma_M3_120)};
M.experiment(i).lambda = {exp(log_lambda_wt),exp(log_lambda_M3_120)};
M.experiment(i).xi     = {exp(log_xi_wt),exp(log_xi_M3_120)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};



% EXPERIMENT
i = i + 1;
data_60_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_60,w_M3_60};
M.experiment(i).gamma  = { gamma_wt  , gamma_M3_60 };
M.experiment(i).sigma  = { exp(log_sigma_wt),exp(log_sigma_M3_60) };
M.experiment(i).lambda = {exp(log_lambda_wt),exp(log_lambda_M3_60) };
M.experiment(i).xi     = {exp(log_xi_wt),exp(log_xi_M3_60)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_30,w_M3_30};
M.experiment(i).gamma  = { gamma_wt ,gamma_M3_30};
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M3_30)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M3_30) };
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M3_30)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_Mad3_delta;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_0,w_M3_0};
M.experiment(i).gamma  = {gamma_wt ,  gamma_M3_0 };
M.experiment(i).sigma  = { exp(log_sigma_M3_0), exp(log_sigma_wt)};
M.experiment(i).lambda = {exp(log_lambda_M3_0),exp(log_lambda_wt)};
M.experiment(i).xi     = {exp(log_xi_M3_0), exp(log_xi_wt)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i+1;
data_WT_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).gamma  = {      gamma_wt };
M.experiment(i).sigma  = { exp(log_sigma_wt)};
M.experiment(i).lambda = {exp(log_lambda_wt)};
M.experiment(i).xi     = {    exp(log_xi_wt)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};


% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
Mc = getMixtureModel(Mc,parameters.sym);