% Data:
%   
%   all Mad2 data (no double mutants)
% Model:
%   mixture of two Johnson SU distribtions
%   one component fixed across experiments = wild type
%   one component variable across experiments = strain specific
 
%% SPECIFY MODEL PARAMETERS 
syms gamma_wt         log_sigma_wt           log_lambda_wt          log_xi_wt              ...
     gamma_M2_200     log_sigma_M2_200       log_lambda_M2_200      log_xi_M2_200        w_M2_200  ...
     gamma_M2_80      log_sigma_M2_80        log_lambda_M2_80       log_xi_M2_80         w_M2_80  ...
     gamma_M2_65_P50  log_sigma_M2_65_P50    log_lambda_M2_65_P50   log_xi_M2_65_P50     w_M2_65_P50  ...
     gamma_M2_65_P188 log_sigma_M2_65_P188   log_lambda_M2_65_P188  log_xi_M2_65_P188    w_M2_65_P188  ...
     gamma_M2_40      log_sigma_M2_40        log_lambda_M2_40       log_xi_M2_40         w_M2_40  ...
     gamma_M2_20      log_sigma_M2_20        log_lambda_M2_20       log_xi_M2_20         w_M2_20  ...
     gamma_M2_10      log_sigma_M2_10        log_lambda_M2_10       log_xi_M2_10         w_M2_10  ...
     gamma_M2_0       log_sigma_M2_0         log_lambda_M2_0        log_xi_M2_0          w_M2_0  ...
     gamma_cen        log_sigma_cen          log_lambda_cen         log_xi_cen;
 
 
parameters.sym  = [ gamma_wt ;         log_sigma_wt ;        log_lambda_wt ;         log_xi_wt    ;          ...
                    gamma_M2_200;      log_sigma_M2_200  ;   log_lambda_M2_200    ;  log_xi_M2_200  ;      w_M2_200 ; ...
                    gamma_M2_80;       log_sigma_M2_80  ;    log_lambda_M2_80 ;      log_xi_M2_80;         w_M2_80 ; ...
                    gamma_M2_65_P50  ; log_sigma_M2_65_P50;  log_lambda_M2_65_P50 ;  log_xi_M2_65_P50   ;  w_M2_65_P50 ; ...
                    gamma_M2_65_P188 ; log_sigma_M2_65_P188; log_lambda_M2_65_P188;  log_xi_M2_65_P188 ;   w_M2_65_P188;  ...
                    gamma_M2_40 ;      log_sigma_M2_40  ;    log_lambda_M2_40  ;     log_xi_M2_40 ;        w_M2_40;  ...
                    gamma_M2_20;       log_sigma_M2_20 ;     log_lambda_M2_20  ;     log_xi_M2_20  ;       w_M2_20 ; ...
                    gamma_M2_10 ;      log_sigma_M2_10  ;    log_lambda_M2_10  ;     log_xi_M2_10  ;       w_M2_10 ; ...
                    gamma_M2_0  ;      log_sigma_M2_0 ;      log_lambda_M2_0   ;     log_xi_M2_0    ;      w_M2_0 ; ...
                    gamma_cen ;        log_sigma_cen ;       log_lambda_cen    ;     log_xi_cen];
 
 
parameters.name = {'\gamma_{wt}'         ; 'log(\sigma_{wt})' ;        'log(\lambda_{wt})' ;       'log(\xi_{wt})'    ;                     ...
                   '\gamma_{M2,200}'     ; 'log(\sigma_{M2,200})';     'loglambda_{M2,200})';      'log(\xi_{M2,200})';      'w_{M2,200}'; ...
                   '\gamma_{M2,80}'      ; 'log(\sigma_{M2,80})';      'log(\lambda_{M2,80})';     'log(\xi__{M2,80})';      'w_{M2,80}'; ...
                   '\gamma_{M2,65_P50}'  ; 'log(\sigma_{M2,65_P50})';  'log(\lambda_{M2,65_P50})'; 'log(\xi__{M2,65_P50})';  'w_{M2,65_P50}'; ...
                   '\gamma_{M2,65_P188}' ; 'log(\sigma_{M2,65_P188})'; 'log(\lambda_{M2,65_P188})';'log(\xi__{M2,65_P188})'; 'w_{M2,65_P188}'; ...
                   '\gamma_{M2,40}'      ; 'log(\sigma_{M2,40})';      'log(\lambda_{M2,40})';     'log(\xi__{M2,40})';      'w_{M2,40}'; ...
                   '\gamma_{M2,20}'      ; 'log(\sigma_{M2,20})';      'log(\lambda_{M2,20})';     'log(\xi__{M2,20})';      'w_{M2,20}'; ...
                   '\gamma_{M2,10}'      ; 'log(\sigma_{M2,10})';      'log(\lambda_{M2,10})';     'log(\xi__{M2,10})';      'w_{M2,10}'; ...
                   '\gamma_{M2,0}'       ; 'log(\sigma_{M2,0})';       'log(\lambda_{M2,0})';      'log(\xi__{M2,0})';       'w_{M2,0}'; ...
                   '\gamma_{cen}'        ; 'log(\sigma_{cen})'     ;   'log(\lambda_{cen})' ;      'log(\xi_{cen})' };
               
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  -20   ; log(1e-4); log(5); log(5);  ...
                     -20   ; log(1e-4); log(5); log(5); 0; ...     
                     -20   ; log(1e-4); log(5); log(5); 0; ...     
                     -20   ; log(1e-4); log(5); log(5); 0; ...    
                     -20   ; log(1e-4); log(5); log(5); 0; ...       
                     -20   ; log(1e-4); log(5); log(5); 0; ...       
                     -20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5)];  
                 
parameters.max  = [ 20; log(1e4); log(2e3); log(2e3);    ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 1; ...
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
data_200_Mad2_P259bp;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_200,w_M2_200};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_200 };
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_200)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_200)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_200)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_80_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_80,w_M2_80};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_80 };
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_80)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_80)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_80)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i+1 ;
data_65_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P50,w_M2_65_P50};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_65_P50 };
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_65_P50)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_65_P50)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_65_P50)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_65_P188_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P188,w_M2_65_P188};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_65_P188};
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_65_P188)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_65_P188)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_65_P188)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_40_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_40,w_M2_40};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_40 };
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_40)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_40)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_40)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};



% EXPERIMENT
i = i + 1;
data_20_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_20,w_M2_20};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_20};
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_20)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_20)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_20)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_10_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_10,w_M2_10};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_10};
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_10)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_10)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_10)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i + 1;
data_Mad2_delta;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_0,w_M2_0};
M.experiment(i).gamma  = { gamma_wt,gamma_M2_0};
M.experiment(i).sigma  = { exp(log_sigma_wt), exp(log_sigma_M2_0)};
M.experiment(i).lambda = {exp(log_lambda_wt), exp(log_lambda_M2_0)};
M.experiment(i).xi     = {exp(log_xi_wt), exp(log_xi_M2_0)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i = i+1;
data_WT_Mad2;

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