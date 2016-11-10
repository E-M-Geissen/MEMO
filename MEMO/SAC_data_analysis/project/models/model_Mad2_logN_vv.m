% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   two component variable across experiments, except wild type
 
%% SPECIFY MODEL PARAMETERS
syms mu_wt            log_sigma_wt ... 
     mu_M2_200_1      log_sigma_M2_200_1      mu_M2_200_2      log_sigma_M2_200_2      w_M2_200   ...
     mu_M2_80_1       log_sigma_M2_80_1       mu_M2_80_2       log_sigma_M2_80_2       w_M2_80   ...
     mu_M2_65_P50_1   log_sigma_M2_65_P50_1   mu_M2_65_P50_2   log_sigma_M2_65_P50_2   w_M2_65_P50   ...
     mu_M2_65_P188_1  log_sigma_M2_65_P188_1  mu_M2_65_P188_2  log_sigma_M2_65_P188_2  w_M2_65_P188   ...
     mu_M2_40_1       log_sigma_M2_40_1       mu_M2_40_2       log_sigma_M2_40_2       w_M2_40   ...
     mu_M2_20_1       log_sigma_M2_20_1       mu_M2_20_2       log_sigma_M2_20_2       w_M2_20   ...
     mu_M2_10_1       log_sigma_M2_10_1       mu_M2_10_2       log_sigma_M2_10_2       w_M2_10   ...
     mu_M2_0_1        log_sigma_M2_0_1        mu_M2_0_2        log_sigma_M2_0_2        w_M2_0 ...
     gamma_cen        log_sigma_cen           elambda_cen      exi_cen;


parameters.sym  = [mu_wt;            log_sigma_wt; ... 
                   mu_M2_200_1;      log_sigma_M2_200_1;      mu_M2_200_2;      log_sigma_M2_200_2;      w_M2_200;   ...
                   mu_M2_80_1;       log_sigma_M2_80_1;       mu_M2_80_2;       log_sigma_M2_80_2;       w_M2_80;   ...
                   mu_M2_65_P50_1;   log_sigma_M2_65_P50_1;   mu_M2_65_P50_2;   log_sigma_M2_65_P50_2;   w_M2_65_P50;   ...
                   mu_M2_65_P188_1;  log_sigma_M2_65_P188_1;  mu_M2_65_P188_2;  log_sigma_M2_65_P188_2;  w_M2_65_P188;   ...
                   mu_M2_40_1;       log_sigma_M2_40_1;       mu_M2_40_2;       log_sigma_M2_40_2;       w_M2_40;   ...
                   mu_M2_20_1;       log_sigma_M2_20_1;       mu_M2_20_2;       log_sigma_M2_20_2;       w_M2_20;   ...
                   mu_M2_10_1;       log_sigma_M2_10_1;       mu_M2_10_2;       log_sigma_M2_10_2;       w_M2_10;   ...
                   mu_M2_0_1;        log_sigma_M2_0_1;        mu_M2_0_2;        log_sigma_M2_0_2;        w_M2_0; ...
                   gamma_cen;        log_sigma_cen;           elambda_cen;      exi_cen];


parameters.name = { '\mu_{wt}';            'log(\sigma_{wt})' ; ...
                    '\mu 1_{M2,200}';      'log(\sigma 1_{M2,200})' ;    '\mu_2_{M2,200}';      'log(\sigma 2_{M2,200})'   ;   'w_{M2,200}'   ; ...
                    '\mu 1_{M2,80}'   ;    'log(\sigma 1_{M2,80})'   ;   '\mu 2_{M2,80}'   ;    'log(\sigma 2_{M2,80})'   ;    'w_{M2,80}'   ; ...
                    '\mu 1_{M2,65_P50}' ;  'log(\sigma 1_{M2,65_P50})';  '\mu 2_{M2,65_P50}';   'log(\sigma 2_{M2,65_P50})' ;  'w_{M2,65_P50}'   ; ...
                    '\mu 1_{M2,65_P188}';  'log(\sigma 1_{M2,65_P188})'; '\mu 2_{M2,65_P188}';  'log(\sigma 2_{M2,65_P188})';  'w_{M2,65_P188}'   ; ...
                    '\mu 1_{M2,40}'   ;    'log(\sigma 1_{M2,40})'   ;   '\mu 2_{M2,40}'   ;    'log(\sigma 2_{M2,40})'   ;    'w_{M2,40}'   ; ...
                    '\mu 1_{M2,20}'   ;    'log(\sigma 1_{M2,20})'   ;   '\mu 2_{M2,20}'   ;    'log(\sigma 2_{M2,20})'   ;    'w_{M2,20}'   ; ...
                    '\mu 1_{M2,10}'    ;   'log(\sigma 1_{M2,10})'    ;  '\mu 2_{M2,10}'    ;   'log(\sigma 2_{M2,10})'    ;   'w_{M2,10}'    ; ...
                    '\mu 1_{M2,0}'    ;    'log(\sigma 1_{M2,0})'    ;   '\mu 2_{M2,0}'    ;    'log(\sigma 2_{M2,0})'    ;    'w_{M2,0}'    ; ...
                    'gamma_{cen}'  ;       'log(\sigma 1_{cen})'     ;   'log(\lambda_{cen})' ; 'log(\xi_{cen})' };


parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);

parameters.min  = [  log(5); log(1e-1);
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                     log(5); log(1e-1); log(5); log(1e-1); 0; ...
                        -20; log(1e-4); log(5); log(5)];

parameters.max  = [  log(2e3); log(1e1); 
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                     log(2e3); log(1e1); log(2e3); log(1e1); 1; ...
                           20; log(1e4); log(2e3); log(2e3)];


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
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_0,w_M2_0};
M.experiment(i).mu    = {mu_M2_0_1,mu_M2_0_2};
M.experiment(i).sigma = {exp(log_sigma_M2_0_1),exp(log_sigma_M2_0_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
M.experiment(i).mu    = {mu_M2_40_1,mu_M2_40_2};
M.experiment(i).sigma = {exp(log_sigma_M2_40_1),exp(log_sigma_M2_40_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


% EXPERIMENT
i = i + 1;
data_65_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P50,w_M2_65_P50};
M.experiment(i).mu    = {mu_M2_65_P50_1,mu_M2_65_P50_2};
M.experiment(i).sigma = {exp(log_sigma_M2_65_P50_1),exp(log_sigma_M2_65_P50_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
M.experiment(i).mu    = {mu_M2_65_P188_1,mu_M2_65_P188_2};
M.experiment(i).sigma = {exp(log_sigma_M2_65_P188_1),exp(log_sigma_M2_65_P188_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
M.experiment(i).mu    = {mu_M2_80_1,mu_M2_80_2};
M.experiment(i).sigma = {exp(log_sigma_M2_80_1),exp(log_sigma_M2_80_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_200_Mad2_P259bp;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_200,w_M2_200};
M.experiment(i).mu    = {mu_M2_200_1,mu_M2_200_2};
M.experiment(i).sigma = {exp(log_sigma_M2_200_1),exp(log_sigma_M2_200_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};




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
M.experiment(i).mu    = {mu_M2_10_1,mu_M2_10_2};
M.experiment(i).sigma = {exp(log_sigma_M2_10_1),exp(log_sigma_M2_10_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
M.experiment(i).mu    = {mu_M2_20_1,mu_M2_20_2};
M.experiment(i).sigma = {exp(log_sigma_M2_20_1),exp(log_sigma_M2_20_2)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};




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
M.experiment(i).mu    = {mu_wt};
M.experiment(i).sigma = {exp(log_sigma_wt)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
Mc = getMixtureModel(Mc,parameters.sym);

