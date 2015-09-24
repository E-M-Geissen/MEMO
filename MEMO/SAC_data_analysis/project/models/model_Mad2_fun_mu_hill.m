% Data:
%   control
%   Mad2 data no 200% data and no double mutants
% Model:
%   log-normal
%   one component fixed across experiments = wild type
%   one component variable across experiments = strain specific
%   parametrization of mu of log-normal distribution via hill function mu=mu_null+vmax*(Mad2+error)^nh/(KM2^nh+(Mad2+error)^nh)
 
%% SPECIFY MODEL PARAMETERS 
syms mu_WT       log_sigma_WT                 ...
                 log_sigma_M2_80       w_M2_80      e_M2_80   ...
                 log_sigma_M2_65_P50   w_M2_65_P50  e_M2_65_P50   ...
                 log_sigma_M2_65_P188  w_M2_65_P188 e_M2_65_P188  ...
                 log_sigma_M2_40       w_M2_40      e_M2_40 ...
                 log_sigma_M2_20       w_M2_20      e_M2_20 ...
                 log_sigma_M2_10       w_M2_10      e_M2_10  ...
                 log_sigma_M2_0        w_M2_0    ...
     mu_null KM2 nh vmax...
     gamma_cen  log_sigma_cen     log_lambda_cen     log_xi_cen;
 
 
parameters.sym  = [  mu_WT       ; log_sigma_WT       ;             ...
                                   log_sigma_M2_80       ; w_M2_80       ; e_M2_80;...
                                   log_sigma_M2_65_P50   ; w_M2_65_P50   ; e_M2_65_P50  ; ...
                                   log_sigma_M2_65_P188  ; w_M2_65_P188  ; e_M2_65_P188  ; ...
                                   log_sigma_M2_40       ; w_M2_40       ; e_M2_40 ;...
                                   log_sigma_M2_20       ; w_M2_20       ; e_M2_20 ; ...
                                   log_sigma_M2_10       ; w_M2_10       ; e_M2_10;...
                                    log_sigma_M2_0        ; w_M2_0    ; ...
                     mu_null; KM2 ; nh; vmax; ...
                     gamma_cen  ; log_sigma_cen ; log_lambda_cen ; log_xi_cen];
                 
parameters.name = {'mu_WT'         ; 'log_sigma_WT'         ;                 ...
                                     'log(\sigma_{M2,80})'      ; 'w_{M2,80}'       ; 'e_{M2_80}';...
                                     'log(\sigma_{M2,65_P50})'  ; 'w_{M2,65_P50}'   ; 'e{M2_65_P50}'  ; ...
                                     'log(\sigma_{M2,65_P188})' ; 'w_{M2,65_P188}'  ; 'e{M2_65_P188}'  ; ...
                                     'log(\sigma_{M2,40})'      ; 'w_{M2,40}'       ; 'e_{M2_40}'; ...
                                     'log(\sigma_{M2,20})'      ; 'w_{M2,20}'       ; 'e_{M2_20}';...
                                     'log(\sigma_{M2,10})'      ; 'w_{M2,10}'       ; 'e_{M2_10}'; ...
                                     'log(\sigma_{M2,0})'       ; 'w_{M2,0}'        ; ...
                    'mu_null'; 'KM_2'; 'nh'; 'v_{max}';... 
                    '\gamma_{cen}'  ; 'log(\sigma_{cen})'     ; 'log(\lambda_{cen}9' ; 'log(\xi_{cen})' };
                
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [   log(5);log(1e-2);   ...     
                      log(1e-2); 0; -0.1;...     
                      log(1e-2); 0;-0.1; ...    
                      log(1e-2); 0; -0.1;...       
                      log(1e-2); 0; -0.1;...       
                      log(1e-2); 0; -0.1;...
                      log(1e-2); 0; -0.1;...
                      log(1e-2); 0; ...  
                      1; 0 ; 1;0.1; ...
                     -20   ; log(1e-4); log(5); log(5)];  
                 
parameters.max  = [ log(2e3); log(1e1);    ...
                    log(1e1); 1; 0.1;...
                     log(1e1); 1; 0.1;...
                     log(1e1); 1; 0.1;...
                     log(1e1); 1;0.1; ...
                     log(1e1); 1; 0.1;...
                     log(1e1); 1; 0.1;...
                     log(1e1); 1; ...
                     5;5;5;5;
                     20; log(1e4); log(2e3); log(2e3)];
                      


%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'Johnson SU';
Mc.label.x = 'time [min]';
Mc.label.y = 'probability density';



% EXPERIMENT
i =  1;
data_Mad2_delta;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_0,w_M2_0};
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2)/(KM2+Mad2)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_0)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};



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
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2+e_M2_10)^nh/(KM2^nh+(Mad2+e_M2_10)^nh)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_10)};

M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_10};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};

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
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2+e_M2_20)^nh/(KM2^nh+(Mad2+e_M2_20)^nh)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_20)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_20};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};

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
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2+e_M2_40)^nh/(KM2^nh+(Mad2+e_M2_40)^nh)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_40)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_40};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(log_sigma_cen)};
Mc.experiment(i).lambda = {exp(log_lambda_cen)};
Mc.experiment(i).xi     = {    exp(log_xi_cen)};

% EXPERIMENT
i =i+ 1 ;
data_65_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_65_P50,w_M2_65_P50};
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2+e_M2_65_P50)^nh/(KM2^nh+(Mad2+e_M2_65_P50)^nh)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_65_P50)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_65_P50};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};

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
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2+e_M2_65_P188)^nh/(KM2^nh+(Mad2+e_M2_65_P188)^nh)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_65_P188)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_65_P188};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};

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
M.experiment(i).mu    = {mu_WT,mu_null+vmax*(Mad2+e_M2_80)^nh/(KM2^nh+(Mad2+e_M2_80)^nh)};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_80)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_80};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'mu_null', 'KM2', 'nh', 'vmax'};

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
M.experiment(i).mu    = {mu_WT};
M.experiment(i).sigma = {exp(log_sigma_WT)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {};

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