% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_c       esigma_c                 ...
     mu_M2_200   esigma_M2_200   w_M2_200  ...
     mu_M2_80   esigma_M2_80   w_M2_80   ...
     mu_M2_60_P50   esigma_M2_60_P50   w_M2_60_P50   ...
     mu_M2_60_P188   esigma_M2_60_P188    w_M2_60_P188   ...
     mu_M2_40   esigma_M2_40   w_M2_40   ...
     mu_M2_20   esigma_M2_20   w_M2_20   ...
     mu_M2_10    esigma_M2_10    w_M2_10;
 
 
parameters.sym  = [ mu_c       ; esigma_c       ;             ...
                    mu_M2_200   ; esigma_M2_200   ; w_M2_200   ; ...
                    mu_M2_80   ; esigma_M2_80   ; w_M2_80   ; ...
                    mu_M2_60_P50   ; esigma_M2_60_P50   ; w_M2_60_P50   ; ...
                    mu_M2_60_P188    ; esigma_M2_60_P188    ; w_M2_60_P188   ; ...
                    mu_M2_40   ; esigma_M2_40   ; w_M2_40   ; ...
                    mu_M2_20   ; esigma_M2_20   ; w_M2_20   ; ...
                    mu_M2_10    ; esigma_M2_10    ; w_M2_10    ];
parameters.name = {'mu_c'         ; 'esigma_c'         ;                 ...
                   'mu_{M2,200}'   ; 'esigma_{M2,200}'   ; 'w_{M2,200}'   ; ...
                   'mu_{M2,80}'   ; 'esigma_{M2,80}'   ; 'w_{M2,80}'   ; ...
                   'mu_{M2,60_P50}'   ; 'esigma_{M2,60_P50}'   ; 'w_{M2,60_P50}'   ; ...
                   'mu_{M2,60_P188}'   ; 'esigma_{M2,60_P188}'   ; 'w_{M2,60_P188}'   ; ...
                   'mu_{M2,40}'   ; 'esigma_{M2,40}'   ; 'w_{M2,40}'   ; ...
                   'mu_{M2,20}'   ; 'esigma_{M2,20}'   ; 'w_{M2,20}'   ; ...
                   'mu_{M2,10}'    ; 'esigma_{M2,10}'    ; 'w_{M2,10}'     };
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-5);   ...
                     log(5); log(1e-5); 0; ...     
                     log(5); log(1e-5); 0; ...     
                     log(5); log(1e-5); 0; ...    
                     log(5); log(1e-5); 0; ...       
                     log(5); log(1e-5); 0; ...       
                     log(5); log(1e-5); 0; ...
                     log(5); log(1e-5); 0];  
                 
parameters.max  = [ log(2e3); log(1e1);    ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1];
                      


%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';







% EXPERIMENT
i = 1 ;
data_60_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60_P50,w_M2_60_P50};
M.experiment(i).mu    = {mu_c,mu_M2_60_P50};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_60_P50)};



% EXPERIMENT
i = i + 1;
data_60_P188_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60_P188,w_M2_60_P188};
M.experiment(i).mu    = {mu_c,mu_M2_60_P188};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_60_P188)};


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
M.experiment(i).mu    = {mu_c,mu_M2_80};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_80)};



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
M.experiment(i).mu    = {mu_c,mu_M2_200};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_200)};














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
M.experiment(i).mu    = {mu_c,mu_M2_10};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_10)};



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
M.experiment(i).mu    = {mu_c,mu_M2_20};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_20)};




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
M.experiment(i).mu    = {mu_c,mu_M2_40};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_40)};



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
M.experiment(i).mu    = {mu_c};
M.experiment(i).sigma = {exp(esigma_c)};




% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
