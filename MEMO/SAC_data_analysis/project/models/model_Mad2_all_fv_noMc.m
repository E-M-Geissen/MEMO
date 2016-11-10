% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments
%   no censoring model Mc specified
 
%% SPECIFY MODEL PARAMETERS 
syms mu_WT          log_sigma_WT                 ...
     mu_M2_200      log_sigma_M2_200      w_M2_200  ...
     mu_M2_80       log_sigma_M2_80       w_M2_80   ...
     mu_M2_60_P50   log_sigma_M2_60_P50   w_M2_60_P50   ...
     mu_M2_60_P188  log_sigma_M2_60_P188  w_M2_60_P188   ...
     mu_M2_40       log_sigma_M2_40       w_M2_40   ...
     mu_M2_20       log_sigma_M2_20       w_M2_20   ...
     mu_M2_10       log_sigma_M2_10       w_M2_10    ...
     mu_M2_0        log_sigma_M2_0        w_M2_0;
 
 
parameters.sym  = [  mu_WT;          log_sigma_WT  ;               ...
                     mu_M2_200 ;     log_sigma_M2_200 ;     w_M2_200;  ...
                     mu_M2_80 ;      log_sigma_M2_80  ;     w_M2_80 ;  ...
                     mu_M2_60_P50;   log_sigma_M2_60_P50;   w_M2_60_P50;   ...
                     mu_M2_60_P188 ; log_sigma_M2_60_P188;  w_M2_60_P188;   ...
                     mu_M2_40 ;      log_sigma_M2_40 ;      w_M2_40 ;  ...
                     mu_M2_20 ;      log_sigma_M2_20 ;      w_M2_20 ;  ...
                     mu_M2_10  ;     log_sigma_M2_10 ;      w_M2_10 ;   ...
                     mu_M2_0 ;       log_sigma_M2_0;        w_M2_0 ];
                
                
parameters.name = {'mu_WT'         ; 'log_sigma_WT'         ;                 ...
                   'mu_{M2,200}'   ; 'log_sigma_{M2,200}'   ; 'w_{M2,200}'   ; ...
                   'mu_{M2,80}'   ; 'log_sigma_{M2,80}'   ; 'w_{M2,80}'   ; ...
                   'mu_{M2,60_P50}'   ; 'log_sigma_{M2,60_P50}'   ; 'w_{M2,60_P50}'   ; ...
                   'mu_{M2,60_P188}'   ; 'log_sigma_{M2,60_P188}'   ; 'w_{M2,60_P188}'   ; ...
                   'mu_{M2,40}'   ; 'log_sigma_{M2,40}'   ; 'w_{M2,40}'   ; ...
                   'mu_{M2,20}'   ; 'log_sigma_{M2,20}'   ; 'w_{M2,20}'   ; ...
                   'mu_{M2,10}'    ; 'log_sigma_{M2,10}'    ; 'w_{M2,10}'    ; ...
                   'mu_{M2,0}'    ; 'log_sigma_{M2,0}'    ; 'w_{M2,0}'     };
               
               
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-2);   ...
                     log(5); log(1e-2); 0; ...     
                     log(5); log(1e-2); 0; ...     
                     log(5); log(1e-2); 0; ...    
                     log(5); log(1e-2); 0; ...       
                     log(5); log(1e-2); 0; ...       
                     log(5); log(1e-2); 0; ...
                     log(5); log(1e-2); 0; ...
                     log(5); log(1e-2); 0];  
                 
parameters.max  = [ log(2e3); log(1e1);    ...
                    log(2e3); log(1e1); 1; ...
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
data_65_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60_P50,w_M2_60_P50};
M.experiment(i).mu    = {mu_WT,mu_M2_60_P50};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_60_P50)};



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
M.experiment(i).w     = {1-w_M2_60_P188,w_M2_60_P188};
M.experiment(i).mu    = {mu_WT,mu_M2_60_P188};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_60_P188)};



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
M.experiment(i).mu    = {mu_WT,mu_M2_80};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_80)};



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
M.experiment(i).mu    = {mu_WT,mu_M2_200};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_200)};



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
M.experiment(i).mu    = {mu_WT,mu_M2_0};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_0)};

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
M.experiment(i).mu    = {mu_WT,mu_M2_10};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_10)};


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
M.experiment(i).mu    = {mu_WT,mu_M2_20};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_20)};



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
M.experiment(i).mu    = {mu_WT,mu_M2_40};
M.experiment(i).sigma = {exp(log_sigma_WT),exp(log_sigma_M2_40)};



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



% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
