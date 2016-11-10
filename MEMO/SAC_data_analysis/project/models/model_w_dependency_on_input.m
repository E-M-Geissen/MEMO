% Data:
%   all Mad2 and Mad3 data sets with reduced Mad abundance, data sets of
%   120% Mad3 and 200% Mad2 modeled as wild type 
%   
% Model:
%   log-normal mixtures of 2 components
%   Johnson SU distribution model for censoring times
%   one component fixed across experiments = wild type
%   one component variable across experiments
%   fraction of wild type like component as hill fuction of input = relative abundance of protein  (different parameters for dependency on Mad2 and Mad3)


%% SPECIFY MODEL PARAMETERS 
syms mu_WT       esigma_WT                 ...
     mu_M2_80   esigma_M2_80 e_M2_80 ...
     mu_M2_65_P50   esigma_M2_65_P50   e_M2_65_P50   ... 
     mu_M2_65_P188   esigma_M2_65_P188   e_M2_65_P188  ...
     mu_M2_40   esigma_M2_40     e_M2_40  ...
     mu_M2_20   esigma_M2_20     e_M2_20  ...
     mu_M2_10    esigma_M2_10     e_M2_10  ...
     mu_M2_0    esigma_M2_0        ...
     mu_M3_60   esigma_M3_60   e_M3_60   ...
     mu_M3_30   esigma_M3_30   e_M3_30   ...
     mu_M3_0 esigma_M3_0  ...
     gamma_cen  esigma_cen     elambda_cen     exi_cen ...
     K2 n2 K3 n3;

   

parameters.sym  = [ mu_WT       ; esigma_WT       ;             ...
                    mu_M2_80   ; esigma_M2_80;  e_M2_80; ...
                    mu_M2_65_P50   ; esigma_M2_65_P50   ;  e_M2_65_P50  ; ... 
                    mu_M2_65_P188    ; esigma_M2_65_P188    ;  e_M2_65_P188  ; ...    
                    mu_M2_40 ;  esigma_M2_40 ;    e_M2_40 ; ...
                    mu_M2_20;   esigma_M2_20;     e_M2_20 ; ...
                    mu_M2_10 ;   esigma_M2_10;     e_M2_10;  ...
                    mu_M2_0    ; esigma_M2_0    ;  ...
                    mu_M3_60 ;  esigma_M3_60;   e_M3_60;   ...
                    mu_M3_30 ;  esigma_M3_30 ;  e_M3_30;   ...
                    mu_M3_0; esigma_M3_0 ; ...
                    gamma_cen  ; esigma_cen     ; elambda_cen ; exi_cen; ...
                    K2; n2; K3; n3;];
                
parameters.name = {'mu_WT'          ; 'esigma_WT'          ;  ...
                   'mu_{M2,80}'     ; 'esigma_{M2,80}'     ; 'e_{M2_80}'; ...   
                   'mu_{M2,65_P50}' ; 'esigma_{M2,65_P50}' ; 'e{M2_65_P50}'  ;  ...
                   'mu_{M2,65_P188}'; 'esigma_{M2,65_P188}'; 'e{M2_65_P188}'  ;  ...
                   'mu_{M2,40}'     ; 'esigma_{M2,40}'     ; 'e_{M2_40}'; ...
                   'mu_{M2,20}'     ; 'esigma_{M2,20}'     ; 'e_{M2_20}'; ...
                   'mu_{M2,10}'     ; 'esigma_{M2,10}'     ; 'e_{M2_10}'; ...
                   'mu_{M2,0}'      ; 'esigma_{M2,0}'      ;  ...
                   'mu_{M3_60}'     ; 'esigma_{M3_60}'     ; 'e_{M3_60}';   ...
                   'mu_{M3_30'      ; 'esigma_{M3_30}'     ; 'e_{M3_30}';   ...
                   'mu_{M3_0}'      ; 'esigma_{M3_0}'      ;  ...
                   'gamma_{cen}'    ; 'esigma_{cen}'       ; 'elambda_{cen}' ; 'exi_{cen}' ; ...
                     'K_2' ; 'n_2'; 'K_3' ; 'n_3'  };
                 
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-2);   ...     
                     log(5); log(1e-2); -0.1;  ...     
                     log(5); log(1e-2); -0.1;  ... 
                     log(5); log(1e-2); -0.1;  ... 
                     log(5); log(1e-2); -0.1;  ...       
                     log(5); log(1e-2); -0.1;  ...
                     log(5); log(1e-2); -0.1; ...
                     log(5); log(1e-2);   ...     
                     log(5); log(1e-2); -0.1;  ...     
                     log(5); log(1e-2); -0.1;  ... 
                     log(5); log(1e-2);   ... 
                     -20   ; log(1e-4); log(5); log(5);...
                     0; 0.1;0; 0.1];  
                 
parameters.max  = [ log(2e3); log(1e1);    ...
                    log(2e3); log(1e1); 0.1; ...
                    log(2e3); log(1e1); 0.1;  ...
                    log(2e3); log(1e1); 0.1;  ...
                    log(2e3); log(1e1); 0.1;  ...
                    log(2e3); log(1e1); 0.1;  ...
                    log(2e3); log(1e1); 0.1;  ...
                    log(2e3); log(1e1); ...
                    log(2e3); log(1e1); 0.1; ...
                    log(2e3); log(1e1); 0.1;  ...
                    log(2e3); log(1e1); ...
                     20; log(1e4); log(2e3); log(2e3);...
                     1;50;  1;50];
                      


%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'Johnson SU';
Mc.label.x = 'time [min]';
Mc.label.y = 'probability density';





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
M.experiment(i).w     = {(Mad2+e_M2_65_P50)^n2*(1+K2^n2)/((Mad2+e_M2_65_P50)^n2+K2^n2),1-(Mad2+e_M2_65_P50)^n2*(1+K2^n2)/((Mad2+e_M2_65_P50)^n2+K2^n2)};
M.experiment(i).mu    = {mu_WT,mu_M2_65_P50};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_65_P50)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_65_P50};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).w     =  {(Mad2+e_M2_65_P188)^n2*(1+K2^n2)/((Mad2+e_M2_65_P188)^n2+K2^n2),1-(Mad2+e_M2_65_P188)^n2*(1+K2^n2)/((Mad2+e_M2_65_P188)^n2+K2^n2)};
M.experiment(i).mu    = {mu_WT,mu_M2_65_P188};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_65_P188)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_65_P188};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).w     = {(Mad2+e_M2_80)^n2*(1+K2^n2)/((Mad2+e_M2_80)^n2+K2^n2),1-(Mad2+e_M2_80)^n2*(1+K2^n2)/((Mad2+e_M2_80)^n2+K2^n2)};
M.experiment(i).mu    = {mu_WT,mu_M2_80};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_80)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_80};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).mu    = {mu_WT};
M.experiment(i).sigma = {exp(esigma_WT)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


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
M.experiment(i).w     = {0,1};
M.experiment(i).mu    = {mu_WT,mu_M2_0};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_0)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).w     = {(Mad2+e_M2_10)^n2*(1+K2^n2)/((Mad2+e_M2_10)^n2+K2^n2),1-(Mad2+e_M2_10)^n2*(1+K2^n2)/((Mad2+e_M2_10)^n2+K2^n2)};
M.experiment(i).mu    = {mu_WT,mu_M2_10};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_10)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_10};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).w     = {(Mad2+e_M2_20)^n2*(1+K2^n2)/((Mad2+e_M2_20)^n2+K2^n2),1-(Mad2+e_M2_20)^n2*(1+K2^n2)/((Mad2+e_M2_20)^n2+K2^n2)};
M.experiment(i).mu    = {mu_WT,mu_M2_20};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_20)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_20};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).w     = {(Mad2+e_M2_40)^n2*(1+K2^n2)/((Mad2+e_M2_40)^n2+K2^n2),1-(Mad2+e_M2_40)^n2*(1+K2^n2)/((Mad2+e_M2_40)^n2+K2^n2)};
M.experiment(i).mu    = {mu_WT,mu_M2_40};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M2_40)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {e_M2_40};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K2','n2'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


% EXPERIMENT
i = i + 1;
data_Mad3_delta;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;
M.experiment(i).cond.y = {Mad3};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {};

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {0,1};
M.experiment(i).mu    = {mu_WT,mu_M3_0};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M3_0)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
M.experiment(i).w     = {(Mad3+e_M3_30)^n3*(1+K3^n3)/((Mad3+e_M3_30)^n3+K3^n3),1-(Mad3+e_M3_30)^n3*(1+K3^n3)/((Mad3+e_M3_30)^n3+K3^n3)};
M.experiment(i).mu    = {mu_WT,mu_M3_30};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M3_30)};
M.experiment(i).cond.y = {Mad3};
M.experiment(i).cond.e = {e_M3_30};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K3','n3'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
M.experiment(i).w     =  {(Mad3+e_M3_60)^n3*(1+K3^n3)/((Mad3+e_M3_60)^n3+K3^n3),1-(Mad3+e_M3_60)^n3*(1+K3^n3)/((Mad3+e_M3_60)^n3+K3^n3)};
M.experiment(i).mu    = {mu_WT,mu_M3_60};
M.experiment(i).sigma = {exp(esigma_WT),exp(esigma_M3_60)};
M.experiment(i).cond.y = {Mad3};
M.experiment(i).cond.e = {e_M3_60};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K3','n3'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_120_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).mu    = {mu_WT};
M.experiment(i).sigma = {exp(esigma_WT)};
M.experiment(i).cond.y = {Mad3};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {'K3','n3'};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};




% EXPERIMENT
i = i+1;
data_WT_fus;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).mu    = {mu_WT};
M.experiment(i).sigma = {exp(esigma_WT)};
M.experiment(i).cond.y = {Mad2};
M.experiment(i).cond.e = {0};
M.experiment(i).cond.sigma = {0.02};
M.experiment(i).cond.param = {};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};


% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
Mc = getMixtureModel(Mc,parameters.sym);