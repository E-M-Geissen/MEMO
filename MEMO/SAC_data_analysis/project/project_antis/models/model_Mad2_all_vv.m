% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   
%  two components variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_c     esigma_c             ...
     mu_M2_200_c   esigma_M2_200_c mu_M2_200   esigma_M2_200   w_M2_200  ...
     mu_M2_80_c esigma_M2_80_c mu_M2_80 esigma_M2_80 w_M2_80 ...
     mu_M2_60_P50_c esigma_M2_60_P50_c mu_M2_60_P50 esigma_M2_60_P50 w_M2_60_P50 ...
     mu_M2_60_P188_c esigma_M2_60_P188_c mu_M2_60_P188 esigma_M2_60_P188 w_M2_60_P188 ...
     mu_M2_40_c esigma_M2_40_c mu_M2_40 esigma_M2_40 w_M2_40 ...
     mu_M2_20_c esigma_M2_20_c mu_M2_20 esigma_M2_20 w_M2_20 ...
     mu_M2_10_c esigma_M2_10_c mu_M2_10 esigma_M2_10 w_M2_10 ...
     gamma_cen  esigma_cen     elambda_cen     exi_cen;
 
parameters.sym  = [ mu_c     ; esigma_c     ; ...
                    mu_M2_200_c;   esigma_M2_200_c; mu_M2_200;   esigma_M2_200;   w_M2_200;  ...
                    mu_M2_80_c ; esigma_M2_80_c ;mu_M2_80 ; esigma_M2_80 ; w_M2_80 ; ...
                    mu_M2_60_P50_c; esigma_M2_60_P50_c; mu_M2_60_P50; esigma_M2_60_P50; w_M2_60_P50; ...
                    mu_M2_60_P188_c; esigma_M2_60_P188_c; mu_M2_60_P188; esigma_M2_60_P188; w_M2_60_P188; ...
                    mu_M2_40_c ; esigma_M2_40_c ;mu_M2_40 ; esigma_M2_40 ; w_M2_40 ; ...
                    mu_M2_20_c ; esigma_M2_20_c ;mu_M2_20 ; esigma_M2_20 ; w_M2_20 ; ...
                    mu_M2_10_c ; esigma_M2_10_c ;mu_M2_10 ; esigma_M2_10 ; w_M2_10 ; ...
                    gamma_cen  ; esigma_cen     ; elambda_cen ; exi_cen];

parameters.name = {'mu_c'       ; 'esigma_c'       ; ...
                   'mu_c_{M2,200}' ; 'esigma_c_{M2,200}' ;'mu_{M2,200}' ; 'esigma_{M2,200}' ;'w_{M2,200}'; ...
                   'mu_c_{M2,80}' ; 'esigma_c_{M2,80}' ;'mu_{M2,80}' ; 'esigma_{M2,80}' ;'w_{M2,80}'; ...
                   'mu_c_{M2,60_P50}' ; 'esigma_c_{M2,60_P50}' ;'mu_{M2,60_P50}' ; 'esigma_{M2,60_P50}' ;'w_{M2,60_P50}'; ...
                   'mu_c_{M2,60_P188}' ; 'esigma_c_{M2,60_P188}' ;'mu_{M2,60_P188}' ; 'esigma_{M2,60_P188}' ;'w_{M2,60_P188}'; ... 
                   'mu_c_{M2,40}' ; 'esigma_c_{M2,40}' ;'mu_{M2,40}' ; 'esigma_{M2,40}' ;'w_{M2,40}'; ...
                   'mu_c_{M2,20}' ; 'esigma_c_{M2,20}' ;'mu_{M2,20}' ; 'esigma_{M2,20}' ;'w_{M2,20}'; ...
                    'mu_c_{M2,10}' ; 'esigma_c_{M2,10}' ;'mu_{M2,10}' ; 'esigma_{M2,10}' ;'w_{M2,10}'; ...
                    'gamma_{cen}'  ; 'esigma_{cen}'     ; 'elambda_{cen}' ; 'exi_{cen}' };

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [ log(5); log(1e-1);   ...
                    log(5); log(1e-1);log(5); log(1e-1);0; ...       
                    log(5); log(1e-1);log(5); log(1e-1);0; ...    
                    log(5); log(1e-1);log(5); log(1e-1);0; ...     
                    log(5); log(1e-1);log(5); log(1e-1);0; ...       
                    log(5); log(1e-1);log(5); log(1e-1);0; ...       
                    log(5); log(1e-1);log(5); log(1e-1);0; ...      
                    log(5); log(1e-1);log(5); log(1e-1);0; ...   
                     -20   ; log(1e-4); log(5); log(5)];      
                % As lower bound for the and mean we used the
                % inter-observation time.
parameters.max  = [ log(2e3); log(1e1);   ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    log(2e3); log(1e1);log(2e3); log(1e1);1; ...
                    20; log(1e4); log(2e3); log(2e3)];
        % As lower bound for the median we used twice
                % the observation time.
parameters.init_fun = @(theta_0,theta_min,theta_max,n) ...
            bsxfun(@max,bsxfun(@min,...
                        bsxfun(@plus,theta_0,[...
                                    randn(2,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                                                        randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                                                        randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);0.2*randn(1,n);...
                                    randn(4,n);...
                               ]),theta_max),theta_min);

%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
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
M.experiment(i).mu    = {mu_M2_200_c,mu_M2_200};
M.experiment(i).sigma = {exp(esigma_M2_200_c),exp(esigma_M2_200)};

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
M.experiment(i).w     = {1-w_M2_80,w_M2_80};
M.experiment(i).mu    = {mu_M2_80_c,mu_M2_80};
M.experiment(i).sigma = {exp(esigma_M2_80_c),exp(esigma_M2_80)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};



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
M.experiment(i).mu    = {mu_M2_60_P188_c,mu_M2_60_P188};
M.experiment(i).sigma = {exp(esigma_M2_60_P188_c),exp(esigma_M2_60_P188)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_60_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60_P50,w_M2_60_P50};
M.experiment(i).mu    = {mu_M2_60_P50_c,mu_M2_60_P50};
M.experiment(i).sigma = {exp(esigma_M2_60_P50_c),exp(esigma_M2_60_P50)};

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
M.experiment(i).w     = {1-w_M2_40,w_M2_40};
M.experiment(i).mu    = {mu_M2_40_c,mu_M2_40};
M.experiment(i).sigma = {exp(esigma_M2_40_c),exp(esigma_M2_40)};

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
M.experiment(i).w     = {1-w_M2_20,w_M2_20};
M.experiment(i).mu    = {mu_M2_20_c,mu_M2_20};
M.experiment(i).sigma = {exp(esigma_M2_20_c),exp(esigma_M2_20)};

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
M.experiment(i).w     = {1-w_M2_10,w_M2_10};
M.experiment(i).mu    = {mu_M2_10_c,mu_M2_10};
M.experiment(i).sigma = {exp(esigma_M2_10_c),exp(esigma_M2_10)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).mu    = {mu_c};
M.experiment(i).sigma = {exp(esigma_c)};

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
