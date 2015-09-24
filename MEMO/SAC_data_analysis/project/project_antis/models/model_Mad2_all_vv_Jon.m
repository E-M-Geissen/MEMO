% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms gamma_wt  esigma_wt    elambda_wt     exi_wt              ...
     gamma1_M2_200   esigma1_M2_200     elambda1_M2_200      exi1_M2_200  gamma2_M2_200   esigma2_M2_200     elambda2_M2_200      exi2_M2_200    w_M2_200  ...
     gamma1_M2_80   esigma1_M2_80     elambda1_M2_80      exi1_M2_80 gamma2_M2_80   esigma2_M2_80     elambda2_M2_80      exi2_M2_80    w_M2_80  ...
     gamma1_M2_60_P50   esigma1_M2_60_P50     elambda1_M2_60_P50      exi1_M2_60_P50  gamma2_M2_60_P50   esigma2_M2_60_P50     elambda2_M2_60_P50      exi2_M2_60_P50    w_M2_60_P50  ...
     gamma1_M2_60_P188   esigma1_M2_60_P188    elambda1_M2_60_P188      exi1_M2_60_P188      gamma2_M2_60_P188   esigma2_M2_60_P188    elambda2_M2_60_P188      exi2_M2_60_P188   w_M2_60_P188  ...
     gamma1_M2_40   esigma1_M2_40     elambda1_M2_40      exi1_M2_40  gamma2_M2_40   esigma2_M2_40     elambda2_M2_40      exi2_M2_40    w_M2_40  ...
     gamma1_M2_20   esigma1_M2_20     elambda1_M2_20      exi1_M2_20   gamma2_M2_20   esigma2_M2_20     elambda2_M2_20      exi2_M2_20      w_M2_20  ...
     gamma1_M2_10   esigma1_M2_10     elambda1_M2_10      exi1_M2_10   gamma2_M2_10   esigma2_M2_10     elambda2_M2_10      exi2_M2_10   w_M2_10  ...
     gamma1_M2_0   esigma1_M2_0     elambda1_M2_0      exi1_M2_0 gamma2_M2_0   esigma2_M2_0     elambda2_M2_0      exi2_M2_0    w_M2_0  ...
     gamma_cen  esigma_cen     elambda_cen     exi_cen;
 
 
parameters.sym  = [ gamma_wt;  esigma_wt;    elambda_wt ;    exi_wt    ;          ...
     gamma1_M2_200 ;  esigma1_M2_200  ;   elambda1_M2_200  ;    exi1_M2_200;  gamma2_M2_200 ;  esigma2_M2_200  ;   elambda2_M2_200   ;   exi2_M2_200   ; w_M2_200 ; ...
     gamma1_M2_80  ; esigma1_M2_80   ;  elambda1_M2_80 ;     exi1_M2_80; gamma2_M2_80  ; esigma2_M2_80;     elambda2_M2_80 ;     exi2_M2_80  ;  w_M2_80 ; ...
     gamma1_M2_60_P50 ;  esigma1_M2_60_P50 ;    elambda1_M2_60_P50   ;   exi1_M2_60_P50;  gamma2_M2_60_P50  ; esigma2_M2_60_P50   ;  elambda2_M2_60_P50;      exi2_M2_60_P50 ;   w_M2_60_P50;  ...
     gamma1_M2_60_P188 ;  esigma1_M2_60_P188 ;   elambda1_M2_60_P188  ;    exi1_M2_60_P188   ;   gamma2_M2_60_P188 ;  esigma2_M2_60_P188 ;   elambda2_M2_60_P188   ;   exi2_M2_60_P188 ;  w_M2_60_P188 ; ...
     gamma1_M2_40 ;  esigma1_M2_40  ;   elambda1_M2_40 ;     exi1_M2_40 ; gamma2_M2_40 ;  esigma2_M2_40 ;    elambda2_M2_40   ;   exi2_M2_40 ;   w_M2_40;  ...
     gamma1_M2_20 ;  esigma1_M2_20  ;   elambda1_M2_20  ;    exi1_M2_20 ;  gamma2_M2_20;   esigma2_M2_20 ;    elambda2_M2_20  ;    exi2_M2_20    ;  w_M2_20 ; ...
     gamma1_M2_10 ;  esigma1_M2_10  ;   elambda1_M2_10   ;   exi1_M2_10  ; gamma2_M2_10;   esigma2_M2_10  ;   elambda2_M2_10   ;   exi2_M2_10 ;  w_M2_10 ; ...
     gamma1_M2_0  ; esigma1_M2_0 ;    elambda1_M2_0  ;    exi1_M2_0 ;gamma2_M2_0 ;  esigma2_M2_0  ;   elambda2_M2_0  ;    exi2_M2_0 ;   w_M2_0 ; ...
     gamma_cen ; esigma_cen ;    elambda_cen    ; exi_cen];
 
 
parameters.name = {'gamma_wt' ; 'esigma_wt' ;   'elambda_wt' ;    'exi_wt'    ;                     ...
                   'gamma1_{M2,200}'; 'esigma1_{M2,200}'; 'elambda1_{M2,200}'; 'exi1_{M2,200}'; 'gamma1_{M2,200}'; 'esigma1_{M2,200}'; 'elambda1_{M2,200}'; 'exi1_{M2,200}'; 'w_{M2,200}'; ...
                   'gamma1_{M2,80}'; 'esigma1_{M2,80}'; 'elambda1_{M2,80}'; 'exi1_{M2,80}'; 'gamma2_{M2,80}'; 'esigma2_{M2,80}'; 'elambda2_{M2,80}'; 'exi2_{M2,80}';'w_{M2,80}'; ...
                   'gamma1_{M2,60_P50}'; 'esigma1_{M2,60_P50}'; 'elambda1_{M2,60_P50}'; 'exi1_{M2,60_P50}'; 'gamma2_{M2,60_P50}'; 'esigma2_{M2,60_P50}'; 'elambda2_{M2,60_P50}'; 'exi2_{M2,60_P50}'; 'w_{M2,60_P50}'; ...
                   'gamma1_{M2,60_P188}'; 'esigma1_{M2,60_P188}'; 'elambda1_{M2,60_P188}'; 'exi1_{M2,60_P188}'; 'gamma2_{M2,60_P188}'; 'esigma2_{M2,60_P188}'; 'elambda2_{M2,60_P188}'; 'exi2_{M2,60_P188}'; 'w_{M2,60_P188}'; ...
                   'gamma1_{M2,40}'; 'esigma1_{M2,40}'; 'elambda1_{M2,40}'; 'exi1_{M2,40}'; 'gamma2_{M2,40}'; 'esigma2_{M2,40}'; 'elambda2_{M2,40}'; 'exi2_{M2,40}'; 'w_{M2,40}'; ...
                   'gamma1_{M2,20}'; 'esigma1_{M2,20}'; 'elambda1_{M2,20}'; 'exi1_{M2,20}';'gamma2_{M2,20}'; 'esigma2_{M2,20}'; 'elambda2_{M2,20}'; 'exi2_{M2,20}'; 'w_{M2,20}'; ...
                   'gamma1_{M2,10}'; 'esigma1_{M2,10}'; 'elambd1a_{M2,10}'; 'exi1_{M2,10}'; 'gamma2_{M2,10}'; 'esigma2_{M2,10}'; 'elambda2_{M2,10}'; 'exi2_{M2,10}';'w_{M2,10}'; ...
                   'gamma1_{M2,0}'; 'esigma1_{M2,0}'; 'elambda1_{M2,0}'; 'exi1_{M2,0}';'gamma2_{M2,0}'; 'esigma2_{M2,0}'; 'elambda2_{M2,0}'; 'exi2_{M2,0}';  'w_{M2,0}'; ...
                   'gamma_{cen}'  ; 'esigma_{cen}'     ; 'elambda_{cen}' ; 'exi_{cen}' };
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  -20   ; log(1e-4); log(5); log(5);  ...
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...     
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...     
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...    
                     -20   ; log(1e-4); log(5); log(5); -20   ; log(1e-4); log(5); log(5);0; ...       
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...       
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5); -20   ; log(1e-4); log(5); log(5);0; ...
                     -20   ; log(1e-4); log(5); log(5)];  
                 
parameters.max  = [ 20; log(1e4); log(2e3); log(2e3);    ...
                    20; log(1e4); log(2e3); log(2e3);20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3);20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 20; log(1e4); log(2e3); log(2e3);1; ...
                    20; log(1e4); log(2e3); log(2e3);20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3);20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 20; log(1e4); log(2e3); log(2e3);1; ...
                    20; log(1e4); log(2e3); log(2e3);20; log(1e4); log(2e3); log(2e3); 1; ...
                    20; log(1e4); log(2e3); log(2e3); 20; log(1e4); log(2e3); log(2e3);1; ...
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
M.experiment(i).gamma  = { gamma2_M2_200,gamma1_M2_200 };
M.experiment(i).sigma  = { exp(esigma2_M2_200), exp(esigma1_M2_200)};
M.experiment(i).lambda = {exp(elambda2_M2_200), exp(elambda1_M2_200)};
M.experiment(i).xi     = {exp(exi2_M2_200), exp(exi1_M2_200)};


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
M.experiment(i).gamma  = { gamma2_M2_80,gamma1_M2_80};
M.experiment(i).sigma  = { exp(esigma2_M2_80), exp(esigma1_M2_80)};
M.experiment(i).lambda = {exp(elambda2_M2_80), exp(elambda1_M2_80)};
M.experiment(i).xi     = {exp(exi2_M2_80), exp(exi1_M2_80)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i+1 ;
data_60_P50_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60_P50,w_M2_60_P50};
M.experiment(i).gamma  = { gamma2_M2_60_P50,gamma1_M2_60_P50};
M.experiment(i).sigma  = { exp(esigma2_M2_60_P50), exp(esigma1_M2_60_P50)};
M.experiment(i).lambda = {exp(elambda2_M2_60_P50), exp(elambda1_M2_60_P50)};
M.experiment(i).xi     = {exp(exi2_M2_60_P50), exp(exi1_M2_60_P50)};


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
M.experiment(i).gamma  = { gamma2_M2_200,gamma1_M2_200 };
M.experiment(i).sigma  = { exp(esigma2_M2_60_P188), exp(esigma1_M2_60_P188)};
M.experiment(i).lambda = {exp(elambda2_M2_60_P188), exp(elambda1_M2_60_P188)};
M.experiment(i).xi     = {exp(exi2_M2_60_P188), exp(exi1_M2_60_P188)};

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
M.experiment(i).gamma  = { gamma2_M2_40,gamma1_M2_40 };
M.experiment(i).sigma  = { exp(esigma2_M2_40), exp(esigma1_M2_40)};
M.experiment(i).lambda = {exp(elambda2_M2_40), exp(elambda1_M2_40)};
M.experiment(i).xi     = {exp(exi2_M2_40), exp(exi1_M2_40)};

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
M.experiment(i).gamma  = { gamma2_M2_20,gamma1_M2_20 };
M.experiment(i).sigma  = { exp(esigma2_M2_20), exp(esigma1_M2_20)};
M.experiment(i).lambda = {exp(elambda2_M2_20), exp(elambda1_M2_20)};
M.experiment(i).xi     = {exp(exi2_M2_20), exp(exi1_M2_20)};

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
M.experiment(i).gamma  = { gamma2_M2_10,gamma1_M2_10};
M.experiment(i).sigma  = { exp(esigma2_M2_10), exp(esigma1_M2_10)};
M.experiment(i).lambda = {exp(elambda2_M2_10), exp(elambda1_M2_10)};
M.experiment(i).xi     = {exp(exi2_M2_10), exp(exi1_M2_10)};

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
M.experiment(i).w     = {1-w_M2_0,w_M2_0};
M.experiment(i).gamma  = { gamma2_M2_0,gamma1_M2_0 };
M.experiment(i).sigma  = { exp(esigma2_M2_0), exp(esigma1_M2_0)};
M.experiment(i).lambda = {exp(elambda2_M2_0), exp(elambda1_M2_0)};
M.experiment(i).xi     = {exp(exi2_M2_0), exp(exi1_M2_0)};


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
M.experiment(i).gamma  = {      gamma_wt };
M.experiment(i).sigma  = { exp(esigma_wt)};
M.experiment(i).lambda = {exp(elambda_wt)};
M.experiment(i).xi     = {    exp(exi_wt)};

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