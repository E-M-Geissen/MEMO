% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms gamma_wt  esigma_wt    elambda_wt     exi_wt              ...
         gamma1_M3_120   esigma1_M3_120     elambda1_M3_120      exi1_M3_120  gamma2_M3_120   esigma2_M3_120     elambda2_M3_120      exi2_M3_120    w_M3_120  ...
     gamma1_M3_60   esigma1_M3_60     elambda1_M3_60      exi1_M3_60   gamma2_M3_60   esigma2_M3_60     elambda2_M3_60      exi2_M3_60      w_M3_60  ...
     gamma1_M3_30   esigma1_M3_30     elambda1_M3_30      exi1_M3_30   gamma2_M3_30   esigma2_M3_30     elambda2_M3_30      exi2_M3_30   w_M3_30  ...
     gamma1_M3_0   esigma1_M3_0     elambda1_M3_0      exi1_M3_0 gamma2_M3_0   esigma2_M3_0     elambda2_M3_0      exi2_M3_0    w_M3_0  ...
     gamma_cen  esigma_cen     elambda_cen     exi_cen;
 
 
parameters.sym  = [ gamma_wt;  esigma_wt;    elambda_wt ;    exi_wt    ;          ...
     gamma1_M3_120 ;  esigma1_M3_120  ;   elambda1_M3_120 ;     exi1_M3_120 ; gamma2_M3_120 ;  esigma2_M3_120 ;    elambda2_M3_120   ;   exi2_M3_120 ;   w_M3_120;  ...
     gamma1_M3_60 ;  esigma1_M3_60  ;   elambda1_M3_60  ;    exi1_M3_60 ;  gamma2_M3_60;   esigma2_M3_60 ;    elambda2_M3_60  ;    exi2_M3_60    ;  w_M3_60 ; ...
     gamma1_M3_30 ;  esigma1_M3_30  ;   elambda1_M3_30   ;   exi1_M3_30  ; gamma2_M3_30;   esigma2_M3_30  ;   elambda2_M3_30   ;   exi2_M3_30 ;  w_M3_30 ; ...
     gamma1_M3_0  ; esigma1_M3_0 ;    elambda1_M3_0  ;    exi1_M3_0 ;gamma2_M3_0 ;  esigma2_M3_0  ;   elambda2_M3_0  ;    exi2_M3_0 ;   w_M3_0 ; ...
     gamma_cen ; esigma_cen ;    elambda_cen    ; exi_cen];
 
 
parameters.name = {'gamma_wt' ; 'esigma_wt' ;   'elambda_wt' ;    'exi_wt'    ;                     ...
                  'gamma1_{M3,120}'; 'esigma1_{M3,120}'; 'elambda1_{M3,120}'; 'exi1_{M3,120}'; 'gamma2_{M3,120}'; 'esigma2_{M3,120}'; 'elambda2_{M3,120}'; 'exi2_{M3,120}'; 'w_{M3,120}'; ...
                   'gamma1_{M3,60}'; 'esigma1_{M3,60}'; 'elambda1_{M3,60}'; 'exi1_{M3,60}';'gamma2_{M3,60}'; 'esigma2_{M3,60}'; 'elambda2_{M3,60}'; 'exi2_{M3,60}'; 'w_{M3,60}'; ...
                   'gamma1_{M3,30}'; 'esigma1_{M3,30}'; 'elambd1a_{M3,30}'; 'exi1_{M3,30}'; 'gamma2_{M3,30}'; 'esigma2_{M3,30}'; 'elambda2_{M3,30}'; 'exi2_{M3,30}';'w_{M3,30}'; ...
                   'gamma1_{M3,0}'; 'esigma1_{M3,0}'; 'elambda1_{M3,0}'; 'exi1_{M3,0}';'gamma2_{M3,0}'; 'esigma2_{M3,0}'; 'elambda2_{M3,0}'; 'exi2_{M3,0}';  'w_{M3,0}'; ...
                   'gamma_{cen}'  ; 'esigma_{cen}'     ; 'elambda_{cen}' ; 'exi_{cen}' };
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  -20   ; log(1e-4); log(5); log(5);  ...
                      -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...       
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5);-20   ; log(1e-4); log(5); log(5); 0; ...
                     -20   ; log(1e-4); log(5); log(5); -20   ; log(1e-4); log(5); log(5);0; ...
                     -20   ; log(1e-4); log(5); log(5)];  
                 
parameters.max  = [ 20; log(1e4); log(2e3); log(2e3);    ...
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
data_120_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_120,w_M3_120};
M.experiment(i).gamma  = { gamma2_M3_120,gamma1_M3_120 };
M.experiment(i).sigma  = { exp(esigma2_M3_120), exp(esigma1_M3_120)};
M.experiment(i).lambda = {exp(elambda2_M3_120), exp(elambda1_M3_120)};
M.experiment(i).xi     = {exp(exi2_M3_120), exp(exi1_M3_120)};

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
M.experiment(i).w     = {1-w_M3_60,w_M3_60};
M.experiment(i).gamma  = { gamma2_M3_60,gamma1_M3_60 };
M.experiment(i).sigma  = { exp(esigma2_M3_60), exp(esigma1_M3_60)};
M.experiment(i).lambda = {exp(elambda2_M3_60), exp(elambda1_M3_60)};
M.experiment(i).xi     = {exp(exi2_M3_60), exp(exi1_M3_60)};

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
M.experiment(i).w     = {1-w_M3_30,w_M3_30};
M.experiment(i).gamma  = { gamma2_M3_30,gamma1_M3_30};
M.experiment(i).sigma  = { exp(esigma2_M3_30), exp(esigma1_M3_30)};
M.experiment(i).lambda = {exp(elambda2_M3_30), exp(elambda1_M3_30)};
M.experiment(i).xi     = {exp(exi2_M3_30), exp(exi1_M3_30)};

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

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_0,w_M3_0};
M.experiment(i).gamma  = { gamma2_M3_0,gamma1_M3_0 };
M.experiment(i).sigma  = { exp(esigma2_M3_0), exp(esigma1_M3_0)};
M.experiment(i).lambda = {exp(elambda2_M3_0), exp(elambda1_M3_0)};
M.experiment(i).xi     = {exp(exi2_M3_0), exp(exi1_M3_0)};


Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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