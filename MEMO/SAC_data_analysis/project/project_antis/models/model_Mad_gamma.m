%% SPECIFY MODEL PARAMETERS 
syms log_alpha_c     log_beta_c                 ...
     log_alpha_M2_60 log_beta_M2_60   w_M2_60   ...
     log_threshold;
parameters.sym  = [ log_alpha_c       ; log_beta_c       ;             ...
                    log_alpha_M2_60   ; log_beta_M2_60   ; w_M2_60   ; ...
                    log_threshold   ];
parameters.name = {'log(\alpha_c)'         ; 'log(\beta_c)'         ;                 ...
                   'log(\alpha_{M2,60})'   ; 'log(\beta_{M2,60})'   ; 'w_{M2,60}'   ; ...
                   'log(threshold)'    };
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-1);   ...
                     log(5); log(1e-1); 0; ...       
                     log(5)];       
parameters.max  = [ log(2e3); log(1e1);    ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3)];
parameters.guess = parameters.max - parameters.min;

%% SPECIFY MODEL AND DATA
M.mixture.type = 'gamma';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'delta';
Mc.label.x = 'time [min]';
Mc.label.y = 'probability density';

% EXPERIMENT
i = 1;
data_WT_fus;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).alpha = {exp(log_alpha_c)};
M.experiment(i).beta  = {exp(log_beta_c)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).location = {exp(log_threshold)};

% EXPERIMENT
i = i + 1;
data_60_Mad2_new;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 0;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60,w_M2_60};
M.experiment(i).alpha = {exp(log_alpha_c),exp(log_alpha_M2_60)};
M.experiment(i).beta  = {exp(log_beta_c) ,exp(log_beta_M2_60) };

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).location = {exp(log_threshold)};

% Compile model
% (This generates the functional expression of parameters and derivatives.)
M  = getMixtureModel(M ,parameters.sym);
Mc = getMixtureModel(Mc,parameters.sym);
