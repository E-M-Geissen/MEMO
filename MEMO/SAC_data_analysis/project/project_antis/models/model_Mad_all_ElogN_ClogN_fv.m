% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_c       esigma_c                 ...
     mu_M1_30   esigma_M1_30   w_M1_30   ...
     mu_M1_10   esigma_M1_10   w_M1_10   ...
     mu_M2_80   esigma_M2_80   w_M2_80   ...
     mu_M2_60   esigma_M2_60   w_M2_60   ...
     mu_M2_25   esigma_M2_25   w_M2_25   ...
     mu_M2_0    esigma_M2_0    w_M2_0    ...
     mu_M3_60   esigma_M3_60   w_M3_60   ...
     mu_M3_30_1 esigma_M3_30_1 w_M3_30_1 ...
     mu_M3_30_2 esigma_M3_30_2 w_M3_30_2 ...
     mu_M3_0    esigma_M3_0    w_M3_0    ...
     mu_cen     esigma_cen;
parameters.sym  = [ mu_c       ; esigma_c       ;             ...
                    mu_M1_30   ; esigma_M1_30   ; w_M1_30   ; ...
                    mu_M1_10   ; esigma_M1_10   ; w_M1_10   ; ...
                    mu_M2_80   ; esigma_M2_80   ; w_M2_80   ; ...
                    mu_M2_60   ; esigma_M2_60   ; w_M2_60   ; ...
                    mu_M2_25   ; esigma_M2_25   ; w_M2_25   ; ...
                    mu_M2_0    ; esigma_M2_0    ; w_M2_0    ; ...
                    mu_M3_60   ; esigma_M3_60   ; w_M3_60   ; ...
                    mu_M3_30_1 ; esigma_M3_30_1 ; w_M3_30_1 ; ...
                    mu_M3_30_2 ; esigma_M3_30_2 ; w_M3_30_2 ; ...
                    mu_M3_0    ; esigma_M3_0    ; w_M3_0    ; ...
                    mu_cen     ; esigma_cen];
parameters.name = {'mu_c'         ; 'esigma_c'         ;                 ...
                   'mu_{M1,30}'   ; 'esigma_{M1,30}'   ; 'w_{M1,30}'   ; ...
                   'mu_{M1,10}'   ; 'esigma_{M1,10}'   ; 'w_{M1,10}'   ; ...
                   'mu_{M2,80}'   ; 'esigma_{M2,80}'   ; 'w_{M2,80}'   ; ...
                   'mu_{M2,60}'   ; 'esigma_{M2,60}'   ; 'w_{M2,60}'   ; ...
                   'mu_{M2,25}'   ; 'esigma_{M2,25}'   ; 'w_{M2,25}'   ; ...
                   'mu_{M2,0}'    ; 'esigma_{M2,0}'    ; 'w_{M2,0}'    ; ...
                   'mu_{M3,60}'   ; 'esigma_{M3,60}'   ; 'w_{M3,60}'   ; ...
                   'mu_{M3,30,1}' ; 'esigma_{M3,30,1}' ; 'w_{M3,30,1}' ; ...
                   'mu_{M3,30,2}' ; 'esigma_{M3,30,2}' ; 'w_{M3,30,2}' ; ...
                   'mu_{M3,0}'    ; 'esigma_{M3,0}'    ; 'w_{M3,0}'    ; ...
                   'mu_{cen}'     ; 'esigma_{cen}'    };
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-1);   ...
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1); 0; ...       
                     log(5); log(1e-1)];       
parameters.max  = [ log(2e3); log(1e1);    ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1); 1; ...
                    log(2e3); log(1e1)];
parameters.init_fun = @(theta_0,theta_min,theta_max,n) ...
            bsxfun(@max,bsxfun(@min,...
                        bsxfun(@plus,theta_0,[...
                                    randn(2,n);...
                                    randn(2,n);0.2*randn(1,n);...
                                    randn(2,n);0.2*randn(1,n);...
                                    randn(2,n);0.2*randn(1,n);...
                                    randn(2,n);0.2*randn(1,n);...
                                    randn(2,1);0.2*randn(1,1);...
                                    randn(2,1);0.2*randn(1,1);...
                                    randn(2,1);0.2*randn(1,1);...
                                    randn(2,1);0.2*randn(1,1);...
                                    randn(2,1);0.2*randn(1,1);...
                                    randn(2,1);0.2*randn(1,1);...
                                    randn(2,n);...
                               ]),theta_max),theta_min);

%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'log-normal';
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
M.experiment(i).mu    = {mu_c};
M.experiment(i).sigma = {exp(esigma_c)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad1;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M1_30,w_M1_30};
M.experiment(i).mu    = {mu_c,mu_M1_30};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M1_30)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% EXPERIMENT
i = i + 1;
data_10_Mad1;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M1_10,w_M1_10};
M.experiment(i).mu    = {mu_c,mu_M1_10};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M1_10)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

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

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% EXPERIMENT
i = i + 1;
data_60_Mad2_new;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_60,w_M2_60};
M.experiment(i).mu    = {mu_c,mu_M2_60};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_60)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% EXPERIMENT
i = i + 1;
data_25_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M2_25,w_M2_25};
M.experiment(i).mu    = {mu_c,mu_M2_25};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_25)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

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
M.experiment(i).mu    = {mu_c,mu_M2_0};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_0)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

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
M.experiment(i).mu    = {mu_c,mu_M3_60};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M3_60)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad3_1;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_30_1,w_M3_30_1};
M.experiment(i).mu    = {mu_c,mu_M3_30_1};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M3_30_1)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad3_2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w_M3_30_2,w_M3_30_2};
M.experiment(i).mu    = {mu_c,mu_M3_30_2};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M3_30_2)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

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
M.experiment(i).mu    = {mu_c,mu_M3_0};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M3_0)};

Mc.experiment(i).name  = tit;
Mc.experiment(i).size  = 1;
Mc.experiment(i).w     = {1};
Mc.experiment(i).mu    = {mu_cen};
Mc.experiment(i).sigma = {exp(esigma_cen)};

% Compile model
% (This generates the functional expression of parameters and derivatives.)
M = getMixtureModel(M,parameters.sym);
Mc = getMixtureModel(Mc,parameters.sym);

%% GET INITIAL ESTIMATE FOR MODEL PARAMETERS
% Loop: Experiments
for j = 1:length(D)
    % Construct data vector
    X = D{j}.data.uncensored(:);
    if j == 1
        % Perform estimation for this dataset
        gm = gmdistribution.fit(log(X),1);
        mu = gm.mu;
        esigma = log(sqrt(gm.Sigma));
        % Assign estimates
        parameters.guess([1:2]) = [mu;esigma];
    else
        % Perform estimation for this dataset
        try
            % Compute Gaussian mixture model for given dataset
            gm = gmdistribution.fit(log(X),2);
            % Assign parameters
            [mu,ind] = sort(gm.mu);
            esigma = log(sqrt(gm.Sigma(ind)));
            w = gm.PComponents(ind);
        catch
            gm = gmdistribution.fit(log(X),1);
            mu = gm.mu + [-1,+1];
            esigma = log(sqrt(gm.Sigma)*[1,1]);
            w = [0.5 0.5];
        end
        % Assign estimates
        parameters.guess([3*(j-1):3*j-1]) = [mu(1);esigma(1);w(1)];
    end
end

% Censoring times
X = [];
for j = 1:length(D)
    X = [X;D{j}.data.censored(:)];
end
% Perform estimation for this dataset
gm = gmdistribution.fit(log(X),1);
mu = gm.mu;
esigma = log(sqrt(gm.Sigma));
% Assign estimates
parameters.guess([end-1:end]) = [mu;esigma];

% Restrict to feasible set
parameters.guess = max(min(parameters.guess,parameters.max),parameters.min);

