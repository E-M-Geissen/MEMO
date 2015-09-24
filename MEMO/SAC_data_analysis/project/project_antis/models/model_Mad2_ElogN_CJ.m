% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments

%% SPECIFY MODEL PARAMETERS
syms mu_c       esigma_c                 ...
    mu_M2_80   esigma_M2_80   w_M2_80   ...
    mu_M2_60_P50   esigma_M2_60_P50   w_M2_60_P50   ...
    mu_M2_60_P188   esigma_M2_60_P188    w_M2_60_P188   ...
    mu_M2_40   esigma_M2_40   w_M2_40   ...
    mu_M2_0 esigma_M2_0  w_M2_0 ...
    gamma_cen  esigma_cen     elambda_cen     exi_cen;


parameters.sym  = [ mu_c       ; esigma_c       ; ...
    mu_M2_80   ; esigma_M2_80   ; w_M2_80   ; ...
    mu_M2_60_P50   ; esigma_M2_60_P50   ; w_M2_60_P50   ; ...
    mu_M2_60_P188    ; esigma_M2_60_P188    ; w_M2_60_P188   ; ...
    mu_M2_40   ; esigma_M2_40   ; w_M2_40   ; ...
    mu_M2_0; esigma_M2_0;  w_M2_0; ...
    gamma_cen  ; esigma_cen     ; elambda_cen ; exi_cen];


parameters.name = {'\mu_c'         ; 'log(\sigma_c)'         ;                 ...
    '\mu_{M2,80}'   ; 'log(\sigma_{M2,80})'   ; 'w_{M2,80}'   ; ...
    '\mu_{M2,60_P50}'   ; 'log(\sigma_{M2,60_P50})'   ; 'w_{M2,60_P50}'   ; ...
    '\mu_{M2,60_P188}'   ; 'log(\sigma_{M2,60_P188})'   ; 'w_{M2,60_P188}'   ; ...
    '\mu_{M2,40}'   ; 'log(\sigma_{M2,40})'   ; 'w_{M2,40}'   ; ...
    '\mu_{M2,0}'    ; 'log(\sigma_{M2,0})'    ; 'w_{M2,0}'    ; ...
    '\gamma_{cen}'  ; 'log(\sigma_{cen})'     ; 'log(\lambda_{cen})' ; 'log(\xi_{cen})' };


parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-1);   ...
log(5); log(1e-1); 0; ...
    log(5); log(1e-1); 0; ...
    log(5); log(1e-1); 0; ...
    log(5); log(1e-1); 0; ...
    log(5); log(1e-1); 0; ...
    -20   ; log(1e-4); log(5); log(5)];

parameters.max  = [ log(2e3); log(1e1);    ...
    log(2e3); log(1e1); 1; ...
    log(2e3); log(1e1); 1; ...
    log(2e3); log(1e1); 1; ...
    log(2e3); log(1e1); 1; ...
    log(2e3); log(1e1); 1; ...
    20; log(1e4); log(2e3); log(2e3)];

parameters.init_fun = @(theta_0,theta_min,theta_max,n) ...
    bsxfun(@max,bsxfun(@min,...
    bsxfun(@plus,theta_0,[...
    randn(2,n);...
    randn(2,n);0.2*randn(1,n);...
    randn(2,n);0.2*randn(1,n);...
    randn(2,n);0.2*randn(1,n);...
    randn(2,n);0.2*randn(1,n);...
    randn(2,n);0.2*randn(1,n);...
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
M.experiment(i).mu    = {mu_c,mu_M2_40};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_40)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
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
M.experiment(i).w     = {1-w_M2_60_P50,w_M2_60_P50};
M.experiment(i).mu    = {mu_c,mu_M2_60_P50};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_60_P50)};

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
M.experiment(i).w     = {1-w_M2_60_P188,w_M2_60_P188};
M.experiment(i).mu    = {mu_c,mu_M2_60_P188};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_60_P188)};

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
M.experiment(i).mu    = {mu_c,mu_M2_80};
M.experiment(i).sigma = {exp(esigma_c),exp(esigma_M2_80)};

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

% %% GET INITIAL ESTIMATE FOR MODEL PARAMETERS
% % Loop: Experiments
% for j = 1:length(D)
%     % Construct data vector
%     X = D{j}.data.uncensored(:);
%     if j == 1
%         % Perform estimation for this dataset
%         gm = gmdistribution.fit(log(X),1);
%         mu = gm.mu;
%         esigma = log(sqrt(gm.Sigma));
%         % Assign estimates
%         parameters.guess([1:2]) = [mu;esigma];
%     else
%         % Perform estimation for this dataset
%         try
%             % Compute Gaussian mixture model for given dataset
%             gm = gmdistribution.fit(log(X),2);
%             % Assign parameters
%             [mu,ind] = sort(gm.mu);
%             esigma = log(sqrt(gm.Sigma(ind)));
%             w = gm.PComponents(ind);
%         catch
%             gm = gmdistribution.fit(log(X),1);
%             mu = gm.mu + [-1,+1];
%             esigma = log(sqrt(gm.Sigma)*[1,1]);
%             w = [0.5 0.5];
%         end
%         % Assign estimates
%         parameters.guess([3*(j-1):3*j-1]) = [mu(1);esigma(1);w(1)];
%     end
% end
% 
% % Censored data
% Xc = [];
% for j = 1:length(D)
%     % Construct data vector
%     Xc = [Xc;D{j}.data.censored(:)];
% end
% % Perform estimation for this dataset
% gm = gmdistribution.fit(Xc,1);
% mu = gm.mu;
% sigma = sqrt(gm.Sigma);
% % Assign estimates
% parameters.guess(end-3:end) = [0;log(10);log(10*sigma);log(mu)];

% Restrict to feasible set
parameters.guess = max(min(parameters.guess,parameters.max),parameters.min);

