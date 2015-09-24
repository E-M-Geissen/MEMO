% Data:
%   control
%   all Mad2 data (no double mutants)
% Model:
%   log-normal
%   one component fixed across experiments
%   one component variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms gamma_c       esigma_c       elambda_c       exi_c                 ...
     gamma_M1_30   esigma_M1_30   elambda_M1_30   exi_M1_30   w_M1_30   ...
     gamma_M1_10   esigma_M1_10   elambda_M1_10   exi_M1_10   w_M1_10   ...
     gamma_M2_80   esigma_M2_80   elambda_M2_80   exi_M2_80   w_M2_80   ...
     gamma_M2_60   esigma_M2_60   elambda_M2_60   exi_M2_60   w_M2_60   ...
     gamma_M2_25   esigma_M2_25   elambda_M2_25   exi_M2_25   w_M2_25   ...
     gamma_M2_0    esigma_M2_0    elambda_M2_0    exi_M2_0    w_M2_0    ...
     gamma_M3_60   esigma_M3_60   elambda_M3_60   exi_M3_60   w_M3_60   ...
     gamma_M3_30_1 esigma_M3_30_1 elambda_M3_30_1 exi_M3_30_1 w_M3_30_1 ...
     gamma_M3_30_2 esigma_M3_30_2 elambda_M3_30_2 exi_M3_30_2 w_M3_30_2 ...
     gamma_M3_0    esigma_M3_0    elambda_M3_0    exi_M3_0    w_M3_0    ...
     gamma_cen     esigma_cen     elambda_cen     exi_cen;
parameters.sym  = [  gamma_c       ; esigma_c       ; elambda_c       ; exi_c       ;             ...
                     gamma_M1_30   ; esigma_M1_30   ; elambda_M1_30   ; exi_M1_30   ; w_M1_30   ; ...
                     gamma_M1_10   ; esigma_M1_10   ; elambda_M1_10   ; exi_M1_10   ; w_M1_10   ; ...
                     gamma_M2_80   ; esigma_M2_80   ; elambda_M2_80   ; exi_M2_80   ; w_M2_80   ; ...
                     gamma_M2_60   ; esigma_M2_60   ; elambda_M2_60   ; exi_M2_60   ; w_M2_60   ; ...
                     gamma_M2_25   ; esigma_M2_25   ; elambda_M2_25   ; exi_M2_25   ; w_M2_25   ; ...
                     gamma_M2_0    ; esigma_M2_0    ; elambda_M2_0    ; exi_M2_0    ; w_M2_0    ; ...
                     gamma_M3_60   ; esigma_M3_60   ; elambda_M3_60   ; exi_M3_60   ; w_M3_60   ; ...
                     gamma_M3_30_1 ; esigma_M3_30_1 ; elambda_M3_30_1 ; exi_M3_30_1 ; w_M3_30_1 ; ...
                     gamma_M3_30_2 ; esigma_M3_30_2 ; elambda_M3_30_2 ; exi_M3_30_2 ; w_M3_30_2 ; ...
                     gamma_M3_0    ; esigma_M3_0    ; elambda_M3_0    ; exi_M3_0    ; w_M3_0    ; ...
                     gamma_cen     ; esigma_cen     ; elambda_cen     ; exi_cen];
parameters.name = { 'gamma_c'         ; 'esigma_c'         ; 'elambda_c'         ; 'exi_c'         ;                 ...
                    'gamma_{M1,30}'   ; 'esigma_{M1,30}'   ; 'elambda_{M1,30}'   ; 'exi_{M1,30}'   ; 'w_{M1,30}'   ; ...
                    'gamma_{M1,10}'   ; 'esigma_{M1,10}'   ; 'elambda_{M1,10}'   ; 'exi_{M1,10}'   ; 'w_{M1,10}'   ; ...
                    'gamma_{M2,80}'   ; 'esigma_{M2,80}'   ; 'elambda_{M2,80}'   ; 'exi_{M2,80}'   ; 'w_{M2,80}'   ; ...
                    'gamma_{M2,60}'   ; 'esigma_{M2,60}'   ; 'elambda_{M2,60}'   ; 'exi_{M2,60}'   ; 'w_{M2,60}'   ; ...
                    'gamma_{M2,25}'   ; 'esigma_{M2,25}'   ; 'elambda_{M2,25}'   ; 'exi_{M2,25}'   ; 'w_{M2,25}'   ; ...
                    'gamma_{M2,0}'    ; 'esigma_{M2,0} '   ; 'elambda_{M2,0}'    ; 'exi_{M2,0}'    ; 'w_{M2,0}'    ; ...
                    'gamma_{M3,60}'   ; 'esigma_{M3,60}'   ; 'elambda_{M3,60}'   ; 'exi_{M3,60}'   ; 'w_{M3,60}'   ; ...
                    'gamma_{M3,30,1}' ; 'esigma_{M3,30,1}' ; 'elambda_{M3,30,1}' ; 'exi_{M3,30,1}' ; 'w_{M3,30,1}' ; ...
                    'gamma_{M3,30,2}' ; 'esigma_{M3,30,2}' ; 'elambda_{M3,30,2}' ; 'exi_{M3,30,2}' ; 'w_{M3,30,2}' ; ...
                    'gamma_{M3,0}'    ; 'esigma_{M3,0}'    ; 'elambda_{M3,0}'    ; 'exi_{M3,0}'    ; 'w_{M3,0}'    ; ...
                    'gamma_{cen}'     ; 'esigma_{cen}'     ; 'elambda_{cen}'     ; 'exi_{cen}'   };
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  -20; log(1e-4);    log(5);    log(5);    ...
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5); 0; ...       
                     -20; log(1e-4);    log(5);    log(5)];      
                % As lower bound for the variance and mean we used the
                % inter-observation time.
parameters.max  = [   20;  log(1e4);  log(2e3);  log(2e3);    ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3); 1; ...
                      20;  log(1e4);  log(2e3);  log(2e3)];
                % As lower bound for the variance and mean we used twice
                % the observation time.
parameters.init_fun = @(theta_0,theta_min,theta_max) ...
            max(min( theta_0 + [...
                                    2*randn(4,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);0.3*randn(1,1);...
                                    2*randn(4,1);...
                               ],theta_max),theta_min);

%% SPECIFY MODEL AND DATA
M.mixture.type = 'Johnson SU';
M.label.x = 'time [min]';
M.label.y = 'probability density';

Mc.mixture.type = 'Johnson SU';
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

M.experiment(i).name   = tit;
M.experiment(i).size   = 1;
M.experiment(i).w      = {            1 };
M.experiment(i).gamma  = {      gamma_c };
M.experiment(i).sigma  = { exp(esigma_c)};
M.experiment(i).lambda = {exp(elambda_c)};
M.experiment(i).xi     = {    exp(exi_c)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad1;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {  1 - w_M1_30 ,           w_M1_30 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M1_30 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M1_30)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M1_30)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M1_30)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_10_Mad1;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {  1 - w_M1_10 ,           w_M1_10 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M1_10 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M1_10)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M1_10)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M1_10)};

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

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {  1 - w_M2_80 ,           w_M2_80 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M2_80 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M2_80)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M2_80)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M2_80)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_60_Mad2_new;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {  1 - w_M2_60 ,           w_M2_60 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M2_60 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M2_60)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M2_60)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M2_60)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_25_Mad2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {  1 - w_M2_25 ,           w_M2_25 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M2_25 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M2_25)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M2_25)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M2_25)};

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

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {   1 - w_M2_0 ,           w_M2_0 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M2_0 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M2_0)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M2_0)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M2_0)};

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

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {  1 - w_M3_60 ,           w_M3_60 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M3_60 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M3_60)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M3_60)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M3_60)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad3_1;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {1 - w_M3_30_1 ,           w_M3_30_1 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M3_30_1 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M3_30_1)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M3_30_1)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M3_30_1)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

% EXPERIMENT
i = i + 1;
data_30_Mad3_2;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {1 - w_M3_30_2 ,           w_M3_30_2 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M3_30_2 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M3_30_2)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M3_30_2)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M3_30_2)};

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

M.experiment(i).name   = tit;
M.experiment(i).size   = 2;
M.experiment(i).w      = {   1 - w_M3_0 ,           w_M3_0 };
M.experiment(i).gamma  = {      gamma_c ,       gamma_M3_0 };
M.experiment(i).sigma  = { exp(esigma_c),  exp(esigma_M3_0)};
M.experiment(i).lambda = {exp(elambda_c), exp(elambda_M3_0)};
M.experiment(i).xi     = {    exp(exi_c),     exp(exi_M3_0)};

Mc.experiment(i).name   = tit;
Mc.experiment(i).size   = 1;
Mc.experiment(i).w      = {            1 };
Mc.experiment(i).gamma  = {      gamma_cen };
Mc.experiment(i).sigma  = { exp(esigma_cen)};
Mc.experiment(i).lambda = {exp(elambda_cen)};
Mc.experiment(i).xi     = {    exp(exi_cen)};

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
        gm = gmdistribution.fit(X,1);
        mu = gm.mu;
        sigma = sqrt(gm.Sigma);
        % Assign estimates
        parameters.guess([1:4]) = [0;log(10);log(10*sigma);log(mu)];
    else
        % Perform estimation for this dataset
        try
            % Compute Gaussian mixture model for given dataset
            gm = gmdistribution.fit(X,2);
            % Assign parameters
            [mu,ind] = sort(gm.mu);
            sigma = sqrt(gm.Sigma(ind));
            w = gm.PComponents(ind);
        catch
            gm = gmdistribution.fit(X,1);
            mu = gm.mu + [-1,+1];
            esigma = log(sqrt(gm.Sigma)*[1,1]);
            w = [0.5 0.5];
        end
        % Assign estimates
        parameters.guess([5*(j-1):5*j-1]) = [0;log(10);log(10*sigma(1));log(mu(1));w(1)];
    end
end

% Censored data
Xc = [];
for j = 1:length(D)
    % Construct data vector
    Xc = [Xc;D{j}.data.censored(:)];
end
% Perform estimation for this dataset
gm = gmdistribution.fit(Xc,1);
mu = gm.mu;
sigma = sqrt(gm.Sigma);
% Assign estimates
parameters.guess(end-3:end) = [0;log(10);log(10*sigma);log(mu)];

% Restrict to feasible set
parameters.guess = max(min(parameters.guess,parameters.max),parameters.min);

