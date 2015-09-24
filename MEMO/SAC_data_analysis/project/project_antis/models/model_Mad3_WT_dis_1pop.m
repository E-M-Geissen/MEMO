% Data:
% Mad2 delta 
% Model:
%   log-normal
%   
%  two components variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_M3_WT_1 esigma_M3_WT_1  ...
        gamma_cen  esigma_cen     elambda_cen     exi_cen;
parameters.sym  = [ mu_M3_WT_1 ; esigma_M3_WT_1 ;...
                    gamma_cen  ; esigma_cen     ; elambda_cen ; exi_cen];

parameters.name = {'mu_1_{M3,WT}' ; 'esigma_1_{M3,WT}' ; ...
                   'gamma_{cen}'  ; 'esigma_{cen}'     ; 'elambda_{cen}' ; 'exi_{cen}' };

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-2); ...
                     -20   ; log(1e-4); log(5); log(5)];      
% As lower bound for the and mean we used the
% inter-observation time.
parameters.max  = [ 
                    log(2e3); log(1e1); ...
                          20; log(1e4); log(2e3); log(2e3)];
% As lower bound for the median we used twice
% the observation time.
parameters.init_fun = @(theta_0,theta_min,theta_max,n) ...
            bsxfun(@max,bsxfun(@min,...
                        bsxfun(@plus,theta_0,[randn(2,n);0.2*randn(1,n);...
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
data_WT_Mad3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 5;

M.experiment(i).name  = tit;
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).mu    = {mu_M3_WT_1};
M.experiment(i).sigma = {exp(esigma_M3_WT_1)};

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