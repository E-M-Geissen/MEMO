% Data:
% Mad2 delta 
% Model:
%   log-normal
%   
%  two components variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_M2_0_1 esigma_M2_0_1 ;

parameters.sym  = [ mu_M2_0_1 ; esigma_M2_0_1 ];

parameters.name = {'mu_1_{M2,0}' ; 'esigma_1_{M2,0}' };

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(5); log(1e-2)];      
                % As lower bound for the and mean we used the
                % inter-observation time.
parameters.max  = [ 
                    log(2e3); log(1e1)];
        % As lower bound for the median we used twice
                % the observation time.
% parameters.init_fun = @(theta_0,theta_min,theta_max,n) ...
%             bsxfun(@max,bsxfun(@min,...
%                         bsxfun(@plus,theta_0,[randn(4,n);0.2*randn(1,n)]),theta_max),theta_min);

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
M.experiment(i).size  = 1;
M.experiment(i).w     = {1};
M.experiment(i).mu    = {mu_M2_0_1};
M.experiment(i).sigma = {exp(esigma_M2_0_1)};




% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);

% %% GET INITIAL ESTIMATE FOR MODEL PARAMETERS
% % Loop: Experiments
% for j = 1:length(D)
%     % Construct data vector
%     X = [D{j}.data.censored(:);D{j}.data.uncensored(:)];
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
%         parameters.guess(2+(j-2)*5+[1:5]) = [mu(2);esigma(2); 
%                                             mu(1);esigma(1);w(1)];
%     end
% end
% 
% 
% 
% % Perform estimation for this dataset
% gm = gmdistribution.fit(Xc,1);
% mu = gm.mu;
% sigma = sqrt(gm.Sigma);
% 
% % Assign estimates
% parameters.guess(end-3:end) = [0;log(10);log(10*sigma);log(mu)];
% 
% % Restrict to feasible set
% parameters.guess = max(min(parameters.guess,parameters.max),parameters.min);
