
%  two components variable across experiments

%% SPECIFY MODEL PARAMETERS 
syms mu_batf_1 log_sigma_batf_1 mu_batf_2 log_sigma_batf_2  ...
     mu_gata3_1 log_sigma_gata3_1 mu_gata3_2 log_sigma_gata3_2 ...
     mu_il4ra_1 log_sigma_il4ra_1 mu_il4ra_2 log_sigma_il4ra_2 ...
     mu_il2ra_1 log_sigma_il2ra_1 mu_il2ra_2 log_sigma_il2ra_2 w;

parameters.sym  = [ mu_batf_1; log_sigma_batf_1; mu_batf_2; log_sigma_batf_2;  ...
     mu_gata3_1; log_sigma_gata3_1; mu_gata3_2; log_sigma_gata3_2; ...
     mu_il4ra_1; log_sigma_il4ra_1; mu_il4ra_2; log_sigma_il4ra_2; ...
     mu_il2ra_1; log_sigma_il2ra_1; mu_il2ra_2; log_sigma_il2ra_2; w];

parameters.name = {'\mu_{batf_1}'; 'log(\sigma_{batf_1)}'; '\mu_{batf_2}'; 'log(\sigma_{batf_2})';  ...
     '\mu_{gata3_1}'; 'log(\sigma_{gata3_1)}'; '\mu_{gata3_2}'; 'log(\sigma_{gata3_2)}'; ...
     '\mu_{il4ra_1}'; 'log(\sigma_{il4ra_1)}'; '\mu_{il4ra_2}'; 'log(\sigma_{il4ra_2)}'; ...
     '\mu_{il2ra_1}'; 'log(\sigma_{il2ra_1)}'; '\mu_{il2ra_2}'; 'log(\sigma_{il2ra_2)}'; 'w' };

parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
% parameters.min  = [  log(2); log(1e-2);log(2); log(1e-2); 
%     log(2); log(1e-2);log(2); log(1e-2);
%     log(2); log(1e-2);log(2); log(1e-2);
%     log(2); log(1e-2);log(2); log(1e-2);0 ];      
% % As lower bound for the and mean we used the
% % inter-observation time.
% parameters.max  = [ log(5); log(1e1);log(5); log(1e1);
%     log(5); log(1e1);log(5); log(1e1);
%     log(5); log(1e1);log(5); log(1e1);
%                     log(5); log(1e1);log(5); log(1e1);1];
parameters.min  = [  2; -2;2; -2; 
    2; -2;2; -2;
    2; -2;2; -2;
    2; -2;2; -2;0 ];      
% As lower bound for the and mean we used the
% inter-observation time.
parameters.max  = [ 15; 5;15; 5;
    15; 5;15; 5;
    15; 5;15; 5;
                    10; 5;15; 5;1];

%% SPECIFY MODEL AND DATA
M.mixture.type = 'log-normal';
M.label.x = 'counts';
M.label.y = 'probability density';



% EXPERIMENT
i = 1;
data_Batf;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 1;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w,w};
M.experiment(i).mu    = {mu_batf_1,mu_batf_2};
M.experiment(i).sigma = {exp(log_sigma_batf_1),exp(log_sigma_batf_2)};




i=i+1;
data_Gata3;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 1;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w,w};
M.experiment(i).mu    = {mu_gata3_1,mu_gata3_2};
M.experiment(i).sigma = {exp(log_sigma_gata3_1),exp(log_sigma_gata3_2)};

i=i+1;
data_Il4ra;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 1;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w,w};
M.experiment(i).mu    = {mu_il4ra_1,mu_il4ra_2};
M.experiment(i).sigma = {exp(log_sigma_il4ra_1),exp(log_sigma_il4ra_2)};

i=i+1;
data_Il2ra;

D{i}.name = tit;
D{i}.description = [];
D{i}.data.uncensored = Tm;
D{i}.data.censored = Tc;
D{i}.observation_interval = 1;

M.experiment(i).name  = tit;
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w,w};
M.experiment(i).mu    = {mu_il2ra_1,mu_il2ra_2};
M.experiment(i).sigma = {exp(log_sigma_il2ra_1),exp(log_sigma_il2ra_2)};

% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);
