clc;
clear all;
close all;

%% OPTIONS
plot_opt ='false';% 'true'; % plots are shown


%% Model
model_Mad2_fun_mu_hill; % model with mu of strain specific subpopulation parameterized as hill function

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

%% OPTIMIZATION
% Options
options_fit.n_starts = 200;
options_fit.plot = plot_opt;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...
                           'display','off',...
                           'GradObj','on',...
                           'MaxIter',4000,...
                           'MaxFunEvals',4000*parameters.number);
% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMMc(theta,M,Mc,D,opt),options_fit);
% save
save('optimization','parameters','M','Mc','D','-v7.3');

% EVALUATION OF THE BIC
N = 0;
for j = 1:length(D)
    N = N + length(D{j}.data.censored) + length(D{j}.data.uncensored);
end
parameters.BIC  = -2*parameters.MS.MAP.logPost +   parameters.number*log(N);

% Print result of estimation
printModel(M,parameters);
printModel(Mc,parameters);




