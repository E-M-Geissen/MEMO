clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'true'; % plots are shown
options.compute_P = 'true'; % computation of profiles

%% Model
%model_kinetics_fun_mu; %
model_kinetics_fun_mu_new_param_log10_1024; %

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

Mc=[];
%% OPTIMIZATION
% Options
options_fit.n_starts = 200;
options_fit.plot = plot_opt;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...%'active-set',...
                           'display','off',...
                           'GradObj','on',...
                           'MaxIter',4000,...
                           'MaxFunEvals',4000*parameters.number);
%

%parameters.guess= [-1.31224933892433;-0.800170025950256;0.308686895511533;-1.00664990398034;-0.793472192038934;-0.958362694943307;-0.566965993149609;-0.880727397693913;-0.782505437418861;-0.844955577642641;-0.659188864221808;-0.893688851987427;-0.609383901991252;-0.805528955667771;-0.336775651990680;-0.904264595011374;-0.402351929567410;-0.849032552318242;-0.533867263234241;-0.924174324473709;-0.677188465471327;-1.01589979044468;-0.735719735078803;0.000396284244769350;2.03477266321198;-6.17071496092370;3.13002830013107;-4.52353825267214;-1.64103414327923;2.89605292538153];
parameters.guess=[-2.57375975785411;-4.60515965242695;0.415541633166493;-2.61370522674492;-2.26711101021480;-4.60516322191122;-2.00003759813577;-2.60171587128831;-2.22324547246141;-2.53141041418677;-2.03046529058666;-4.23242127745770;-2.28401141489884;-2.89532106867208;-2.01534525927953;-3.44303737782804;-2.05934014825201;-2.77975809305127;-1.99622716680515;-2.69195469480857;-2.06858578380835;-4.03041868213555;-2.19340645741587;0.0102782091322202;1.99334956026027;-5.78970679738279;3.70860227278220;-2.15487002218453;-2.56965328120499;0.503346753385896];
% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_fit);
% Print result of estimation
printModel(M,parameters);

%%
% save
save('optimization','parameters','M','Mc','D','-v7.3');
