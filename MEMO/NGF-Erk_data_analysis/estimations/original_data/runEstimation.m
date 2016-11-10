clc;
clear all;
close all;

%% OPTIONS
plot_opt = 'false'; % plots are shown
options.compute_P = 'true'; % enable computation of profiles

%% Model
% Select model
model_kinetic_dose_org_data; 

% Print model on screen
printModel(M);

% Process data for faster computations
D = processData(D);

% No model for right censored data since no right censored data present.
% Only needed for plot routines
Mc=[];

%% OPTIMIZATION
% Options
options_fit.n_starts = 500;
options_fit.plot = plot_opt;
options_fit.proposal = 'latin hypercube';
options_fit.fmincon = optimset('algorithm','interior-point',...%'active-set',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',4000,...
    'MaxFunEvals',4000*parameters.number);

% Run estimation
[parameters,M.fh.fit] = optimizeMultiStart(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_fit);

% Print result of estimation
printModel(M,parameters);

% Save optimization results
save('optimization','parameters','M','Mc','D','-v7.3');




%%  COMPUTE LIKELIHOOD PROFILES
if strcmp(options.compute_P,'true')
    % Options
    options_PL.plot = plot_opt;
    options_PL.parameter_index=[24:30];
    options_PL.plot = 'true';
    options_PL.P_next_step.min = 1e-1;
    options_PL.parameter_index=[24:30];
    
    % Compute profile likelihoods
    [parameters,M.fh.PL] = computeProfiles(parameters,@(theta,opt) logLikelihoodMM(theta,M,D,opt),options_PL);
end
save('profile','parameters','-v7.3');
opt.interval = 'static';

plotP(parameters,[],[24:30],opt);

%% VIZUALISATION OF RESULTS

% COMPARISION OF MIXTURE MODEL AND DATA 
% plot cdfs
options.plot_type='cdf';
[M.fh.model_data] = plotMixtureModel3(parameters,M,Mc,D);

% plot pdfs
options.width=18;
options.height= 12;

options.plot_type='pdf';
options.plot_data.type= 'hist';
options.plot_data.bins=50;
options.plot_subpop = 'true';


[M.fh.model_data] = plotMixtureModel_paper(parameters,M,Mc,D,options);
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 options.width options.height])
print('-depsc2','-r1000',['./figs/model_datafit_pdfs']);


%% PLOT MEAN RESPONSE OF POPULATION AND SUBPOPULATIONS
% specify width and height [cm] of plots
width=6;
height=6;


TextSizes.DefaultAxesFontSize = 8;
TextSizes.DefaultTextFontSize =8;
set(0,TextSizes);


%% Visualization -  kinetic response
e = 1;

% Simulation
t = linspace(DD(e).t(1),DD(e).t(end),100);
theta=parameters.MS.MAP.par;

% Calculation of model derived values
for k=1:2
    for i=1:5
        mu(k,i)=M.experiment(1,i).mu_fun{k}(theta);
        sigm(k,i)=M.experiment(1,i).sigma_fun{k}(theta);
        mean(k,i)= mu(k,i);
    end
end

mean_wei= mean(1,:).*(M.experiment(1,1).w_fun{1}(theta))+mean(2,:).*(M.experiment(1,1).w_fun{2}(theta));

% Plot
f1=figure('name','NGF kinetic - mean');
set(f1, 'Units','centimeters','Position', [5 5 width height])
indt = [1,1:length(t),length(t)];
t_indt = [-3,t,t(end)+3];



for k = 1:length(DD(e).t)
    plot(DD(e).t(k)*[1,1],nanmean(DD(e).Ey(k,:)) + nanstd(DD(e).Ey(k,:))*[-1,1],'-','color',[0.5,0.5 ,0.5],'linewidth',1);
    hold on
end

lh(1) = plot(DD(e).t,mean(2,:),'*','color',[0.75,0.75 ,0],'linewidth',1); hold on
lh(2) = plot(DD(e).t,mean(1,:),'*','color',[0.35,0.35 ,0],'linewidth',1); hold on
lh(3) = plot(DD(e).t,mean_wei,'x','color',[0,0.5 ,0],'linewidth',1); hold on
lh(4) = plot(DD(e).t,nanmean(DD(e).Ey,2),'o','color',[0.5,0.5 ,0.5],'linewidth',1);
hold on

% Legend / Label
xlabel('time {\itt} [min]');
ylabel(DD(1).measurand);
title('Kinetic for [NGF]_0 = 1 nM')

ylim([-0.2 1.2])

set(gca,'ytick',[-1:1],'yticklabel',[0.1,1,10]);
set(f1, 'Units','centimeters','Position', [5 5 width height])

% Save plot
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
print('-depsc2',['./figs/kin_resp']);

%% Visualization -  dose response
e = 2;
% Simulation
u = [0,10.^linspace(-3.1,1.1,21)];
u_plot = [10^-4,10.^linspace(-3.1,1.1,21)];
u_plot_data = 10.^[-4:1];
for k=1:2
    for i=6:11
        mu_dose(k,i-5)=M.experiment(1,i).mu_fun{k}(theta);
        sigm_dose(k,i-5)=M.experiment(1,i).sigma_fun{k}(theta);
        mean_dose(k,i-5)= mu_dose(k,i-5);
    end
end
mean_wei_dose= mean_dose(1,:).*(M.experiment(1,1).w_fun{1}(theta))+mean_dose(2,:).*(M.experiment(1,1).w_fun{2}(theta));
% Plot
f2=figure('name','NGF dose response - mean');
set(f2, 'Units','centimeters','Position', [5 5 width height])


plot(u_plot_data([1,2]).*[2.5,1/2],1.2*[1,1],'k:','markersize',8,'linewidth',1);hold on
for d = 1:length(DD(e).u)
    plot(u_plot_data(d)*[1,1],nanmean(DD(e).y(d,:)) + nanstd(DD(e).Ey(d,:))*[-1,1],'-','color',[0.5,0.5 ,0.5],'linewidth',1);hold on
end

lh(1) = plot(u_plot_data,mean_dose(2,:),'*','color',[0.75,0.75 ,0],'linewidth',1); hold on
lh(2) = plot(u_plot_data,mean_dose(1,:),'*','color',[0.35,0.35 ,0],'linewidth',1); hold on
lh(3) = plot(u_plot_data,mean_wei_dose,'x','color',[0,0.5 ,0],'linewidth',1); hold on
lh(4) = plot(u_plot_data(2:end),nanmean(DD(e).y(2:end,:),2),'o','color',[0.5,0.5 ,0.5],'linewidth',1);


xlim([u_plot(1)/1.5,u_plot(end)]);

% Legend / Label
legend(lh,{'model: mean s1','model: mean s2','model: mean','data: mean'},'location','northwest');
xlabel('NGF concentration');
ylabel(DD(1).measurand);
title('Dose reponse for t = 30 min')

set(gca,'ytick',[-1:1],'yticklabel',[0.1,1,10]);
set(gca,'xscale','log','xtick',10.^[-4:1],'xticklabel',[0,0.001,0.01,0.1,1,10]);

% Save plot
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
print('-depsc2',['./figs/dose_resp']);



