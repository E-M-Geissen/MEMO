clc;
close all
clear all
load optimization.mat

%% set figure size
width=6;
height=5;

TextSizes.DefaultAxesFontSize = 10;
TextSizes.DefaultTextFontSize =10;
set(0,TextSizes);
%% parameter values from optimization

xi = 0:0.01:1; % Mad2
xj = 0:0.01:0.95; % Mad3

mu_null = parameters.MS.MAP.par(find(parameters.sym=='mu_null'));
KM2 = parameters.MS.MAP.par(find(parameters.sym=='KM2'));
nh = parameters.MS.MAP.par(find(parameters.sym=='nh'));
vmax = parameters.MS.MAP.par(find(parameters.sym=='vmax'));



e_M2_80=parameters.MS.MAP.par(find(parameters.sym=='e_M2_80'));
e_M2_65_P50=parameters.MS.MAP.par(find(parameters.sym=='e_M2_65_P50'));
e_M2_65_P188=parameters.MS.MAP.par(find(parameters.sym=='e_M2_65_P188'));
e_M2_40=parameters.MS.MAP.par(find(parameters.sym=='e_M2_40'));
e_M2_10=parameters.MS.MAP.par(find(parameters.sym=='e_M2_10'));
e_M2_20=parameters.MS.MAP.par(find(parameters.sym=='e_M2_20'));

e_M3_30=parameters.MS.MAP.par(find(parameters.sym=='e_M3_30'));
e_M3_60=parameters.MS.MAP.par(find(parameters.sym=='e_M3_60'));

si_0=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_0')));
si_10=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_10')));
si_20=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_20')));
si_40=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_40')));
si_65_P50=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_65_P50')));
si_65_P188=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_65_P188')));
si_80=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_M2_80')));
si_WT=exp(parameters.MS.MAP.par(find(parameters.sym=='esigma_WT')));



mu_WT=parameters.MS.MAP.par(find(parameters.sym=='mu_WT'));

expname={M.experiment.name};
w10=M.experiment(find(ismember(expname,'10% Mad2'))).w_fun{1}(parameters.MS.MAP.par);
w20=M.experiment(find(ismember(expname,'20% Mad2'))).w_fun{1}(parameters.MS.MAP.par);
w40=M.experiment(find(ismember(expname,'40% Mad2'))).w_fun{1}(parameters.MS.MAP.par);
w65_P50=M.experiment(find(ismember(expname,'65% Mad2 P50'))).w_fun{1}(parameters.MS.MAP.par);
w65_P188=M.experiment(find(ismember(expname,'65% Mad2 P188'))).w_fun{1}(parameters.MS.MAP.par);
w80=M.experiment(find(ismember(expname,'80% Mad2'))).w_fun{1}(parameters.MS.MAP.par);
%w_WT=M.experiment(find(ismember(expname,'WT fusion'))).w_fun{1}(parameters.MS.MAP.par);
% w30=M.experiment(find(ismember(expname,'30% Mad3 (2)'))).w_fun{1}(parameters.MS.MAP.par);
% w60=M.experiment(find(ismember(expname,'60% Mad3'))).w_fun{1}(parameters.MS.MAP.par);
% w120=M.experiment(find(ismember(expname,'120% Mad3'))).w_fun{1}(parameters.MS.MAP.par);


 
 wfun_mu = @(xi) (mu_null+vmax.*(xi).^nh./(KM2^nh+(xi).^nh)); 
 wfun_mux = @(xx) (mu_null+vmax.*(xx).^nh./(KM2^nh+(xx).^nh)); 
w_fun_mean=@(mu,sig)(exp(mu+((sig.^2)./2)));
%% fit
%Mad2
hFig=figure;
set(hFig, 'Units','centimeters','Position', [10 10 width height])

set(gca, 'LooseInset', [0,0.1,0,0]);
xlabel(['Mad2'])
ylabel(['\mu'])
hold on
box off
x_ar=[0, 0.10+e_M2_10,0.20+e_M2_20,0.40+e_M2_40,0.65+e_M2_65_P50,0.65+e_M2_65_P188,0.8+e_M2_80];
plot([0,0.10,0.20,0.40,0.65,0.65,0.8], wfun_mux(x_ar),'o','MarkerSize',6,'Color',[0,0,0],'LineWidth',1.5);

plot(xi,wfun_mu(xi),'LineWidth',1.5); hold on
plot([0.10+e_M2_10,0.20+e_M2_20,0.40+e_M2_40,0.65+e_M2_65_P50,0.65+e_M2_65_P188,0.8+e_M2_80], wfun_mux(x_ar(2:end)),'x','MarkerSize',6,'Color',[87/255,119/255,1],'LineWidth',1.5);
hold on
%plot(1, mu_WT, 'rx','MarkerSize',10);

set(gca,'Xtick', [0  0.5 1],'Xticklabel',{'0','50','100'})
% set(gca,'Ytick', [0  0.5 1],'Yticklabel',{'0%','50%','100%'},'Layer','bottom')
 %xlim([0 1])
 xlabel(['%Mad2'])
% ylim([0 1.01])

set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])
print('-depsc',['figs/mu(u)']);
%Mad2
hFig1=figure;
mus=wfun_mux(x_ar);
sigs=[si_0,si_10,si_20,si_40,si_65_P50,si_65_P188,si_80];

si_WT=parameters.MS.MAP.par(find(parameters.sym=='esigma_WT'));
set(hFig1, 'Units','centimeters','Position', [10 10 width height])
plot(x_ar, w_fun_mean(mus,sigs),'x');hold on

%plot(1,(exp(mu_WT+((si_WT.^2)./2))),'rx');
ylim([0,340])
%ylabel([',WT fraction'])
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 width height])




