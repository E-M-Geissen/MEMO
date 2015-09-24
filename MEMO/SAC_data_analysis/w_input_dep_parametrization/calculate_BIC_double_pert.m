% calculates the BIC of the double perturbation data, given the two model
% hyptheses as shown in manuscript Figure 6

clc;
close all
clear all
load optimization.mat

%%
Ki = parameters.MS.MAP.par(end-3);
Kj = parameters.MS.MAP.par(end-1);

ni = parameters.MS.MAP.par(end-2);
nj = parameters.MS.MAP.par(end);


xj_k=[0.3,0.6,1.1];

% Uncertainty of wild type fraction in double perturbations. Values from MCMC sampling of double perturbations
std_65_30=0.058048;
std_65_60=0.056668;
sig_65_120=0.067816;


sig=[std_65_30 std_65_60 sig_65_120];

w_meas=[0.22 0.66 0.86]; % Values from analysis of double perturbations


wfun_independent = @(xi,xj) (min((xi.^ni*(1+Ki^ni)./(Ki^ni+xi.^ni)),1+zeros(size(xi,1))).*min(((xj.^nj*(1+Kj^nj)./(Kj^nj+xj.^nj))),1+zeros(size(xj,1))));

wfun_dependent = @(xi,xj,a) (min((xi.^ni.*(1+(a*Ki./(xj+(a-1))).^ni)./((a*Ki./(xj+(a-1))).^ni+xi.^ni)),1+zeros(size(xi,1))) .* min((xj.^nj.*(1+(a*Kj./(xi+(a-1))).^nj))./((a*Kj./(xi+(a-1))).^nj+xj.^nj),1+zeros(size(xj,1))));

%% Calculate BIC asuming independent interaction
logL_in=0;
for i=1:length(w_meas)
    logL_in= logL_in+(log(1/(sqrt(2*pi)*sig(i)))-0.5*((w_meas(i)-wfun_independent(0.71,xj_k(i)))^2/sig(i)^2));
end

BIC_in_dep=-2*logL_in+0*log(3)

%% Find optimal additional parameter / BIC for interacting effects using a line search
ind =1;

% calculate BIC for additional parameter m from 1 to 10
for m=1:0.05:10
    aa(ind)=m;
    logL_de=0;
    for i=1:length(w_meas)
        logL_de= logL_de+(log(1/(sqrt(2*pi)*sig(i)))-0.5*((w_meas(i)-wfun_dependent(0.71,xj_k(i),m))^2/sig(i)^2));
    end
    Log_m(ind)=logL_de;
    BIC_ind(ind)=-2*logL_de+1*log(3);
    ind=ind+1;
end

[C,I] = min(BIC_ind)% optimal BIC
aa(I) % optimal parameter
