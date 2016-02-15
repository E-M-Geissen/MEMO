% Data:
%   natural logarithm of original data
% Model:
% new model normal 

%% SPECIFY MODEL PARAMETERS 
syms        log_sigma_0_1         log_sigma_0_2       w      ...
           log_sigma_5_1        log_sigma_5_2           ...
           log_sigma_15_1         log_sigma_15_2           ...
           log_sigma_30_1        log_sigma_30_2          ...
         log_sigma_60_1          log_sigma_60_2 ...
           log_sigma_30_0_1       log_sigma_30_0_2...
          log_sigma_30_0001_1       log_sigma_30_0001_2...
          log_sigma_30_001_1       log_sigma_30_001_2...
          log_sigma_30_01_1       log_sigma_30_01_2...
           log_sigma_30_1_1       log_sigma_30_1_2...
           log_sigma_30_10_1        log_sigma_30_10_2 ...    
      k2 KD k3T B k4 k5 SE  t ;
parameters.sym  =  [     log_sigma_0_1 ;    log_sigma_0_2 ;      w  ;   ...
          log_sigma_5_1 ;       log_sigma_5_2 ;           ...
          log_sigma_15_1 ;      log_sigma_15_2  ;      ...
         log_sigma_30_1 ;       log_sigma_30_2 ;        ...
        log_sigma_60_1 ;       log_sigma_60_2  
            log_sigma_30_0_1;    log_sigma_30_0_2;...
          log_sigma_30_0001_1;       log_sigma_30_0001_2;...
          log_sigma_30_001_1;      log_sigma_30_001_2;...
          log_sigma_30_01_1;      log_sigma_30_01_2;...
         log_sigma_30_1_1;    log_sigma_30_1_2;...
          log_sigma_30_10_1;      log_sigma_30_10_2; ...
        k2; KD; k3T;B; k4 ;k5; SE];
parameters.name = {     'log(\sigma_{0_1})' ;       'log(\sigma_{0_2})' ;      'w'   ;   ...
       'log(\sigma_{5_1})' ;       'log(\sigma_{5_2})' ;         ...
        'log(\sigma_{15_1})' ;        'log(\sigma_{15_2})' ;        ...
        'log(\sigma_{30_1})' ;       'log(\sigma_{30_2})' ;        ...
        'log(\sigma_{60_1})' ;       'log(\sigma_{60_2})' ; ...
          'log(\sigma_{30_0_1})' ;         'log(\sigma_{30_0_2})'; ... 
                'log(\sigma_{30_0001_1})' ;        'log(\sigma_{30_0001_2})'; ...
                   'log(\sigma_{30_001_1})' ;       'log(\sigma_{30_001_2})'; ...
                            'log(\sigma_{30_01_1})' ;       'log(\sigma_{30_01_2})'; ...
                                'log(\sigma_{30_1_1})' ;        'log(\sigma_{30_1_2})'; ...
                                         'log(\sigma_{30_10_1})' ;       'log(\sigma_{30_10_2})';...
            'k_2'; 'KD'; 'k3T'; 'B';'k_4' ;'k_5'; 'SE'};
parameters.number = length(parameters.sym);
parameters.guess = zeros(parameters.number,1);
parameters.min  = [  log(1e-2);   log(1e-2);0;  ...
                      log(1e-2);  log(1e-2); ...       
                     log(1e-2);log(1e-2); ...       
                     log(1e-2); log(1e-2);...       
                      log(1e-2);  log(1e-2);...
                     log(1e-2);   log(1e-2);  ...
                      log(1e-2); log(1e-2); ...       
                     log(1e-2); log(1e-2); ...       
                      log(1e-2); log(1e-2);...     
                     log(1e-2);  log(1e-2);...     
                      log(1e-2);  log(1e-2); ...
                     -10;-10; -10; 0; -10; -10; -10];    
                 

parameters.max  = [  log(1e1);log(1e1); 1;   ...
                     log(1e1);log(1e1);  ...
                     log(1e1);log(1e1);   ...
                     log(1e1);log(1e1);   ...
                     log(1e1);log(1e1);  ...
                    log(1e1);log(1e1);    ...
                    log(1e1);log(1e1);     ...
                    log(1e1);log(1e1); ...
                    log(1e1);log(1e1);   ...
                   log(1e1);log(1e1);   ...
                    log(1e1);log(1e1); ...
                    10;10; 10; 10; 10; 10; 10];

%%
bins=256;
upper=2;
lower=-1;

cen=(upper-lower)/bins;
%% Measurement data
% Load data and asign to struct
ExpC = load_NGF_kinetic;
DD(1).name = 'NGF kinetic';
DD(1).type = 'kinetic';
DD(1).measurand = 'relative pErk concentration';
DD(1).t = [0,5,15,30,60];
DD(1).u = 1;
DD(1).y = nan(length(DD(1).t),6000);

DD(1).Ey = nan(length(DD(1).t),10);
s = [0.7978,0.9586,1.1317,1.1119]; % scaling computed in data_scaling.m
for k = 1:length(DD(1).t)
    n0 = 0;
    for r = 1:length(ExpC(k).experiment)
        nx = size(ExpC(k).experiment(r).data,1);
        DD(1).y(k,n0+1:n0+nx) = s(r)*ExpC(k).experiment(r).data(1:nx,1)';
        DD(1).Ey(k,r) = s(r)*mean(ExpC(k).experiment(r).data(1:nx,1)');
        n0 = n0 + nx;
    end
end
DD(1).y = DD(1).y(:,1:find(sum(isnan(DD(1).y))==length(DD(1).t),1,'first')-1);
DD(1).y =log10(DD(1).y );
DD(1).y=cen*ceil(DD(1).y./cen);
DD(1).Ey = DD(1).Ey(:,1:find(sum(isnan(DD(1).Ey))==length(DD(1).t),1,'first')-1);

% Load data and asign to struct
ExpC = load_NGF_dose_response;
DD(2).name = 'NGF dose response';
DD(2).type = 'dose response';
DD(2).measurand = 'relative pErk conc.';
DD(2).t = [30];
DD(2).u = [0,0.001,0.01,0.1,1,10];
DD(2).y = nan(length(DD(2).u),10000);
DD(2).Ey = nan(length(DD(2).u),10);
s = [0.6795,1.1002,0.8158,1.1120,1.2926]; % scaling computed in data_scaling.m
for d = 1:length(DD(2).u)
    n0 = 0;
    for r = 1:length(ExpC(d).experiment)
        nx = size(ExpC(d).experiment(r).data,1);
        DD(2).y(d,n0+1:n0+nx) = s(r)*ExpC(d).experiment(r).data(1:nx,1)';
        DD(2).Ey(d,r) = s(r)*mean(ExpC(d).experiment(r).data(1:nx,1)');
        n0 = n0 + nx;
    end
end
DD(2).y = DD(2).y(:,1:find(sum(isnan(DD(2).y))==length(DD(2).u),1,'first')-1);
DD(2).y =log10(DD(2).y );
DD(2).y=cen*ceil(DD(2).y./cen);
DD(2).Ey = DD(2).Ey(:,1:find(sum(isnan(DD(2).Ey))==length(DD(2).u),1,'first')-1);                
% SPECIFY MODEL AND DATA
M.mixture.type = 'normal';
M.label.x = 'pErk [UI]';
M.label.y = 'probability density';



% EXPERIMENT
i = 1;


D{i}.name = '0 min';
D{i}.description = [];
data=DD(1).y(1,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=0;
N=1;
M.experiment(i).name  = '0 min';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_0_1),1\log(10)*exp(log_sigma_0_2)};

% EXPERIMENT
i = i+1;


D{i}.name = '5 min';
D{i}.description = [];
data=DD(1).y(2,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=5;
N=1;
M.experiment(i).name  = '5 min';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_5_1),1\log(10)*exp(log_sigma_5_2)};

% EXPERIMENT
i = i+1;


D{i}.name = '15 min';
D{i}.description = [];
data=DD(1).y(3,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;
t=15;
N=1;

M.experiment(i).name  = '15 min';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_15_1),1\log(10)*exp(log_sigma_15_2)};


% EXPERIMENT
i = i+1;


D{i}.name = '30 min';
D{i}.description = [];
data=DD(1).y(4,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=1;

M.experiment(i).name  = '30 min';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_1),1\log(10)*exp(log_sigma_30_2)};

% EXPERIMENT
i = i+1;


D{i}.name = '60 min';
D{i}.description = [];
data=DD(1).y(5,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=60;
N=1;
M.experiment(i).name  = '60 min';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_60_1),1\log(10)*exp(log_sigma_60_2)};


%Dose response experiments
% EXPERIMENT
i = i+1;


D{i}.name = '0 nM';
D{i}.description = [];
data=DD(2).y(1,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=0;
M.experiment(i).name  = '0 nM';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_0_1),1\log(10)*exp(log_sigma_30_0_2)};


% EXPERIMENT
i = i+1;

D{i}.name = '0.001 nM';
D{i}.description = [];
data=DD(2).y(2,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=0.001;
M.experiment(i).name  = '0.001 nM';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_0001_1),1\log(10)*exp(log_sigma_30_0001_2)};

% EXPERIMENT
i = i+1;


D{i}.name = '0.01 nM';
D{i}.description = [];
data=DD(2).y(3,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=0.01;
M.experiment(i).name  = '0.01 nM';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_001_1),1\log(10)*exp(log_sigma_30_001_2)};

% EXPERIMENT
i = i+1;

D{i}.name = '0.1 nM';
D{i}.description = [];
data=DD(2).y(4,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=0.1;
M.experiment(i).name  = '0.1 nM';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_01_1),1\log(10)*exp(log_sigma_30_01_2)};

% EXPERIMENT
i = i+1;

D{i}.name = '1 nM';
D{i}.description = [];
data=DD(2).y(5,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=1;
M.experiment(i).name  = '1 nM';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_1_1),1\log(10)*exp(log_sigma_30_1_2)};


% EXPERIMENT
i = i+1;

D{i}.name = '10 nM';
D{i}.description = [];
data=DD(2).y(6,:)';
D{i}.data.uncensored =  data(~isnan(data));
D{i}.data.censored = [];
D{i}.observation_interval = cen;

t=30;
N=10;
M.experiment(i).name  = '10 nM';
M.experiment(i).size  = 2;
M.experiment(i).w     = {1-w, w};
M.experiment(i).mu    = {1\log(10)*log(exp(SE)*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2)))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+((exp(k2)*exp(KD)*N*exp(k3T))/(exp(k2)*exp(KD)*N+exp(k2))))),1\log(10)*log(exp(SE)*((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))*(1-exp(k5)*exp(-t*(((exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))) + exp(k4) + exp(k5)))/(exp(k4) + exp(k5)))+exp(k4))/(exp(k4)+exp(k5)+(exp(k2)*exp(KD)*N*exp(k3T)*exp(B))/(exp(k2)*exp(KD)*N+exp(k2))))};
M.experiment(i).sigma = {1\log(10)*exp(log_sigma_30_10_1),1\log(10)*exp(log_sigma_30_10_2)};




% Compile model
% (This generates the functional expression of parameters and derivatives.)
[M,parameters.constraints] = getMixtureModel(M,parameters.sym);


