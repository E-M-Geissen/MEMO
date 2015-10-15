% performModelSelection performs a backward model selection on the provided mixture
%    model. The routine tries to eliminate mixture components while
%    achieving a good model-data agreement. As selection criterion the
%    likelihood-ration-test with user-supplied p-value, the Bayesian Information 
%    Criterion (BIC), and the Akaike Information Criterion are employed.
%    As the number of model alternatives grows exponentially instaead of using a
%    brute-force approach, we applied a sequential procedure.
%    First all possible eliminations of a single mixture component are
%    tried. The reductions which fail to provide an acceptable model are
%    disregarded and only the combination of two reductions which were
%    acceptable individual is tried. Following this procedure, for the next
%    iteration all double deletion of mixture components which did not
%    result in an acceptable model are disregarded, ... and so on.
%    As many component deletions can be excluded early on, this approach
%    results in a tremendous reduction of the computational complexity.
%
% USAGE:
% ======
% [...] = performModelSelection(theta,M,Mc,D)
% [...] = performModelSelection(theta,M,Mc,D,options)
% [M_opt,Mc_opt,parameters_opt,S,R] = performModelSelection(...)
%
% INPUTS:
% =======
% theta ... parameter value for model
% M ... mixture model for given data
% D ... measurement data
% options ... options of the algorithm:
%   .criterion ... criterion based on which to decide on goodness of fit 'BIC', 'AIC', 'likelihood ratio' ( default = 'BIC')
%   .p_value ... p-value for model acceptance, only consiedred if .criterion = 'likelihood ratio' (default = 0.05)
%   .print_model ... printing options:
%       = 'none' (default) ... no output on the screen.
%       = 'all' ... the properties of all models which have been analyzed
%           is printed on the screen.
%       = 'accepted' ... the properties of all accepted models
%           is printed on the screen.
%   .optim_opt options for optimization via optimizeMultiStart.m
%        .n_starts number of starts of the local optimization (default = 0)
%        .plot plot the results during the computation (default = 'false').
%        .mode output mode 'normal' or 'silent' (default = 'silent')
%        .fmincon ... options for local optimization.
%
% Outputs:
% ========
% M_opt ... optimal model according to the selcted model criterion.
% Mc_opt ... optimal censoring model, given M_opt, according to the selcted model criterion.
% parameters_opt ... parameters of the optimal model.
% S ... set of all analyzed models and the corresponding parameters, sorted according to specified criterion.
% R ... component deletion considered.
%
% 2012/06/15 Jan Hasenauer
% modified Eva-Maria Geissen
%
% function [M_opt,parameters_opt,S,R] = performModelSelection(parameters,M,D)
function [M_opt,Mc_opt,parameters_opt,S,R] = performModelSelection(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 3
    parameters = varargin{1};
    M = varargin{2};
    Mc = varargin{3};
    D = varargin{4};
else
    error('performModelSelection.m requires at least three inputs.')
end

% Set defaults and assign user-supplied values
options.p_value = 0.05;
options.criterion = 'BIC'; % 'likelihood ratio', 'AIC', 'BIC'
options.print_model = 'none'; %'all', 'none', 'accepted'
options.optim_opt.n_starts = 0;
options.optim_opt.plot = 'false';
options.optim_opt.mode = 'silent';
options.optim_opt.fmincon = optimset('algorithm','active-set',...
    'display','off',...
    'GradObj','on',...
    'MaxIter',1000,...
    'MaxFunEvals',1000*parameters.number);

if nargin == 5
    options = setdefault(varargin{5},options);
end

%% DETERMINE PARAMETERS OF CENSORING MODEL
par_cen = [];
if ~isempty(Mc)
    par_cen = getSymbolicParameterVector(Mc);
end

%% EVALUATION OF THE NUMBER OF DATA POINTS
N = 0;
for j = 1:length(D)
    N = N + length(D{j}.data.censored) + length(D{j}.data.uncensored);
end

%% MODEL REDUCTION
% COnstruct set of possible modl reductions
R = {};
for j = 1:length(M.experiment)
    % Construct new set
    Rj = {};
    if M.experiment(j).size == 2
        Rj{1} = [j,1];
        Rj{2} = [j,2];
    end
        if M.experiment(j).size == 3
        Rj{1}  = [j,1];
        Rj{2}  = [j,2];
        Rj{3}  = [j,3];
        Rj{4}  = [j,1;j,2];
        Rj{5}  = [j,1;j,3];
        Rj{6}  = [j,2;j,3];
    end
    if M.experiment(j).size >= 4
        error('This is not implemented yet.');
    end
    % Assign perturbations
    if isempty(R)
        R = Rj;
    else
        for k = 1:length(Rj)
            R{end+1} = Rj{k};
        end
        n = length(R)-length(Rj);
        for k1 = 1:length(Rj)
            for k2 = 1:n
                R{end+1} = [R{k2};Rj{k1}];
            end
        end
    end
end
% Largest number of reductions
Smax = 0;
for j = 1:length(M.experiment)
    Smax = Smax + M.experiment(j).size - 1;
end
% Sort set
Rs = {};
for k1 = 1:Smax
    for k2 = 1:length(M.experiment)
        for j = 1:length(R)
            if (size(R{j},1) == k1) && (min(R{j}(:,1)) == k2)
                Rs{end+1} = sortrows(R{j});
            end
        end
    end
end
R = Rs;

%% PREPARATION OF MODEL
% Calculation of fit quality
parameters.AIC  = -2*parameters.MS.MAP.logPost + 2*parameters.number;
parameters.BIC  = -2*parameters.MS.MAP.logPost +   parameters.number*log(N);
% Constucton of component index vector
for j = 1:length(M.experiment)
    if ~isfield(M.experiment(j),'component_index')
        M.experiment(j).component_index = 1:M.experiment(j).size;
    else
        if isempty(M.experiment(j).component_index)
            M.experiment(j).component_index = 1:M.experiment(j).size;
        end
    end
end

%% INITIAL MODEL SET
S{1}.model = M;
S{1}.parameters = parameters;

%% MODEL REDUCTION
disp(['Progress']);
% Active reductions
ind = 1:length(R);
% Loop: As long as model reductions are possible
for i = 1:length(R)
    % Check whether analysis is promising
    if max(ind == i)
        
        %% Construction of reduced model
        % Reduction
        m = R{i}(:,1);
        n = R{i}(:,2);
        % Initialize reduced model
        rM = M;
        orM = M;
        % Interative reduction of model
        for j = 1:length(m)
            % Component index
            c = find(orM.experiment(m(j)).component_index == n(j));
            % Reduced component number
            rM.experiment(m(j)).size = orM.experiment(m(j)).size - 1;
            % Reinitialization of component properties
            switch rM.mixture.type
                case 'normal'
                    rM.experiment(m(j)).w      = {};
                    rM.experiment(m(j)).mu     = {};
                    rM.experiment(m(j)).sigma  = {};
                    
                case 'log-normal'
                    rM.experiment(m(j)).w      = {};
                    rM.experiment(m(j)).mu     = {};
                    rM.experiment(m(j)).sigma  = {};
                case 'Johnson SU'
                    rM.experiment(m(j)).w      = {};
                    rM.experiment(m(j)).gamma  = {};
                    rM.experiment(m(j)).sigma  = {};
                    rM.experiment(m(j)).lambda = {};
                    rM.experiment(m(j)).xi     = {};
                    
                case 'gamma'
                    rM.experiment(m(j)).w     = {};
                    rM.experiment(m(j)).alpha = {};
                    rM.experiment(m(j)).beta  = {};
            end
            % Reduction of model
            wc = orM.experiment(m(j)).w{c};
            for l = 1:rM.experiment(m(j)).size
                switch rM.mixture.type
                    case 'normal'
                        rM.experiment(m(j)).w{l}      = orM.experiment(m(j)).w{l + (l >= c)} + wc*(l==1);
                        rM.experiment(m(j)).mu{l}     = orM.experiment(m(j)).mu{l + (l >= c)};
                        rM.experiment(m(j)).sigma{l}  = orM.experiment(m(j)).sigma{l + (l >= c)};
                    case 'log-normal'
                        rM.experiment(m(j)).w{l}      = orM.experiment(m(j)).w{l + (l >= c)} + wc*(l==1);
                        rM.experiment(m(j)).mu{l}     = orM.experiment(m(j)).mu{l + (l >= c)};
                        rM.experiment(m(j)).sigma{l}  = orM.experiment(m(j)).sigma{l + (l >= c)};
                    case 'Johnson SU'
                        rM.experiment(m(j)).w{l}      = orM.experiment(m(j)).w{l + (l >= c)} + wc*(l==1);
                        rM.experiment(m(j)).gamma{l}  = orM.experiment(m(j)).gamma{l + (l >= c)};
                        rM.experiment(m(j)).sigma{l}  = orM.experiment(m(j)).sigma{l + (l >= c)};
                        rM.experiment(m(j)).lambda{l} = orM.experiment(m(j)).lambda{l + (l >= c)};
                        rM.experiment(m(j)).xi{l}     = M.experiment(m(j)).xi{l + (l >= c)};
                    case 'gamma'
                        rM.experiment(m(j)).w{l}      = orM.experiment(m(j)).w{l + (l >= c)} + wc*(l==1);
                        rM.experiment(m(j)).alpha{l}  = orM.experiment(m(j)).mu{l + (l >= c)};
                        rM.experiment(m(j)).beta{l}   = orM.experiment(m(j)).sigma{l + (l >= c)};
                        
                end
            end
            % Assign component index to map back to original model
            rM.experiment(m(j)).component_index = setdiff(rM.experiment(m(j)).component_index,rM.experiment(m(j)).component_index(c));
            % Assignment of current model to old model for next reduction step
            orM = rM;
        end
        
        %% Construction of parameter vector
        % Determine parameter vector of reduced model
        rparameters.sym = getSymbolicParameterVector(rM);
        rparameters.sym = [rparameters.sym;par_cen];
        % Determine index mapping
        index_mapping = [];
        for l = 1:length(rparameters.sym)
            index_mapping(l) = find(parameters.sym == rparameters.sym(l));
        end
        index_mapping = sort(index_mapping);
        % Construct reduced parameter object
        rparameters.sym  = parameters.sym(index_mapping);
        rparameters.name = parameters.name(index_mapping);
        rparameters.number = length(index_mapping);
        rparameters.min = parameters.min(index_mapping);
        rparameters.max = parameters.max(index_mapping);
        rparameters.guess = parameters.MS.MAP.par(index_mapping);
        
        %% Compiling of model, optimization, and fit analysis
        % Compile models
        rM  = getMixtureModel(rM ,rparameters.sym);
        % Assigment of objective function
        if isempty(Mc)
            logPosterior = @(theta,opt) logLikelihoodMM(theta,rM,D,opt);
        else
            rMc = getMixtureModel( Mc,rparameters.sym);
            logPosterior = @(theta,opt) logLikelihoodMMc(theta,rM,rMc,D,opt);
        end
        % Optimization
        rparameters = optimizeMultiStart(rparameters,logPosterior,options.optim_opt);
        % Performance of reduced model
        rparameters.LR.Lambda = exp(rparameters.MS.MAP.logPost - parameters.MS.MAP.logPost);
        rparameters.LR.p = 1 - chi2cdf(-2*(rparameters.MS.MAP.logPost - parameters.MS.MAP.logPost),parameters.number - rparameters.number);
        rparameters.AIC  = -2*rparameters.MS.MAP.logPost + 2*rparameters.number;
        rparameters.BIC  = -2*rparameters.MS.MAP.logPost + rparameters.number*log(N);
        
        %% Store of model and updating of model indexes
        S{i+1}.model  = rM;
        if ~isempty(Mc)
            S{i+1}.modelc = rMc;
        end
        S{i+1}.parameters = rparameters;
        
        %% CHECK FOR ACCEPTATBILITY OF MODEL
        acc = 'false';
        switch options.criterion
            % Likelihood ratio test
            case 'likelihood ratio'
                if rparameters.LR.p <= options.p_value
                    acc = 'true';
                end
                % Test based on Akaike information criterion
            case 'AIC'
                if rparameters.AIC <= parameters.AIC
                    acc = 'true';
                end
                % Test based on Bayesian information criterion
            case 'BIC'
                if rparameters.BIC <= parameters.BIC
                    acc = 'true';
                end
        end
        
        %% Update index of considered models
        if strcmp(acc,'false')
            % Loop: Remaining model reductions
            J = [];
            for j = ind(ind > i)
                % Determine whether rejected reduction is elemet of other
                % reduction, which can than be eliminated
                d = 1;
                for k = 1:size(R{i},1)
                    d = d*max((R{j}(:,1)==R{i}(k,1)).*(R{j}(:,2)==R{i}(k,2)));
                end
                % Eliminate model is possible
                if d == 1
                    J = [J,j];
                end
            end
            ind = setdiff(ind,J);
            
        end
        
        % Progess
        str = [num2str(round(1000*(length(R)-sum(ind>i))/length(R))/10,'%.1f %%') ' - Elimination: '];
        for k = 1:size(R{i},1)
            str = [str 'e = ' num2str(R{i}(k,1),'%d') ', c = ' num2str(R{i}(k,2),'%d')];
            if k ~= size(R{i},1)
                str = [str ' / '];
            end
        end
        if strcmp(acc,'true')
            str = [str ' -> accepted'];
        else
            str = [str ' -> rejected'];
        end
        
        disp(str);
        
    else
        S{i+1}.model = [];
        S{i+1}.modelc = [];
        S{i+1}.parameters = [];
    end
    
end

%% PRINT MODELS ON SCREEN
switch options.print_model
    case 'all'
        % Print model on screen
        printModel(S{1}.model,S{1}.parameters);
        for i = 1:length(ind)
            printModel(S{ind(i)+1}.model,S{ind(i)+1}.parameters);
        end
    case 'accepted'
        % Print model on screen
        printModel(S{1}.model,S{1}.parameters);
        for i = 1:length(ind)
            if S{ind(i)+1}.parameters.LR.p < options.p_value
                printModel(S{ind(i)+1}.model,S{ind(i)+1}.parameters);
            end
        end
end

%% DETERMINE MOST PLAUSIBLE MODEL
% Initialization
i_opt  = ind(1);
M_opt  = S{i_opt+1}.model;
if ~isempty(Mc)
    Mc_opt = S{i_opt+1}.modelc;
else
    Mc_opt = [];
end
parameters_opt = S{i_opt+1}.parameters;
% Loop over models
for i = 1:length(ind)
    % Better?
    better = 'false';
    switch options.criterion
        % Likelihood ratio test
        case 'likelihood ratio'
            if S{ind(i)+1}.parameters.LR.p <= S{i_opt+1}.parameters.LR.p
                better = 'true';
            end
            % Test based on Akaike information criterion
        case 'AIC'
            if S{ind(i)+1}.parameters.AIC <= S{i_opt+1}.parameters.AIC
                better = 'true';
            end
            % Test based on Bayesian information criterion
        case 'BIC'
            if S{ind(i)+1}.parameters.BIC <= S{i_opt+1}.parameters.BIC
                better = 'true';
            end
    end
    % Criterion
    if strcmp(better,'true')
        i_opt = ind(i);
        M_opt  = S{ind(i)+1}.model;
        if ~isempty(Mc)
            Mc_opt = S{ind(i)+1}.modelc;
        end
        parameters_opt = S{ind(i)+1}.parameters;
    end
end


