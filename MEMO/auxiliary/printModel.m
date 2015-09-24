% printModel prints the model properties on the screen in a human readable
%   form. If the parameters are provided, the expressions of the model
%   properties are evaluated and the numerical values are shown.
%
% USAGE:
% ======
% printModel(M)
% printModel(M,parameters)
%
% INPUTS:
% =======
% M ... mixture model containing symbolic experssions
%    .mixture.type ... type of mixture
%    .experiment{j} ... j-th experiment
%       .name ... name of experiment
%       .size ... number of mixture components
%       For mixture of log-normals:
%           f(x) = sum_k w_k*Phi((x-mu_k)/sigma_k)
%          .w{k}
%          .mu{k}
%          .sigma{k}
%             ... symbolic expression for mixture parameters w_k, mu_k and
%                 sigma_k as function of the symbolic elements of theta.
%          .w_fun{k}
%          .mu_fun{k}
%          .sigma_fun{k}
%             ... mixture parameters w_k, mu_k and sigma_k as function
%                 of theta (e.g., w_k = w_fun{k}(theta).
%          .dwdtheta_fun{k}
%          .dmudtheta_fun{k}
%          .dsigmadtheta_fun{k}
%              ... derivatives of mixture parameters w_k, mu_k and sigma_k
%                  with respect to theta as function of theta.
%       For Johnson SU distribution:
%           f(x) = sum_k w_k*Phi(gamma_k + sigma_k*arcsinh((x-xi_k)/lambda_k))
%          .w{k}
%          .gamma{k}
%          .sigma{k}
%          .sigma{k}
%          .lambda{k}
%          .xi{k}
%             ... symbolic expression for mixture parameters w_k, gamma_k,
%                 sigma_k, lamda_k, and xi_k as function of symbolic
%                 elements of theta.
%          .w_fun{k}
%          .gamma_fun{k}
%          .sigma_fun{k}
%          .sigma_fun{k}
%          .lambda_fun{k}
%          .xi_fun{k}
%             ... mixture parameters w_k, gamma_k, sigma_k, lamda_k, and xi_k
%                 as function of theta (e.g., w_k = w_fun{k}(theta).
%          .w_fun{k}
%          .dgammadtheta_fun{k}
%          .dsigmadtheta_fun{k}
%          .dsigmadtheta_fun{k}
%          .dlambdadtheta_fun{k}
%          .dxidtheta_fun{k}
%              ... derivatives of mixture parameters w_k, gamma_k, sigma_k,
%                  lamda_k, and xi_k with respect to theta as function of theta.
%                  of theta.
% parameters ... struct containg information about parameters, at least:
%   .MS.MAP.par ... maximum a posterior estimate
%   .MS.MAP.logPost ... log-posterior at the maximum a posterior estimate
%
% 2012/05/17 Jan Hasenauer

function printModel(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 1
    M = varargin{1};
else
    error('printModel a model as input.');
end

if nargin == 2
    parameters = varargin{2};
    est_par = 1;
else
    est_par = 0;
end

%% INITIALIZE
space = ' ';
line  = '-';

%% CONSTRUCT COMPONENT INDEX
for j = 1:length(M.experiment)
    if isfield(M.experiment(j),'component_index')
        if isempty(M.experiment(j).component_index)
            M.experiment(j).component_index = 1:M.experiment(j).size;
        end
    else
        M.experiment(j).component_index = 1:M.experiment(j).size;
    end
end

%% CONSTRUCTION OF HEADER
kmax = 0;
% Loop: Experiments
for j = 1:length(M.experiment)
    kmax = max([kmax,M.experiment(j).size,max(M.experiment(j).component_index)]);
end
% Loop: Mixture components
for k = 1:kmax
    strh{k} = [' component ' num2str(k,'%d') ' '];
end

strh{kmax+1}=['input dependency'];


%% PARMATER LABEL
% Loop: Mixture components
for k = 1:kmax
    switch M.mixture.type
        case 'normal'
            strtn{k}{1} = [' w  '];
            strtn{k}{2} = [' mu '];
            strtn{k}{3} = [' sigma '];
        case 'log-normal'
            strtn{k}{1} = [' w  '];
            strtn{k}{2} = [' mu '];
            strtn{k}{3} = [' sigma '];
        case 'Johnson SU'
            strtn{k}{1} = [' w  '];
            strtn{k}{2} = [' gamma '];
            strtn{k}{3} = [' sigma '];
            strtn{k}{4} = [' lambda '];
            strtn{k}{5} = [' xi '];
        case 'gamma'
            strtn{k}{1} = [' w  '];
            strtn{k}{2} = [' alpha '];
            strtn{k}{3} = [' beta  '];
        otherwise
            error('Unknown type of distribution.')
    end
    % Keep track of string length
    for l = 1:length(strtn{k})
        ll{k}{l} = length(strtn{k}{l});
    end
end

if isfield(M.experiment(1),'cond')
    
    strtn{k+1}{1}=['u'];
    strtn{k+1}{2}=['est. error'];
    strtn{k+1}{3}=['sigma error'];
    %    % Keep track of string length
    %     for l = 1:length(strtn{k+1})
    %         ll{k+1}{l} = length(strtn{k+1}{l});
    %     end
    if isfield (M.experiment(1).cond, 'param')
        fun_param={};
        for i=1:length(M.experiment)
            fun_param=[fun_param,M.experiment(i).cond.param ];
        end
        fun_param=unique(fun_param);
        
        for p=1:length(fun_param)
            strtn{k+1}{3+p}=[' ' fun_param{p} ' '];
        end
    else
        error('Please add ''M.experiment(i).cond.param'' to the input parametrization.')
    end
    % Keep track of string length
    for l = 1:length(strtn{k+1})
        ll{k+1}{l} = length(strtn{k+1}{l});
    end
    
end
%% CONSTRUCTION OF DATA
ld = 0;
% Loop: Experiments
for j = 1:length(M.experiment)
    % Name of dataset
    strd{j} = [' ' M.experiment(j).name ' '];
    % Keep track of string length
    ld = max([ld,length(strd{j})]);
end



if isfield(M.experiment(1),'cond')
    %% CONSTRUCTION OF TABLE
    % Loop: Experiments
    for j = 1:length(M.experiment)
        % Mixture component index
        ci = M.experiment(j).component_index;
        % Loop: Mixture components
        for k = 1:kmax
            K = find(k == ci);
            if ~isempty(K)
                K = sum(ci <= ci(K));
                switch M.mixture.type
                    case 'normal'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})     ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).mu{K})    ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).sigma{K}) ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).mu_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).sigma_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                        end
                    case 'log-normal'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})     ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).mu{K})    ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).sigma{K}) ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).mu_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).sigma_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                        end
                    case 'Johnson SU'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})      ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).gamma{K})  ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).sigma{K})  ' '];
                            strt{j}{k}{4} = [' ' char(M.experiment(j).lambda{K}) ' '];
                            strt{j}{k}{5} = [' ' char(M.experiment(j).xi{K})     ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')      ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).gamma_fun{K}(parameters.MS.MAP.par),'%.2f')  ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).sigma_fun{K}(parameters.MS.MAP.par),'%.2f')  ' '];
                            strt{j}{k}{4} = [' ' num2str(M.experiment(j).lambda_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                            strt{j}{k}{5} = [' ' num2str(M.experiment(j).xi_fun{K}(parameters.MS.MAP.par),'%.2f')     ' '];
                        end
                    case 'gamma'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})     ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).alpha{K})    ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).beta{K}) ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).alpha_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).beta_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                        end
                end
                if (M.experiment(j).size == 1) && (M.experiment(j).component_index == K)
                    if nargin == 1
                        strt{j}{1}{1} = [' 1 '];
                    else
                        strt{j}{1}{1} = [' 1.00 '];
                    end
                end
                
                % For this column no entry, as not contained in 'compenent_index'
            else
                switch M.mixture.type
                    case 'normal'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                    case 'log-normal'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                    case 'Johnson SU'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                        strt{j}{k}{4} = [''];
                        strt{j}{k}{5} = [''];
                    case 'gamma'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                end
                
            end
            % Keep track of string length
            for l = 1:length(strt{j}{k})
                ll{k}{l} = max([ll{k}{l},length(strt{j}{k}{l})]);
            end
        end
        % Process input dependency
        if est_par == 0 || length(M.experiment(j).cond.sigma)==2
            strt{j}{kmax+1}{1} = [' '];
            strt{j}{kmax+1}{2} = [' '];
            strt{j}{kmax+1}{3} = [' '];
            for p=1:length(fun_param)
                strt{j}{kmax+1}{3+p}=[' '];
            end
            
        else
            strt{j}{kmax+1}{1} = [' ' num2str(M.experiment(j).cond.y{1},'%.2f')    ' '];
            strt{j}{kmax+1}{3} = [' ' num2str(M.experiment(j).cond.sigma{1},'%.2f')    ' '];
            strt{j}{kmax+1}{2} = [' ' num2str(M.experiment(j).cond.e_fun{1}(parameters.MS.MAP.par),'%.2f') ' '];
            
            %         % find coresponding parameters for each experiment
            %         ind_save=[];
            %         for p=1:length(strtn{kmax+1})-3
            %             ind=find(ismember(M.experiment(j).cond.param,strtrim(strtn{kmax+1}(3+p))));
            %             if~isempty(ind)
            %                 strt{j}{kmax+1}{3+p}=[' ' num2str(parameters.MS.MAP.par(find(parameters.sym ==char(M.experiment(j).cond.param{ind})))  ,'%.2f') ' '];
            %                 ind_save=[ind_save,ind];
            %             end
            %         end
            %         for p=1:length(strtn{kmax+1})-3
            %             if isempty(find(ind_save==p))
            %                 strt{j}{kmax+1}{3+p}=' ';
            %             end
            %         end
            
            % find coresponding parameters for each experiment
            ind_save=[];
            for p=1:length(strtn{kmax+1})-3
                ind=find(ismember(strtrim(strtn{kmax+1}(3+p)),M.experiment(j).cond.param));
                
                if~isempty(ind)
                    ind_map=find(ismember(M.experiment(j).cond.param,strtrim(strtn{kmax+1}(3+p))));
                    strt{j}{kmax+1}{3+p}=[' ' num2str(parameters.MS.MAP.par(find(parameters.sym ==char(M.experiment(j).cond.param{ind_map})))  ,'%.2f') ' '];
                    ind_save=[ind_save,p];
                end
            end
            for p=1:length(strtn{kmax+1})-3
                if isempty(find(ind_save==p))
                    strt{j}{kmax+1}{3+p}=' ';
                end
            end
            
        end
        
        % Keep track of string length
        for l = 1:length(strt{j}{kmax+1})
            ll{kmax+1}{l} = max([ll{kmax+1}{l},length(strt{j}{kmax+1}{l})]);
        end
    end
    
    %% COMPUTATION OF LENGTH OF DIFFERENT COMPONENTS
    % Loop: Mixture components
    for k = 1:kmax+1
        for l = 1:length(strt{j}{k})
            if l == 1
                lk{k} = ll{k}{l};
            else
                lk{k} = lk{k} + ll{k}{l} + 1;
            end
        end
    end
    
    %% ADAPTATION OF STRING LENGTH
    % 1st column
    
    strds = [' Dataset ' space(ones(1,ld-9))];
    ld=max(ld,length(strds));
    % Loop: Experiments
    for j = 1:length(M.experiment)
        strd{j} = [strd{j} space(ones(1,ld-length(strd{j})))];
    end
    
    % Other columns
    % Loop: Mixture components
    for k = 1:kmax+1
        % Head line
        strh{k} = [space(ones(1,floor((lk{k}-length(strh{k}))/2))) strh{k} ...
            space(ones(1, ceil((lk{k}-length(strh{k}))/2)))];
        strl{k} = [line(ones(1,length(strh{k})))];
        % Names of parameters
        % Loop: Mixture component parameter
        for l = 1:length(strt{j}{k})
            % Name of parameter
            strtn{k}{l} = [space(ones(1,floor((ll{k}{l}-length(strtn{k}{l}))/2))) strtn{k}{l} ...
                space(ones(1, ceil((ll{k}{l}-length(strtn{k}{l}))/2)))];
        end
        % Data
        % Loop: Experiments
        for j = 1:length(M.experiment)
            % Loop: Mixture component parameter
            for l = 1:length(strt{j}{k})
                strt{j}{k}{l} = [space(ones(1,ll{k}{l}-length(strt{j}{k}{l}))) strt{j}{k}{l}];
            end
        end
    end
    
    %% CONSTRUCT STRING TO PRINT
    % Type of mixture
    strp = ['\n\nType of mixture: ' M.mixture.type '\n \n'];
    
    % 1st row
    strp = [strp ' ' line(ones(1,length(strds))) '-'];
    for k = 1:kmax+1
        if k < kmax+1
            strp = [strp strl{k} '-'];
        else
            strp = [strp strl{k} ' '];
        end
    end
    strp = [strp '\n'];
    
    % 2nd row
    strp = [strp '|' strds '|'];
    % Loop: Mixture components
    for k = 1:kmax+1
        strp = [strp strh{k} '|'];
    end
    strp = [strp '\n'];
    
    % 3rd row
    strp = [strp '|' line(ones(1,length(strds))) '|'];
    for k = 1:kmax+1
        strp = [strp strl{k} '|'];
    end
    strp = [strp '\n'];
    
    % 4th row
    strp = [strp '|' space(ones(1,length(strds))) '|'];
    % Loop: Mixture components
    for k = 1:kmax+1
        for l = 1:length(strtn{k})
            strp = [strp strtn{k}{l} '|'];
        end
    end
    strp = [strp '\n'];
    
    % 5th row
    strp = [strp '|' line(ones(1,length(strds))) '|'];
    for k = 1:kmax+1
        strp = [strp strl{k} '|'];
    end
    strp = [strp '\n'];
    
    % Other rows
    % Loop: Experiments
    for j = 1:length(M.experiment)
        % Name of dataset
        strp = [strp '|' strd{j} '|'];
        % Loop: Mixture component parameter
        for k = 1:kmax+1
            for l = 1:length(strt{j}{k})
                strp = [strp strt{j}{k}{l} '|'];
            end
        end
        % Linebreak
        strp = [strp '\n'];
    end
    
    % Last line
    strp = [strp ' ' line(ones(1,length(strds))) '-'];
    for k = 1:kmax+1
        if k < kmax+1
            strp = [strp strl{k} '-'];
        else
            strp = [strp strl{k} ' '];
        end
    end
else
    %% CONSTRUCTION OF TABLE
    % Loop: Experiments
    for j = 1:length(M.experiment)
        % Mixture component index
        ci = M.experiment(j).component_index;
        % Loop: Mixture components
        for k = 1:kmax
            K = find(k == ci);
            if ~isempty(K)
                K = sum(ci <= ci(K));
                switch M.mixture.type
                    case 'normal'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})     ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).mu{K})    ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).sigma{K}) ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).mu_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).sigma_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                        end
                    case 'log-normal'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})     ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).mu{K})    ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).sigma{K}) ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).mu_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).sigma_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                        end
                    case 'Johnson SU'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})      ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).gamma{K})  ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).sigma{K})  ' '];
                            strt{j}{k}{4} = [' ' char(M.experiment(j).lambda{K}) ' '];
                            strt{j}{k}{5} = [' ' char(M.experiment(j).xi{K})     ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')      ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).gamma_fun{K}(parameters.MS.MAP.par),'%.2f')  ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).sigma_fun{K}(parameters.MS.MAP.par),'%.2f')  ' '];
                            strt{j}{k}{4} = [' ' num2str(M.experiment(j).lambda_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                            strt{j}{k}{5} = [' ' num2str(M.experiment(j).xi_fun{K}(parameters.MS.MAP.par),'%.2f')     ' '];
                        end
                    case 'gamma'
                        if est_par == 0
                            strt{j}{k}{1} = [' ' char(M.experiment(j).w{K})     ' '];
                            strt{j}{k}{2} = [' ' char(M.experiment(j).alpha{K})    ' '];
                            strt{j}{k}{3} = [' ' char(M.experiment(j).beta{K}) ' '];
                        else
                            strt{j}{k}{1} = [' ' num2str(M.experiment(j).w_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{2} = [' ' num2str(M.experiment(j).alpha_fun{K}(parameters.MS.MAP.par),'%.2f')    ' '];
                            strt{j}{k}{3} = [' ' num2str(M.experiment(j).beta_fun{K}(parameters.MS.MAP.par),'%.2f') ' '];
                        end
                end
                if (M.experiment(j).size == 1) && (M.experiment(j).component_index == K)
                    if nargin == 1
                        strt{j}{1}{1} = [' 1 '];
                    else
                        strt{j}{1}{1} = [' 1.00 '];
                    end
                end
                
                % For this column no entry, as not contained in 'compenent_index'
            else
                switch M.mixture.type
                    case 'normal'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                    case 'log-normal'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                    case 'Johnson SU'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                        strt{j}{k}{4} = [''];
                        strt{j}{k}{5} = [''];
                    case 'gamma'
                        strt{j}{k}{1} = [''];
                        strt{j}{k}{2} = [''];
                        strt{j}{k}{3} = [''];
                end
                
            end
            % Keep track of string length
            for l = 1:length(strt{j}{k})
                ll{k}{l} = max([ll{k}{l},length(strt{j}{k}{l})]);
            end
        end
    end
    
    %% COMPUTATION OF LENGTH OF DIFFERENT COMPONENTS
    % Loop: Mixture components
    for k = 1:kmax
        for l = 1:length(strt{j}{k})
            if l == 1
                lk{k} = ll{k}{l};
            else
                lk{k} = lk{k} + ll{k}{l} + 1;
            end
        end
    end
    
    %% ADAPTATION OF STRING LENGTH
    % 1st column
    
    strds = [' Dataset ' space(ones(1,ld-9))];
    ld=max(ld,length(strds));
    % Loop: Experiments
    for j = 1:length(M.experiment)
        strd{j} = [strd{j} space(ones(1,ld-length(strd{j})))];
    end
    
    % Other columns
    % Loop: Mixture components
    for k = 1:kmax
        % Head line
        strh{k} = [space(ones(1,floor((lk{k}-length(strh{k}))/2))) strh{k} ...
            space(ones(1, ceil((lk{k}-length(strh{k}))/2)))];
        strl{k} = [line(ones(1,length(strh{k})))];
        % Names of parameters
        % Loop: Mixture component parameter
        for l = 1:length(strt{j}{k})
            % Name of parameter
            strtn{k}{l} = [space(ones(1,floor((ll{k}{l}-length(strtn{k}{l}))/2))) strtn{k}{l} ...
                space(ones(1, ceil((ll{k}{l}-length(strtn{k}{l}))/2)))];
        end
        % Data
        % Loop: Experiments
        for j = 1:length(M.experiment)
            % Loop: Mixture component parameter
            for l = 1:length(strt{j}{k})
                strt{j}{k}{l} = [space(ones(1,ll{k}{l}-length(strt{j}{k}{l}))) strt{j}{k}{l}];
            end
        end
    end
    
    %% CONSTRUCT STRING TO PRINT
    % Type of mixture
    strp = ['\n\nType of mixture: ' M.mixture.type '\n \n'];
    
    % 1st row
    strp = [strp ' ' line(ones(1,length(strds))) '-'];
    for k = 1:kmax
        if k < kmax
            strp = [strp strl{k} '-'];
        else
            strp = [strp strl{k} ' '];
        end
    end
    strp = [strp '\n'];
    
    % 2nd row
    strp = [strp '|' strds '|'];
    % Loop: Mixture components
    for k = 1:kmax
        strp = [strp strh{k} '|'];
    end
    strp = [strp '\n'];
    
    % 3rd row
    strp = [strp '|' line(ones(1,length(strds))) '|'];
    for k = 1:kmax
        strp = [strp strl{k} '|'];
    end
    strp = [strp '\n'];
    
    % 4th row
    strp = [strp '|' space(ones(1,length(strds))) '|'];
    % Loop: Mixture components
    for k = 1:kmax
        for l = 1:length(strtn{k})
            strp = [strp strtn{k}{l} '|'];
        end
    end
    strp = [strp '\n'];
    
    % 5th row
    strp = [strp '|' line(ones(1,length(strds))) '|'];
    for k = 1:kmax
        strp = [strp strl{k} '|'];
    end
    strp = [strp '\n'];
    
    % Other rows
    % Loop: Experiments
    for j = 1:length(M.experiment)
        % Name of dataset
        strp = [strp '|' strd{j} '|'];
        % Loop: Mixture component parameter
        for k = 1:kmax
            for l = 1:length(strt{j}{k})
                strp = [strp strt{j}{k}{l} '|'];
            end
        end
        % Linebreak
        strp = [strp '\n'];
    end
    
    % Last line
    strp = [strp ' ' line(ones(1,length(strds))) '-'];
    for k = 1:kmax
        if k < kmax
            strp = [strp strl{k} '-'];
        else
            strp = [strp strl{k} ' '];
        end
    end
    
end
% Likelihood
if est_par == 1
    strp = [strp '\n\nNumber of parameters: ' num2str(parameters.number,'%d') '\n'];
    strp = [strp 'log-posterior: ' num2str(parameters.MS.MAP.logPost,'%10.4e') '\n'];
    if isfield(parameters,'AIC')
        strp = [strp 'AIC: ' num2str(parameters.AIC,'%10.4e') '\n'];
    end
    if isfield(parameters,'BIC')
        strp = [strp 'BIC: ' num2str(parameters.BIC,'%10.4e') '\n'];
    end
    if isfield(parameters,'LR')
        strp = [strp 'Likelihood ration: Lambda = ' num2str(parameters.LR.Lambda,'%10.4e') ...
            ', p-value = ' num2str(parameters.LR.p,'%10.4e') '\n'];
    end
end

strp = [strp '\n\n'];

%% MODIFY STRING TO AVOID ERRORS
strp = regexprep(strp,'%','%%');


%% PRINT
fprintf(strp)

