% getMixtureDensity computes the densities associated to a model.
%
% USAGE:
% ======
% [Dens] = getMixtureDensity(x,theta,M)
% [Dens] = getMixtureDensity(x,theta,M,Mc)
% [Dens] = getMixtureDensity(x,theta,M,Mc,options)
%
% INPUTS:
% =======
% x ... grid on which density is evaluated
% theta ... parameter value for model
% M ... mixture model for given data
% Mc ... mixture model for censoring. If Mc is empty, it is assumed that no
%   censoring occures.
% options ... options of this algorithm
%   .refinement ... refinement of x for the numerical integration
%   (default = 4).
%
% Outputs:
% ========
% Dens ... object containing all necessary densities
%   .experiment(j) ... densities for experiment j.
%       .fe ... probability density of event
%       .Fe ... cumulative probability of event
%       .fc ... probability density of censoring
%       .Fc ... cumulative probability of censoring
%       .fec ... conditional probability density of event
%       .Fec ... conditional cumulative probability of event
%       .fce ... conditional probability density of censoring
%       .Fce ... conditional cumulative probability of censoring
%
% 2012/07/02 Jan Hasenauer
% modified Eva-Maria Geissen
% function Dens = getMixtureDensity(x,theta,M,Mc,options)
function Dens = getMixtureDensity(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 3
    x = varargin{1};
    theta = varargin{2};
    M = varargin{3};
else
    error('Not enought inputs.')
end

% Check for censoring model
if nargin >= 4
    Mc = varargin{4};
else
    Mc = [];
end

% Check options
options.refinement = 20;
if nargin == 5
    options = setdefault(varargin{5},options);
end

%% INITIALIZATION
theta = theta(:);
r = options.refinement;
% Integration grid
xi  = x(1);
dx = (x(2)-x(1))/r;

for i = 2:length(x)
    xi = [xi;[linspace(x(i-1)+dx,x(i),(r))]'];
end
ind =r*[1:length(x)]' - (r-1);

% Loop: experiments
for j = 1:length(M.experiment)
    
    % Initialization
    f  = zeros(size(xi));
    F  = zeros(size(xi));
    fc  = zeros(size(xi));
    Fc  = zeros(size(xi));
    
    % Loop: mixture elements
    for k = 1:M.experiment(j).size
        switch M.mixture.type
            case 'normal'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                mk = M.experiment(j).mu_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                
                yki = (xi-mk)/sk;
                
                %% DENSITY CALCULATION
                fk = wk*exp(-0.5*yki.^2)./(sqrt(2*pi)*sk);
                Fk = 0.5*wk*erfc(-yki ./sqrt(2));
                
            case 'log-normal'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                mk = M.experiment(j).mu_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                
                yki = (log(xi)-mk)/sk;
                
                %% DENSITY CALCULATION
                fk = wk*exp(-0.5*yki.^2)./(sqrt(2*pi)*sk*xi);
                Fk = 0.5*wk*erfc(-yki ./sqrt(2));
                
            case 'Johnson SU'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                gk = M.experiment(j).gamma_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                lk = M.experiment(j).lambda_fun{k}(theta);
                xk = M.experiment(j).xi_fun{k}(theta);
                
                zki  = (xi - xk)/lk;
                yki  = gk + sk*asinh(zki);
                
                %% DENSITY CALCULATION
                fk = wk*sk/lk*exp(-0.5*yki.^2)./(sqrt(2*pi)*sqrt(zki.^2 + 1));
                Fk = 0.5*wk*erfc(-yki ./sqrt(2));
                
            case 'gamma'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                ak = M.experiment(j).alpha_fun{k}(theta);
                bk = M.experiment(j).beta_fun{k}(theta);
                
                %% DENSITY CALCULATION
                fk = wk*gampdf(xi,ak,1/bk);
                Fk = wk*gamcdf(xi,ak,1/bk);
                
            otherwise
                error('Unknown distribution type. Only ''normal'', ''log-normal'', ''Johnson SU'' and ''gamma'' are allowed!');
        end
        % Assignment
        f = f + fk;
        F = F + Fk;
        Dens.experiment(j).fek{k} = fk;
        Dens.experiment(j).Fek{k} = Fk;
    end
    
    if ~isempty(Mc)
        % Loop: mixture elements
        for k = 1:Mc.experiment(j).size
            switch Mc.mixture.type
                case 'log-normal'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk = Mc.experiment(j).w_fun{k}(theta);
                    mk = Mc.experiment(j).mu_fun{k}(theta);
                    sk = Mc.experiment(j).sigma_fun{k}(theta);
                    
                    yki = (log(xi)-mk)/sk;
                    
                    %% DENSITY CALCULATION
                    fck = wk*exp(-0.5*yki.^2)./(sqrt(2*pi)*sk*xi);
                    Fck = 0.5*wk*erfc(-yki ./sqrt(2));
                    fc = fc + fck;
                    Fc = Fc + Fck;
                    
                case 'Johnson SU'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk = Mc.experiment(j).w_fun{k}(theta);
                    gk = Mc.experiment(j).gamma_fun{k}(theta);
                    sk = Mc.experiment(j).sigma_fun{k}(theta);
                    lk = Mc.experiment(j).lambda_fun{k}(theta);
                    xk = Mc.experiment(j).xi_fun{k}(theta);
                    
                    zki  = (xi - xk)/lk;
                    yki  = gk + sk*asinh(zki);
                    
                    %% DENSITY CALCULATION
                    fck = wk*sk/lk*exp(-0.5*yki.^2)./(sqrt(2*pi)*sqrt(zki.^2 + 1));
                    Fck = 0.5*wk*erfc(-yki ./sqrt(2));
                    fc = fc + fck;
                    Fc = Fc + Fck;
                    
                case 'gamma'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk = Mc.experiment(j).w_fun{k}(theta);
                    ak = Mc.experiment(j).alpha_fun{k}(theta);
                    bk = Mc.experiment(j).beta_fun{k}(theta);
                    
                    %% DENSITY CALCULATION
                    fck = wk*gampdf(xi,ak,1/bk);
                    Fck = wk*gamcdf(xi,ak,1/bk);
                    fc = fc + fck;
                    Fc = Fc + Fck;
                    
                case 'delta'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk = 1;
                    ak = Mc.experiment(j).location;
                    %% DENSITY CALCULATION
                    if ~isempty(ak)
                        fck = wk*dirac(xi-ak);
                        Fck = wk*heaviside(xi-ak);
                        else
                        fck = zeros(length(xi),1);
                        Fck = zeros(length(xi),1);
                    end
                    
                    fc = fc + fck;
                    Fc = Fc + Fck;
                otherwise
                    error('Unknown distribution type. Only ''log-normal'' and ''Johnson SU'' are allowed!');
            end
            % Assignment
            Dens.experiment(j).fck{k} = fck;
            Dens.experiment(j).Fck{k} = Fck;
        end
        
        %% ASSIGNMENT
        Dens.experiment(j).fe = f(ind);
        Dens.experiment(j).fc = fc(ind);
        Dens.experiment(j).Fe = F(ind);
        Dens.experiment(j).Fc = Fc(ind);
        
        %% EVALUATION OF MIXTURE PROBABILITY OF DATA POINTS
        % Evaluation of event distribution
        fA  = [0;f(2:end) .*(1-Fc(2:end))];
        fB  = [0;fc(2:end).*(1-F(2:end) )];
        
        
        % Evaluation of cumulative event probability
        FA = [0;cumsum(0.5*dx*(fA(1:end-1)+fA(2:end)))];
        FB = [0;cumsum(0.5*dx*(fB(1:end-1)+fB(2:end)))];
        
        if isequal(Mc.mixture.type,'delta') && ~isempty(ak)
            [~,index1]=min(abs(xi-ak));
            fB(index1)=1-F(index1);
            FB = [0;cumsum(0.5*(fB(1:end-1)+fB(2:end)))];
            
        end
        % Loop: mixture elements
        for k = 1:M.experiment(j).size
            fAk = [0;Dens.experiment(j).fek{k}(2:end) .*(1-Fc(2:end))];
            FAk = [0;cumsum(0.5*dx*(fAk(1:end-1)+fAk(2:end)))];
            Dens.experiment(j).fekc{k} = fAk(ind);
            Dens.experiment(j).Fekc{k} = FAk(ind);
        end
        
        %% ASSIGNMENT
        Dens.experiment(j).fec = fA(ind);
        Dens.experiment(j).fce = fB(ind);
        Dens.experiment(j).Fec = FA(ind);
        Dens.experiment(j).Fce = FB(ind);
        
    else
         % Loop: mixture elements
        for k = 1:M.experiment(j).size

            Dens.experiment(j).fek{k} = Dens.experiment(j).fek{k}(ind);
            Dens.experiment(j).Fek{k} = Dens.experiment(j).Fek{k}(ind);
            Dens.experiment(j).fekc{k} = Dens.experiment(j).fek{k};
            Dens.experiment(j).Fekc{k} = Dens.experiment(j).Fek{k};
        end
        %% ASSIGNMENT
        Dens.experiment(j).fe  = f(ind);
        Dens.experiment(j).Fe  = F(ind);
        Dens.experiment(j).fc  = zeros(size(fc(ind)));
        Dens.experiment(j).Fc  = zeros(size(Fc(ind)));
        Dens.experiment(j).fec = f(ind);
        Dens.experiment(j).Fec = F(ind);
        Dens.experiment(j).fce = zeros(size(f(ind)));
        Dens.experiment(j).Fce = zeros(size(F(ind)));
        
    end
    
end
