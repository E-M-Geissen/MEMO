% logLikelihoodMM evaluates the log-likelihood of a mixture model M
%   with respect to a certain dataset D. The mixture model is parametrized
%   using the parameters theta. This paramterization is rather flexible
%   and allows for complex dependencies of the mixture components on theta.
%   It furthermore calculates the gradient of the log_likelihood with
%   respect to the parameters theta.
%
% USAGE:
% ======
% [output] = logLikelihoodMM(theta,M,D)
% [output] = logLikelihoodMM(theta,M,D,options)
% [logL]   = logLikelihoodMM(...)
% [logL,grad] = logLikelihoodMM(...)
%
% INPUTS:
% =======
% theta ... column vector containing the model parameters
% M ... mixture model for given data
%    .mixture.type ... type of mixture
%    .experiment{j} ... j-th experiment
%       .name ... name of experiment
%       .size ... number of mixture components
%       For mixture of normals:
%           f(x) = sum_k w_k*Phi((log(x)-mu_k)/sigma_k)
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
%       For mixture of log-normals:
%           f(x) = sum_k w_k*Phi((x-mu_k)/sigma_k)
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
%       For mixture of Johnson SU distributions:
%           f(x) = sum_k w_k*Phi(gamma_k + sigma_k*arcsinh((x-xi_k)/lambda_k))
%          .w_fun{k}
%          .gamma_fun{k}
%          .sigma_fun{k}
%          .lambda_fun{k}
%          .xi_fun{k}
%             ... mixture parameters w_k, gamma_k, sigma_k, lamda_k, and xi_k
%                 as function of theta (e.g., w_k = w_fun{k}(theta).
%          .dwdtheta_fun{k}
%          .dgammadtheta_fun{k}
%          .dsigmadtheta_fun{k}
%          .dlambdadtheta_fun{k}
%          .dxidtheta_fun{k}
%              ... derivatives of mixture parameters w_k, gamma_k, sigma_k,
%                  lamda_k, and xi_k with respect to theta as function of theta.
%                  of theta.
%       For mixture of gamma distribution:
%           f(x) = sum_k
%           w_k*beta^alpha/Gamma(alpha)*x^(alpha-1)*exp(-beta*x)
%          .w_fun{k}
%          .alpha_fun{k}
%          .beta_fun{k}
%             ... mixture parameters w_k, alpha_k, and beta_k
%                 as function of theta (e.g., w_k = w_fun{k}(theta).
%          .dwdtheta_fun{k}
%          .dalphadtheta_fun{k}
%          .dbetadtheta_fun{k}
%              ... derivatives of mixture parameters w_k, alpha_k,
%                 and beta_k with respect to theta as function of theta.
%                  of theta.
% D{i} ... information about i-th experiment
%     .name ... name of experiment
%     .data ... data
%         .uncensored ... column vector containing uncencored data
% 	      .censored ... column vector containing cencored data
%     .observation_interval ... inter-observation time
% options ...
%     .scale ... option determining whether positive or negative values
%         are returned (see below). The default is s = 0.
%     .grad_ind ... index of parameters with respect to which the gradient
%         is computed.
%
% Outputs:
% ========
% For options.sign = 'positive':
%   logL ... log-likelihood of the data given the parameters theta
%   grad ... gradient of the log-likelihood of the data given the parameters
%   F    ... Fisher information matrix (or something closely related)
% For options.sign = 'negative':
%   logL ... negative log-likelihood of the data given the parameters theta
%   grad ... gradient of the negative log-likelihood of the data given the
%            parameters
%   F    ... negative Fisher information matrix (or something closely related)
% (The negative variantes are often require for optimizations, as by default
%  minimization problems are considered.)
%
% 2012/05/16 Jan Hasenauer

% function [logL,grad] = logLikelihoodMM(theta,M,D,options)
function [varargout] = logLikelihoodMM(varargin)

if nargin >= 3
    theta = varargin{1};
    M = varargin{2};
    D = varargin{3};
else
    error('Not enought inputs.')
end

% Check and assign options
options.sign = 'positive';
options.grad_ind = [1:length(theta)]';
if nargin == 4
    if isfield(varargin{4},'sign')
        options.sign = varargin{4}.sign;
    end
    if isfield(varargin{4},'grad_ind')
        options.grad_ind = varargin{4}.grad_ind;
    end
    
    %    options = setdefault(varargin{4},options);
end

%% INITIALIZATION
theta = theta(:);
logL = 0;
grad = zeros(length(options.grad_ind),1);

% Loop: experiments
for j = 1:length(D)
    
    %% RESTORE DATA
    x  = D{j}.data.uncensored(:);
    xc = D{j}.data.censored(:);
    dx = D{j}.observation_interval;
    
    Px  = zeros(length(x),1);
    Pxc = zeros(length(xc),1);
    
    % Scaling
    if isfield(M.experiment(j),'scaling_fun')
        x  = M.experiment(j).scaling_fun(theta)*x;
        xc = M.experiment(j).scaling_fun(theta)*xc;
        dx = M.experiment(j).scaling_fun(theta)*dx;
    end
    
    % Loop: number of mixture elements per experimental condition
    for k = 1:M.experiment(j).size
        
        switch M.mixture.type
            case 'normal'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                mk = M.experiment(j).mu_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                
                yk  = (x-mk)/sk;
                if dx > 0
                    yk_ = ((x-dx)-mk)/sk;
                end
                ykc = (xc-mk)/sk;
                
            case 'log-normal'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                mk = M.experiment(j).mu_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                yk  = (log(x)-mk)/sk;
                if dx > 0
                    yk_ = (log(x-dx)-mk)/sk;
                end
                ykc = (log(xc)-mk)/sk;
                
            case 'Johnson SU'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                gk = M.experiment(j).gamma_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                lk = M.experiment(j).lambda_fun{k}(theta);
                xk = M.experiment(j).xi_fun{k}(theta);
                
                zk  = (x-xk)/lk;
                yk  = gk + sk*asinh((x-xk)/lk);
                if dx > 0
                    yk_ = gk + sk*asinh((x-dx-xk)/lk);
                end
                ykc = gk + sk*asinh((xc-xk)/lk);
                
            case 'gamma'
                
                %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                wk = M.experiment(j).w_fun{k}(theta);
                ak = M.experiment(j).alpha_fun{k}(theta);
                bk = M.experiment(j).beta_fun{k}(theta);
                
            otherwise
                error('Unknown distribution type. Only ''normal'', ''log-normal'', ''Johnson SU'' and ''gamma'' are allowed!');
                
        end
        
        % Ensure positivity of weights
        if wk <= 0
            switch  options.sign
                case 'positive'
                    varargout{1} = -inf;
                    varargout{2} = nan(size(options.grad_ind));
                case 'negative'
                    varargout{1} =  inf;
                    varargout{2} = nan(size(options.grad_ind));
            end
            return;
        end
        
        %% EVALUATION OF MIXTURE PROBABILITY OF DATA POINTS
        switch M.mixture.type
            case 'normal'
                % Uncensored data
                if length(x) >= 0
                    if dx > 0 % discrete observations
                        Px  = Px  + wk*(    0.5*erfc(-yk ./sqrt(2))...
                            - 0.5*erfc(-yk_./sqrt(2)));
                    else % continuous observations
                        Px  = Px  + wk * exp(-0.5*yk.^2)./(sqrt(2*pi)*sk);
                    end
                end
                
                % Censored data
                if length(xc) >= 1
                    Pxc = Pxc + wk*(1 - 0.5*erfc(-ykc./sqrt(2)));
                end
                
            case 'log-normal'
                % Uncensored data
                if length(x) >= 0
                    if dx > 0 % discrete observations
                        Px  = Px  + wk*(    0.5*erfc(-yk ./sqrt(2))...
                            - 0.5*erfc(-yk_./sqrt(2)));
                    else % continuous observations
                        Px  = Px  + wk * exp(-0.5*yk.^2)./(sqrt(2*pi)*sk*x);
                    end
                end
                
                % Censored data
                if length(xc) >= 1
                    Pxc = Pxc + wk*(1 - 0.5*erfc(-ykc./sqrt(2)));
                end
                
            case 'Johnson SU'
                % Uncensored data
                if length(x) >= 0
                    if dx > 0 % discrete observations
                        Px  = Px  + wk*(    0.5*erfc(-yk ./sqrt(2))...
                            - 0.5*erfc(-yk_./sqrt(2)));
                    else % continuous observations
                        Px  = Px  + wk * exp(-0.5*yk.^2).*(sk./(lk*sqrt(2*pi)*sqrt(zk.^2+1)));
                    end
                end
                
                % Censored data
                if length(xc) >= 1
                    Pxc = Pxc + wk*(1 - 0.5*erfc(-ykc./sqrt(2)));
                end
                
            case 'gamma'
                % Uncensored data
                if length(x) >= 0
                    if dx > 0 % discrete observations
                        Px  = Px  + wk*(gamcdf(x,ak,1/bk) - gamcdf(x-dx,ak,1/bk));
                    else % continuous observations
                        Px  = Px  + wk* gampdf(x,ak,1/bk);
                    end
                end
                
                % Censored data
                if length(xc) >= 1
                    Pxc = Pxc + wk*(1 - gamcdf(xc,ak,1/bk));
                end
        end
    end
    
    %% EVALUATION OF LOG-LIKELIHOOD FUNCTION
    if length(x) >= 1
        logL = logL + sum(log(Px ));
    end
    if length(xc) >= 1
        logL = logL + sum(log(Pxc));
    end
    
    %% GRADIENT CALCULATION
    if nargout >= 2
        
        % Loop: mixture elements
        for k = 1:M.experiment(j).size
            
            switch M.mixture.type
                case 'normal'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk    = M.experiment(j).w_fun{k}(theta);
                    dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                    dwkdt = dwkdt(options.grad_ind);
                    mk    = M.experiment(j).mu_fun{k}(theta);
                    dmkdt = M.experiment(j).dmudtheta_fun{k}(theta);
                    dmkdt = dmkdt(options.grad_ind);
                    sk    = M.experiment(j).sigma_fun{k}(theta);
                    dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                    dskdt = dskdt(options.grad_ind);
                    
                    yk  = (x-mk)'/sk;
                    if dx > 0
                        yk_ = ((x-dx)-mk)'/sk;
                    end
                    ykc = (xc-mk)'/sk;
                    
                    %% EVALUATION OF GRADIENT
                    if length(x) >= 1
                        if dx > 0
                            S = + dwkdt*(( 0.5*erfc(-yk ./sqrt(2))  ...
                                -0.5*erfc(-yk_./sqrt(2)))./Px') ...
                                + wk/sk*(-bsxfun(@times,bsxfun(@plus,dskdt*yk ,dmkdt),1/sqrt(2*pi)*exp(-0.5*yk .^2)./Px')  ...
                                +bsxfun(@times,bsxfun(@plus,dskdt*yk_,dmkdt),1/sqrt(2*pi)*exp(-0.5*yk_.^2)./Px'));
                        else
                            
                            S = + dwkdt*(exp(-0.5*yk.^2)./(sqrt(2*pi)*sk*(Px)')) ...
                                + wk*bsxfun(@times,bsxfun(@plus,dskdt*(yk.^2-1),dmkdt*yk),exp(-0.5*yk.^2)./(sqrt(2*pi)*sk^2*(Px)'));
                        end
                        grad = grad + sum(S,2);
                    end
                    if length(xc) >= 1
                        S = - dwkdt*(0.5*erfc(-ykc./sqrt(2))./Pxc') ...
                            + wk/sk*bsxfun(@times,bsxfun(@plus,dskdt*ykc,dmkdt),1/sqrt(2*pi)*exp(-0.5*ykc.^2)./Pxc');
                        grad = grad + sum(S,2);
                    end
                    
                case 'log-normal'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk    = M.experiment(j).w_fun{k}(theta);
                    dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                    dwkdt = dwkdt(options.grad_ind);
                    mk    = M.experiment(j).mu_fun{k}(theta);
                    dmkdt = M.experiment(j).dmudtheta_fun{k}(theta);
                    dmkdt = dmkdt(options.grad_ind);
                    sk    = M.experiment(j).sigma_fun{k}(theta);
                    dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                    dskdt = dskdt(options.grad_ind);
                    
                    yk  = (log(x)-mk)'/sk;
                    if dx > 0
                        yk_ = (log(x-dx)-mk)'/sk;
                    end
                    ykc = (log(xc)-mk)'/sk;
                    
                    %% EVALUATION OF GRADIENT
                    if length(x) >= 1
                        if dx > 0
                            S = + dwkdt*(( 0.5*erfc(-yk ./sqrt(2))  ...
                                -0.5*erfc(-yk_./sqrt(2)))./Px') ...
                                + wk/sk*(-bsxfun(@times,bsxfun(@plus,dskdt*yk ,dmkdt),1/sqrt(2*pi)*exp(-0.5*yk .^2)./Px')  ...
                                +bsxfun(@times,bsxfun(@plus,dskdt*yk_,dmkdt),1/sqrt(2*pi)*exp(-0.5*yk_.^2)./Px'));
                        else
                            S = + dwkdt*(exp(-0.5*yk.^2)./(sqrt(2*pi)*sk*(x.*Px)')) ...
                                + wk*bsxfun(@times,bsxfun(@plus,dskdt*(yk.^2-1),dmkdt*yk),exp(-0.5*yk.^2)./(sqrt(2*pi)*sk^2*(x.*Px)'));
                            
                        end
                        S(isnan(S))=0;
                        grad = grad + sum(S,2);
                    end
                    if length(xc) >= 1
                        S = - dwkdt*(0.5*erfc(-ykc./sqrt(2))./Pxc') ...
                            + wk/sk*bsxfun(@times,bsxfun(@plus,dskdt*ykc,dmkdt),1/sqrt(2*pi)*exp(-0.5*ykc.^2)./Pxc');
                        grad = grad + sum(S,2);
                    end
                    
                case 'Johnson SU'
                    
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk    = M.experiment(j).w_fun{k}(theta);
                    dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                    dwkdt = dwkdt(options.grad_ind);
                    gk    = M.experiment(j).gamma_fun{k}(theta);
                    dgkdt = M.experiment(j).dgammadtheta_fun{k}(theta);
                    dgkdt = dgkdt(options.grad_ind);
                    sk    = M.experiment(j).sigma_fun{k}(theta);
                    dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                    dskdt = dskdt(options.grad_ind);
                    lk    = M.experiment(j).lambda_fun{k}(theta);
                    dlkdt = M.experiment(j).dlambdadtheta_fun{k}(theta);
                    dlkdt = dlkdt(options.grad_ind);
                    xk    = M.experiment(j).xi_fun{k}(theta);
                    dxkdt = M.experiment(j).dxidtheta_fun{k}(theta);
                    dxkdt = dxkdt(options.grad_ind);
                    
                    zk  = (x-xk)'/lk;
                    if dx > 0
                        zk_ = (x-dx-xk)'/lk;
                    end
                    zkc = (xc-xk)'/lk;
                    
                    yk  = gk + sk*asinh(zk );
                    if dx > 0
                        yk_ = gk + sk*asinh(zk_);
                    end
                    ykc = gk + sk*asinh(zkc);
                    
                    %% EVALUATION OF GRADIENT
                    if length(x) >= 1
                        if dx > 0
                            S =  dwkdt*(( 0.5*erfc(-yk ./sqrt(2))       ...
                                -0.5*erfc(-yk_./sqrt(2)))./Px') ...
                                + wk*( bsxfun(@times,...
                                bsxfun(@plus,...
                                - dlkdt*(sk*zk ./(lk*sqrt(1+zk .^2))) ...
                                - dxkdt*(sk    ./(lk*sqrt(1+zk .^2))) ...
                                + dskdt*asinh(zk ), ...
                                dgkdt),...
                                1/sqrt(2*pi)*exp(-0.5*yk .^2)./Px') ...
                                -bsxfun(@times,...
                                bsxfun(@plus,...
                                - dlkdt*(sk*zk_./(lk*sqrt(1+zk_.^2))) ...
                                - dxkdt*(sk    ./(lk*sqrt(1+zk_.^2))) ...
                                + dskdt*asinh(zk_), ...
                                dgkdt),...
                                1/sqrt(2*pi)*exp(-0.5*yk_.^2)./Px'));
                        else
                            S = +dwkdt*((sk./(lk*sqrt(1+zk .^2))).*(1/sqrt(2*pi).*exp(-0.5*yk .^2)./Px'))...
                                + wk*( bsxfun(@times,...
                                bsxfun(@plus,...
                                + dlkdt*((sk*zk.^2 +sk^2*zk.*yk.*sqrt(1+zk.^2)-sk*(zk.^2+1))./(lk^2*(1+zk.^2) .^(3/2))) ...
                                + dxkdt*((sk*zk+sk^2*yk.*sqrt(1+zk.^2))./(lk^2*(1+zk.^2).^(3/2))) ...
                                + dskdt*((1-sk*asinh(zk).*yk)./(lk*sqrt(1+zk.^2))), ...
                                - dgkdt*(sk*yk./(lk*sqrt(1+zk.^2)))),...
                                1/sqrt(2*pi)*exp(-0.5*yk .^2)./Px'));
                            
                        end
                        grad = grad + sum(S,2);
                    end
                    if length(xc) >= 1
                        S = - dwkdt*(0.5*erfc(-ykc./sqrt(2))./Pxc') ...
                            - wk*bsxfun(@times,...
                            bsxfun(@plus,...
                            - dlkdt*(sk*zkc./(lk*sqrt(1+zkc.^2))) ...
                            - dxkdt*(sk    ./(lk*sqrt(1+zkc.^2))) ...
                            + dskdt*asinh(zkc), ...
                            dgkdt),...
                            1/sqrt(2*pi)*exp(-0.5*ykc.^2)./Pxc');
                        grad = grad + sum(S,2);
                    end
                    
                case 'gamma'
                    %% ASSIGN PARAMETERS AND AUXILIARY VARIABLES
                    wk    = M.experiment(j).w_fun{k}(theta);
                    dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                    dwkdt = dwkdt(options.grad_ind);
                    ak    = M.experiment(j).alpha_fun{k}(theta);
                    dakdt = M.experiment(j).dalphadtheta_fun{k}(theta);
                    dakdt = dakdt(options.grad_ind);
                    bk    = M.experiment(j).beta_fun{k}(theta);
                    dbkdt = M.experiment(j).dbetadtheta_fun{k}(theta);
                    dbkdt = dbkdt(options.grad_ind);
                    
                    %% EVALUATION OF GRADIENT
                    if length(x) >= 1
                        if dx > 0
                            S =  dwkdt*((gamcdf(x,ak,1/bk)' - gamcdf(x-dx,ak,1/bk)')./Px') ...
                                + wk*((...
                                (dakdt*...
                                (((-(gamma(ak))^2*(bk*x).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*x)/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*x,ak+zeros(length(x),1)).*log(bk*x))...
                                ./gamma(ak)...
                                -gamcdf(x,ak,1/bk)*psi(ak))'./Px'))...
                                +dbkdt*(((1/gamma(ak)*(bk*x).^(ak-1).*exp(-bk*x).*x))'./(Px)'))...
                                - ((dakdt*...
                                (((-(gamma(ak))^2*(bk*(x-dx)).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*(x-dx))/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*(x-dx),ak+zeros(length(x),1)).*log(bk*(x-dx)))./gamma(ak)-gamcdf(x-dx,ak,1/bk)*psi(ak))'./Px'))...
                                +dbkdt*(((1/gamma(ak)*(bk*(x-dx)).^(ak-1).*exp(-bk*(x-dx)).*(x-dx)))'./(Px)')));
                            
                            
                        else
                            S = + dwkdt*(gampdf(x,ak,1/bk)'./(Px)')...% + wk*([(gampdf(x,ak,1/bk)./(Px))';(gampdf(x,ak,1/bk)./(Px))'].*(dakdt*(log(x)+log(bk)-psi(ak))'+dbkdt*((ak\bk)-x)'));
                                + wk*bsxfun(@times,bsxfun(@plus,dakdt*(log(x)+log(bk)-psi(ak))',dbkdt*((ak/bk)-x)'),(gampdf(x,ak,1/bk)./(Px))');
                            
                        end
                        grad = grad + sum(S,2);
                        
                    end
                    if length(xc) >= 1
                        S= - (dwkdt*(gamcdf(xc,ak,1/bk)'./(Pxc)') ...
                            + wk*((dakdt*(((-(gamma(ak))^2*(bk*xc).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*xc)/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*xc,ak+zeros(length(xc),1)).*log(bk*xc))./gamma(ak)-gamcdf(xc,ak,1/bk)*psi(ak))'./Pxc'))...
                            +dbkdt*(((1/gamma(ak)*(bk*xc).^(ak-1).*exp(-bk*xc).*xc))'./(Pxc)')));
                        grad = grad + sum(S,2);
                    end
            end
            
        end
        
    end
    
end

%% CONSTRUCTION OF OUTPUT
switch nargout
    % One output
    case {0,1}
        switch  options.sign
            case 'positive'
                varargout{1} =  logL;
            case 'negative'
                varargout{1} = -logL;
        end
        % Two outputs
    case 2
        switch  options.sign
            case 'positive'
                varargout{1} =  logL;
                varargout{2} =  grad;
            case 'negative'
                varargout{1} = -logL;
                varargout{2} = -grad;
        end
end
