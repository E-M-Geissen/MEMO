% logLikelihoodMMc evaluates the log-likelihood of a mixture model M
%   with respect to a certain dataset D. The mixture model is parametrized
%   using the parameters theta. This paramterization is rather flexible
%   and allows for complex dependencies of the mixture components on theta.
%
% USAGE:
% ======
% [output] = logLikelihoodMMc(theta,M,Mc,D)
% [output] = logLikelihoodMMc(theta,M,Mc,D,options)
% [logL]   = logLikelihoodMMc(...)
% [logL,grad] = logLikelihoodMMc(...)
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
%       For mixture of log-normals:
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
%       For mixture of gamma distributions:
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
% Mc ... model for censoring distribution (same structure as M).
%   the same distributions types as for M are supported.
%   Additionally a single delta distribution is supported.
%           wk=1;
%           Mc.experiment(j).location;
%
% D{i} ... information about i-th experiment
%     .name ... name of experiment
%     .data ... data
%         .uncensored ... column vector containing uncencored data
% 	      .censored ... column vector containing cencored data
%         .cen_type ... type of censoring for data in .censored 'left' or
%                   ... 'right', if not specified default is right
%     .observation_interval ... inter-observation time
% options ...
%     .scale ... option determining whether positive or negative values
%         are return (see below). The default is s = 0.
%     .grad_ind ... index of parameters with respect to which the gradient
%         is computed.
%
% Outputs:
% ========
% For options.sign = 'positive':
%   logL ... log-likelihood of the data given the parameters theta
%   grad ... gradient of the log-likelihood of the data given the parameters
% For options.sign = 'negative':
%   logL ... negative log-likelihood of the data given the parameters theta
%   grad ... gradient of the negative log-likelihood of the data given the
%            parameters
% (The negative variantes are often require for optimizations, as by default
%  minimization problems are considered.)
%
% 2012/07/04 Jan Hasenauer
% modified Eva-Maria Geissen 13.02.2016


% function [logL,grad] = logLikelihoodMMc(theta,M,Mc,D,options)
function [varargout] = logLikelihoodMMc(varargin)

if nargin >= 3
    theta = varargin{1};
    M = varargin{2};
    Mc = varargin{3};
    D = varargin{4};
else
    error('Not enought inputs.')
end

% Check and assign options
options.sign = 'positive';
options.grad_ind = (1:length(theta))';
options.refinement = 4;
options.refinement_max = 1000;
% options.quadrature_type = 'trapeziodal';
options.quadrature_type = 'Simpson';
if nargin == 5
    if isfield(varargin{5},'sign')
        options.sign = varargin{5}.sign;
    end
    if isfield(varargin{5},'grad_ind')
        options.grad_ind = varargin{5}.grad_ind;
    end
    if isfield(varargin{5},'refinement')
        options.refinement = varargin{5}.refinement;
    end
    if isfield(varargin{5},'refinement_max')
        options.refinement_max = varargin{5}.refinement_max;
    end
end

if D{1}.observation_interval > 0
    
    % Select refinement
    options.refinement = max(options.refinement,1);
    switch options.quadrature_type
        case 'trapeziodal'
            if mod(options.refinement,2) ~= 0
                options.refinement = options.refinement + 1;
            end
        case 'Simpson'
            if mod(options.refinement,4) ~= 0
                options.refinement = options.refinement + (4 - mod(options.refinement,4));
            end
        otherwise
            error('This quadrature formular is not available.')
    end
    
    %% INITIALIZATION
    theta = theta(:);
    n_grad_theta = length(options.grad_ind);
    logL = 0;
    if nargout >= 2 % Gradient computed
        grad = zeros(n_grad_theta,1);
    end
    r = options.refinement;
    if r > options.refinement_max
        warning('Numerical problems: Maximal refinement reached.');
    end
    refine = 'false';
    
    
    % Loop: experiments
    for j = 1:length(D)
        
        %% RESTORE DATA
        if ~isfield(D{j},'unique_values')
            % Assignment
            x  = D{j}.data.uncensored(:);
            xc = D{j}.data.censored(:);
            dx = D{j}.observation_interval/r;
            % Integration grid
            [xi0,~,I] = unique([x;xc])';
            ind  = I(1:length(x))';
            indc = I(length(x)+1:end)';
        else
            % Assignment
            dx = D{j}.observation_interval/r;
            % Integration grid
            xi0 = D{j}.unique_values;
            ind  = D{j}.position_uncensored';
            indc = D{j}.position_censored'  ;
        end
        xi  = bsxfun(@plus,xi0,  dx*(-r:0))';   xi  = xi(:)';
        xir = bsxfun(@plus,xi0,  dx*(-r:2:0))'; xir = xir(:)';
        n_xi  = length(xi);
        n_xir = length(xir);
        n_x  = length(xi)/(r+1);
        % Test index mapping
        x  = D{j}.data.uncensored(:);
        xc = D{j}.data.censored(:);
        % [x-xi(ind)']
        % [x-xi(ind_)-D{j}.observation_interval]
        % [xc-xi(indc)]
        % [xc-xi(indc_)-D{j}.observation_interval]
        
        % Initialization
        f   = zeros(1,n_xi );
        fr  = zeros(1,n_xir);
        F   = zeros(1,n_xi );
        Fr  = zeros(1,n_xir);
        fc  = zeros(1,n_xi );
        fcr = zeros(1,n_xir);
        Fc  = zeros(1,n_xi );
        Fcr = zeros(1,n_xir);
        if nargout >= 2 % Gradient computed
            dfdt  = zeros(n_grad_theta,n_xi);
            dFdt  = zeros(n_grad_theta,n_xi);
            dfcdt = zeros(n_grad_theta,n_xi);
            dFcdt = zeros(n_grad_theta,n_xi);
        end
        
        % Loop: mixture elements uncensored data
        for k = 1:M.experiment(j).size
            switch M.mixture.type
                case 'normal'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                    wk = M.experiment(j).w_fun{k}(theta);
                    mk = M.experiment(j).mu_fun{k}(theta);
                    sk = M.experiment(j).sigma_fun{k}(theta);
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dmkdt = M.experiment(j).dmudtheta_fun{k}(theta);
                        dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dmkdt = sparse(dmkdt(options.grad_ind));
                        dskdt = sparse(dskdt(options.grad_ind));
                    end
                    
                    %% ASSIGNMENT OF AUXILIARY VARIABLES
                    yki  = (xi-mk)/sk;
                    ykir = (xir-mk)/sk;
                    
                    %% DENSITY CALCULATION
                    phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                    phikr = exp(-0.5*ykir.^2)./sqrt(2*pi);
                    fk  =   phik ./(sk);
                    fkr =   phikr./(sk);
                    Fk  = 0.5*erfc(-yki ./sqrt(2));
                    Fkr = 0.5*erfc(-ykir./sqrt(2));
                    
                    %% GRADIENT CALCULATION
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                            + dskdt*(wk/sk*((yki.^2-1).*fk)));
                        dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                            + dskdt*(wk/sk*(yki.*phik)));
                    end
                    
                case 'log-normal'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                    wk = M.experiment(j).w_fun{k}(theta);
                    mk = M.experiment(j).mu_fun{k}(theta);
                    sk = M.experiment(j).sigma_fun{k}(theta);
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dmkdt = M.experiment(j).dmudtheta_fun{k}(theta);
                        dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dmkdt = sparse(dmkdt(options.grad_ind));
                        dskdt = sparse(dskdt(options.grad_ind));
                    end
                    
                    %% ASSIGNMENT OF AUXILIARY VARIABLES
                    yki  = (log(xi )-mk)/sk;
                    ykir = (log(xir)-mk)/sk;
                    
                    %% DENSITY CALCULATION
                    phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                    phikr = exp(-0.5*ykir.^2)./sqrt(2*pi);
                    fk  =   phik ./(sk*xi);
                    fkr =   phikr./(sk*xir);
                    Fk  = 0.5*erfc(-yki ./sqrt(2));
                    Fkr = 0.5*erfc(-ykir./sqrt(2));
                    
                    %% GRADIENT CALCULATION
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                            + dskdt*(wk/sk*((yki.^2-1).*fk)));
                        dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                            + dskdt*(wk/sk*(yki.*phik)));
                    end
                    
                case 'Johnson SU'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                    wk = M.experiment(j).w_fun{k}(theta);
                    gk = M.experiment(j).gamma_fun{k}(theta);
                    sk = M.experiment(j).sigma_fun{k}(theta);
                    lk = M.experiment(j).lambda_fun{k}(theta);
                    xk = M.experiment(j).xi_fun{k}(theta);
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dgkdt = M.experiment(j).dgammadtheta_fun{k}(theta);
                        dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                        dlkdt = M.experiment(j).dlambdadtheta_fun{k}(theta);
                        dxkdt = M.experiment(j).dxidtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dgkdt = sparse(dgkdt(options.grad_ind));
                        dskdt = sparse(dskdt(options.grad_ind));
                        dlkdt = sparse(dlkdt(options.grad_ind));
                        dxkdt = sparse(dxkdt(options.grad_ind));
                    end
                    
                    %% ASSIGNMENT OF AUXILIARY VARIABLES
                    zki   = (xi  - xk)/lk;
                    zkir  = (xir - xk)/lk;
                    yki   = gk + sk*asinh(zki );
                    ykir  = gk + sk*asinh(zkir);
                    
                    %% DENSITY CALCULATION
                    phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                    phikr = exp(-0.5*ykir.^2)./sqrt(2*pi);
                    fk  =     sk/lk*phik ./(sqrt(zki .^2 + 1));
                    fkr =     sk/lk*phikr./(sqrt(zkir.^2 + 1));
                    Fk  = 0.5*erfc(-yki ./sqrt(2));
                    Fkr = 0.5*erfc(-ykir./sqrt(2));
                    
                    %% GRADIENT CALCULATION
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + (  dgkdt*(wk*phik.*( - sk*yki )./(lk*sqrt(zki.^2+1))) ... % dgamma/dtheta
                            + dskdt*(wk*phik.*(1 - sk*yki.*asinh(zki))./(lk*sqrt(zki.^2+1))) ... % dsigma/dtheta
                            + dlkdt*(wk*phik.*(sk*zki.^2 + sk^2*zki.*yki.*sqrt(zki.^2+1) - sk*(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2))) ... % dlambda/dtheta
                            + dxkdt*(wk*phik.*(sk*zki + sk^2*yki.*sqrt(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2)))); % dxi/dtheta
                        dwkFkdt = dwkdt*Fk + (- dlkdt*(wk*sk/lk*(phik.*zki./sqrt(1+zki.^2))) ...
                            - dxkdt*(wk*sk/lk*(phik     ./sqrt(1+zki.^2))) ...
                            + dgkdt*(wk*       phik                      ) ...
                            + dskdt*(wk*       phik.*asinh(zki)          ));
                    end
                    
                case 'gamma'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES%% ASSIGN PARAMETERS
                    wk = M.experiment(j).w_fun{k}(theta);
                    ak = M.experiment(j).alpha_fun{k}(theta);
                    bk = M.experiment(j).beta_fun{k}(theta);
                    
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dakdt = M.experiment(j).dalphadtheta_fun{k}(theta);
                        dbkdt = M.experiment(j).dbetadtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dakdt = sparse(dakdt(options.grad_ind));
                        dbkdt = sparse(dbkdt(options.grad_ind));
                    end
                    
                    %% DENSITY CALCULATION
                    fk  =   gampdf(xi,ak,1/bk);
                    fkr =   gampdf(xir,ak,1/bk);
                    Fk  =   gamcdf(xi,ak,1/bk);
                    Fkr =   gamcdf(xir,ak,1/bk);
                    
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + wk*((dakdt*(fk.*(log(bk)+log(xi)-psi(ak)))+dbkdt*(fk.*(ak/bk-xi))));
                        dwkFkdt = dwkdt*Fk +wk*(dakdt*(1/gamma(ak)*(-(gamma(ak))^2*(bk*xi).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*xi)/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*xi,ak+zeros(length(xi),1)').*log(bk*xi))...
                            -Fk*psi(ak))+dbkdt*(1/gamma(ak)*(bk*xi).^(ak-1).*exp(-bk*xi).*xi));
                    end
                    
                otherwise
                    error('Unknown distribution type. Only ''normal'', ''log-normal'', ''gamma'' and ''Johnson SU'' are allowed!');
            end
            
            %% SUMMATION
            f  = f  + wk*fk;
            fr = fr + wk*fkr;
            F  = F  + wk*Fk;
            Fr = Fr + wk*Fkr;
            if nargout >= 2 % Gradient computed
                dfdt = dfdt + dwkfkdt;
                dFdt = dFdt + dwkFkdt;
            end
            
        end
        if ~isempty(xc)
            % Loop: mixture elements censored data
            for k = 1:Mc.experiment(j).size
                switch Mc.mixture.type
                    case 'normal'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        mk = Mc.experiment(j).mu_fun{k}(theta);
                        sk = Mc.experiment(j).sigma_fun{k}(theta);
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dmkdt = Mc.experiment(j).dmudtheta_fun{k}(theta);
                            dskdt = Mc.experiment(j).dsigmadtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dmkdt = sparse(dmkdt(options.grad_ind));
                            dskdt = sparse(dskdt(options.grad_ind));
                        end
                        
                        %% ASSIGNMENT OF AUXILIARY VARIABLES
                        yki  = (xi-mk)/sk;
                        ykir = (xir-mk)/sk;
                        
                        %% DENSITY CALCULATION
                        phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                        phikr = exp(-0.5*ykir.^2)./sqrt(2*pi);
                        fk  =   phik ./sk;
                        fkr =   phikr./sk;
                        Fk  = 0.5*erfc(-yki ./sqrt(2));
                        Fkr = 0.5*erfc(-ykir./sqrt(2));
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                                + dskdt*(wk/sk*((yki.^2-1).*fk)));
                            dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                                + dskdt*(wk/sk*(yki.*phik)));
                        end
                        
                    case 'log-normal'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        mk = Mc.experiment(j).mu_fun{k}(theta);
                        sk = Mc.experiment(j).sigma_fun{k}(theta);
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dmkdt = Mc.experiment(j).dmudtheta_fun{k}(theta);
                            dskdt = Mc.experiment(j).dsigmadtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dmkdt = sparse(dmkdt(options.grad_ind));
                            dskdt = sparse(dskdt(options.grad_ind));
                        end
                        
                        %% ASSIGNMENT OF AUXILIARY VARIABLES
                        yki  = (log(xi )-mk)/sk;
                        ykir = (log(xir)-mk)/sk;
                        
                        %% DENSITY CALCULATION
                        phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                        phikr = exp(-0.5*ykir.^2)./sqrt(2*pi);
                        fk  =   phik ./(sk*xi);
                        fkr =   phikr./(sk*xir);
                        Fk  = 0.5*erfc(-yki ./sqrt(2));
                        Fkr = 0.5*erfc(-ykir./sqrt(2));
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                                + dskdt*(wk/sk*((yki.^2-1).*fk)));
                            dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                                + dskdt*(wk/sk*(yki.*phik)));
                        end
                        
                    case 'Johnson SU'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        gk = Mc.experiment(j).gamma_fun{k}(theta);
                        sk = Mc.experiment(j).sigma_fun{k}(theta);
                        lk = Mc.experiment(j).lambda_fun{k}(theta);
                        xk = Mc.experiment(j).xi_fun{k}(theta);
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dgkdt = Mc.experiment(j).dgammadtheta_fun{k}(theta);
                            dskdt = Mc.experiment(j).dsigmadtheta_fun{k}(theta);
                            dlkdt = Mc.experiment(j).dlambdadtheta_fun{k}(theta);
                            dxkdt = Mc.experiment(j).dxidtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dgkdt = sparse(dgkdt(options.grad_ind));
                            dskdt = sparse(dskdt(options.grad_ind));
                            dlkdt = sparse(dlkdt(options.grad_ind));
                            dxkdt = sparse(dxkdt(options.grad_ind));
                        end
                        
                        %% ASSIGNMENT OF AUXILIARY VARIABLES
                        zki   = (xi  - xk)/lk;
                        zkir  = (xir - xk)/lk;
                        yki   = gk + sk*asinh(zki );
                        ykir  = gk + sk*asinh(zkir);
                        
                        %% DENSITY CALCULATION
                        phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                        phikr = exp(-0.5*ykir.^2)./sqrt(2*pi);
                        fk  = sk/lk*phik ./(sqrt(zki .^2 + 1));
                        fkr = sk/lk*phikr./(sqrt(zkir.^2 + 1));
                        Fk  = 0.5*erfc(-yki ./sqrt(2));
                        Fkr = 0.5*erfc(-ykir./sqrt(2));
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + (  dgkdt*(wk*phik.*( - sk*yki )./(lk*sqrt(zki.^2+1))) ... % dgamma/dtheta
                                + dskdt*(wk*phik.*(1 - sk*yki.*asinh(zki))./(lk*sqrt(zki.^2+1))) ... % dsigma/dtheta
                                + dlkdt*(wk*phik.*(sk*zki.^2 + sk^2*zki.*yki.*sqrt(zki.^2+1) - sk*(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2))) ... % dlambda/dtheta
                                + dxkdt*(wk*phik.*(sk*zki + sk^2*yki.*sqrt(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2)))); % dxi/dtheta
                            dwkFkdt = dwkdt*Fk + (- dlkdt*(wk*sk/lk*(phik.*zki./sqrt(1+zki.^2))) ...
                                - dxkdt*(wk*sk/lk*(phik     ./sqrt(1+zki.^2))) ...
                                + dgkdt*(wk*       phik                      ) ...
                                + dskdt*(wk*       phik.*asinh(zki)          ));
                        end
                        
                    case 'gamma'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES%% ASSIGN PARAMETERS
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        ak = Mc.experiment(j).alpha_fun{k}(theta);
                        bk = Mc.experiment(j).beta_fun{k}(theta);
                        
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dakdt = Mc.experiment(j).dalphadtheta_fun{k}(theta);
                            dbkdt = Mc.experiment(j).dbetadtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dakdt = sparse(dakdt(options.grad_ind));
                            dbkdt = sparse(dbkdt(options.grad_ind));
                        end
                        
                        %% DENSITY CALCULATION
                        fk  =   gampdf(xi,ak,1/bk);
                        fkr =   gampdf(xir,ak,1/bk);
                        Fk  =   gamcdf(xi,ak,1/bk);
                        Fkr =   gamcdf(xir,ak,1/bk);
                        
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + wk*((dakdt*(fk.*(log(bk)+log(xi)-psi(ak)))+dbkdt*(fk.*(ak/bk-xi))));
                            dwkFkdt = dwkdt*Fk +wk*(dakdt*(1/gamma(ak)*(-(gamma(ak))^2*(bk*xi).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*xi)/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*xi,ak+zeros(length(xi),1)').*log(bk*xi))...
                                -Fk*psi(ak))+dbkdt*(1/gamma(ak)*(bk*xi).^(ak-1).*exp(-bk*xi).*xi));
                        end
                        
                    case 'delta'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        %wk = Mc.experiment(j).w_fun{k}(theta);
                        wk=1;
                        Tk = Mc.experiment(j).location;
                        
                        
                        %% DENSITY CALCULATION
                        fk(xi==Tk)=1;
                        fk(xi~=Tk)=0;%1e-12;
                        fkr(xir==Tk)=1;
                        fkr(xir~=Tk)=0;%1e-12;
                        Fk(xi<Tk)  = 0;
                        Fk(xi>=Tk)  = 1;
                        Fkr(xir<Tk)  = 0;
                        Fkr(xir>=Tk)  = 1;
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = 0;
                            dwkFkdt = 0;
                        end
                        
                    otherwise
                        error('Unknown distribution type. Only ''normal'', ''log-normal'', ''gamma'' and ''Johnson SU'' are allowed!');
                end
                
                %% SUMMATION
                fc  = fc  + wk*fk;
                fcr = fcr + wk*fkr;
                Fc  = Fc  + wk*Fk;
                Fcr = Fcr + wk*Fkr;
                
                
                if nargout >= 2 % Gradient computed
                    dfcdt = dfcdt + dwkfkdt;
                    dFcdt = dFcdt + dwkFkdt;
                end
                
            end
        end
        %% LIKELIHOOD
        
        if  isfield(D{j},'data.cen_type')
            switch D{j}.data.cen_type
                case 'right'
                    % Event distribution
                    fA  = f .*(1-Fc);
                    fB  = fc.*(1-F );
                    fAr  = fr .*(1-Fcr);
                    fBr  = fcr.*(1-Fr );
                case 'left'
                    % Event distribution
                    fA  = f .*(Fc);
                    fB  = fc.*(F );
                    fAr  = fr .*(Fcr);
                    fBr  = fcr.*(Fr );
            end
        else
            % Event distribution
            fA  = f .*(1-Fc);
            fB  = fc.*(1-F );
            fAr  = fr .*(1-Fcr);
            fBr  = fcr.*(1-Fr );
        end
        
        
        
        %if dx > 0
        % Cumulative event probability
        switch options.quadrature_type
            case 'trapeziodal'
                % Weight
                w  =   dx*[0.5,ones(1,r-1),0.5];
                wr = 2*dx*[0.5,ones(1,r/2-1),0.5];
            case 'Simpson'
                % Weight
                w = ones(1,r+1);
                w(2:2:end) = 4;
                w(3:2:end-2) = 2;
                w = dx/3*w;
                wr = ones(1,r/2+1);
                wr(2:2:end) = 4;
                wr(3:2:end-2) = 2;
                wr = 2*dx/3*wr;
        end
        % Integration
        FA  = w *reshape(fA ,r  +1,n_x);
        FAr = wr*reshape(fAr,r/2+1,n_x);
        
        if ~isempty(Mc)
            if strcmp(Mc.mixture.type, 'delta')
                w_neu=zeros(length(w),1)';
                w_neu(end)=1;
                wr_neu=zeros(length(wr),1)';
                wr_neu(end)=1;
                FB  = w_neu*reshape(fB ,r  +1,n_x);
                FBr = wr_neu*reshape(fBr,r/2+1,n_x);
            else
                FB  = w *reshape(fB ,r  +1,n_x);
                FBr = wr*reshape(fBr,r/2+1,n_x);
            end
        else
            FB  = w *reshape(fB ,r  +1,n_x);
            FBr = wr*reshape(fBr,r/2+1,n_x);
        end
        
        % Probability of individual data points
        Px   = FA(ind );
        Pxc  = FB(indc);
        Pxr  = FAr(ind );
        Pxcr = FBr(indc);
        
        % Likelihood
        if ~isempty(Px)
            logL = logL + sum(log(Px));
            
        end
        if ~isempty(Pxc)
            logL = logL + sum(log(Pxc));
            
        end
        
        %% ERROR CONTROL
        % Error
        Err_max = max(max(abs(Px  - Pxr )/(Px+1e-10)),max(abs(Pxc - Pxcr)/(Pxc+1e-10)));
        % Check for error threshold
        if Err_max >= 0.1
            refine = 'true';
            break;
        end
        
        % figure;
        %
        % subplot(2,2,1);
        % plot(xi,f,'b-'); hold on;
        % plot(xir,fr,'b-x'); hold on;
        % plot(xi,fc,'r-'); hold on;
        % plot(xir,fcr,'r-x'); hold on;
        % legend('f','f_r','fc','fc_r');
        %
        % subplot(2,2,2);
        % plot(xi,fA,'b-'); hold on;
        % plot(xir,fAr,'b-x'); hold on;
        % plot(xi,fB,'r-'); hold on;
        % plot(xir,fBr,'r-x'); hold on;
        % legend('fA','fA_r','fB','fB_r');
        %
        % subplot(2,2,3);
        % plot(x,Px,'bo'); hold on;
        % plot(x,Pxr,'bs'); hold on;
        % legend('Px','Px_r');
        %
        % subplot(2,2,4);
        % plot(xc,Pxc,'bo'); hold on;
        % plot(xc,Pxcr,'bs'); hold on;
        % legend('Pxc','Pxc_r');
        
        
        %% GRADIENT
        if nargout >= 2 % Gradient computed
            
            
            if  isfield(D{j},'data.cen_type')
                switch D{j}.data.cen_type
                    case 'right'
                        % Event distribution
                        dfAdt = [zeros(n_grad_theta,1),bsxfun(@times,dfdt(:,2:end),(1-Fc(2:end))) - bsxfun(@times,f(2:end),dFcdt(:,2:end))];
                        dfBdt = [zeros(n_grad_theta,1),bsxfun(@times,dfcdt(:,2:end),(1-F(2:end))) - bsxfun(@times,fc(2:end),dFdt(:,2:end))];
                    case 'left'
                        % Event distribution
                        dfAdt = [zeros(n_grad_theta,1),bsxfun(@times,dfdt(:,2:end),(Fc(2:end))) + bsxfun(@times,f(2:end),dFcdt(:,2:end))];
                        dfBdt = [zeros(n_grad_theta,1),bsxfun(@times,dfcdt(:,2:end),(F(2:end))) + bsxfun(@times,fc(2:end),dFdt(:,2:end))];
                end
            else
                % Event distribution
                dfAdt = [zeros(n_grad_theta,1),bsxfun(@times,dfdt(:,2:end),(1-Fc(2:end))) - bsxfun(@times,f(2:end),dFcdt(:,2:end))];
                dfBdt = [zeros(n_grad_theta,1),bsxfun(@times,dfcdt(:,2:end),(1-F(2:end))) - bsxfun(@times,fc(2:end),dFdt(:,2:end))];
            end
            
            
            % Cumulative event probability
            switch options.quadrature_type
                case 'trapeziodal'
                    % Weight
                    w  =   dx*[0.5,ones(1,r-1),0.5];
                    % Integration
                    dFAdt  = squeeze(sum(bsxfun(@times,reshape(dfAdt,[n_grad_theta,r+1,n_x]),w),2));
                    dFBdt  = squeeze(sum(bsxfun(@times,reshape(dfBdt,[n_grad_theta,r+1,n_x]),w),2));
                case 'Simpson'
                    % Weight
                    w = ones(1,r+1);
                    w(2:2:end) = 4;
                    w(3:2:end-2) = 2;
                    w = dx/3*w;
                    % Integration
                    dFAdt  = squeeze(sum(bsxfun(@times,reshape(dfAdt,[n_grad_theta,r+1,n_x]),w),2));
                    dFBdt  = squeeze(sum(bsxfun(@times,reshape(dfBdt,[n_grad_theta,r+1,n_x]),w),2));
            end
            % Probability of individual data points
            dPxdt  = dFAdt(:,ind );
            dPxcdt = dFBdt(:,indc);
            % Likelihood
            if ~isempty(Px)
                grad = grad + dPxdt*(Px.^-1)';
            end
            if ~isempty(Pxc)
                grad = grad + dPxcdt*(Pxc.^-1)';
            end
        end
        
    end
    
    %% CHECK IF REFINEMENT IS NECESSARY
    if strcmp(refine,'true')
        % Assign modified options
        options_refine.refinement = 2*options.refinement;
        options_refine.sign = 'positive';
        options_refine.grad_ind = options.grad_ind;
        % Check output setup:
        switch nargout
            % One output
            case {0,1}
                logL = logLikelihoodMMc(theta,M,Mc,D,options_refine);
                % Two outputs
            case 2
                [logL,grad] = logLikelihoodMMc(theta,M,Mc,D,options_refine);
        end
    end
    
    
else
    
    %% INITIALIZATION
    theta = theta(:);
    n_grad_theta = length(options.grad_ind);
    logL = 0;
    if nargout >= 2 % Gradient computed
        grad = zeros(n_grad_theta,1);
    end
    
    
    % Loop: experiments
    for j = 1:length(D)
        
        %% RESTORE DATA
        if ~isfield(D{j},'unique_values')
            % Assignment
            x  = D{j}.data.uncensored(:);
            xc = D{j}.data.censored(:);
            dx = D{j}.observation_interval;
            
            [xi,~,I] = unique([x;xc]');
            ind  = I(1:length(x))';
            indc = I(length(x)+1:end)';
        else
            % Assignment
            dx = D{j}.observation_interval;
            % Integration grid
            xi= D{j}.unique_values';
            ind  = D{j}.position_uncensored';
            indc = D{j}.position_censored'  ;
        end
        
        n_xi  = length(xi);
        n_x  = length(xi);
        % Test index mapping
        x  = D{j}.data.uncensored(:);
        xc = D{j}.data.censored(:);
        % [x-xi(ind)']
        % [x-xi(ind_)-D{j}.observation_interval]
        % [xc-xi(indc)]
        % [xc-xi(indc_)-D{j}.observation_interval]
        
        % Initialization
        f   = zeros(1,n_x  );
        F   = zeros(1,n_x  );
        fc  = zeros(1,n_x );
        Fc  = zeros(1,n_x );
        
        if nargout >= 2 % Gradient computed
            dfdt  = zeros(n_grad_theta,n_x );
            dFdt  = zeros(n_grad_theta,n_x );
            dfcdt = zeros(n_grad_theta,n_x );
            dFcdt = zeros(n_grad_theta,n_x );
        end
        
        % Loop: mixture elements uncensored data
        for k = 1:M.experiment(j).size
            switch M.mixture.type
                case 'normal'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                    wk = M.experiment(j).w_fun{k}(theta);
                    mk = M.experiment(j).mu_fun{k}(theta);
                    sk = M.experiment(j).sigma_fun{k}(theta);
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dmkdt = M.experiment(j).dmudtheta_fun{k}(theta);
                        dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dmkdt = sparse(dmkdt(options.grad_ind));
                        dskdt = sparse(dskdt(options.grad_ind));
                    end
                    
                    %% ASSIGNMENT OF AUXILIARY VARIABLES
                    yki  = (xi-mk)/sk;
                    
                    %% DENSITY CALCULATION
                    phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                    fk  =   phik ./(sk);
                    Fk  = 0.5*erfc(-yki ./sqrt(2));
                    
                    %% GRADIENT CALCULATION
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                            + dskdt*(wk/sk*((yki.^2-1).*fk)));
                        dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                            + dskdt*(wk/sk*(yki.*phik)));
                    end
                    
                case 'log-normal'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                    wk = M.experiment(j).w_fun{k}(theta);
                    mk = M.experiment(j).mu_fun{k}(theta);
                    sk = M.experiment(j).sigma_fun{k}(theta);
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dmkdt = M.experiment(j).dmudtheta_fun{k}(theta);
                        dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dmkdt = sparse(dmkdt(options.grad_ind));
                        dskdt = sparse(dskdt(options.grad_ind));
                    end
                    
                    %% ASSIGNMENT OF AUXILIARY VARIABLES
                    yki  = (log(xi )-mk)/sk;
                    
                    %% DENSITY CALCULATION
                    phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                    fk  =   phik ./(sk*xi);
                    Fk  = 0.5*erfc(-yki ./sqrt(2));
                    
                    %% GRADIENT CALCULATION
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                            + dskdt*(wk/sk*((yki.^2-1).*fk)));
                        dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                            + dskdt*(wk/sk*(yki.*phik)));
                    end
                    
                case 'Johnson SU'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                    wk = M.experiment(j).w_fun{k}(theta);
                    gk = M.experiment(j).gamma_fun{k}(theta);
                    sk = M.experiment(j).sigma_fun{k}(theta);
                    lk = M.experiment(j).lambda_fun{k}(theta);
                    xk = M.experiment(j).xi_fun{k}(theta);
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dgkdt = M.experiment(j).dgammadtheta_fun{k}(theta);
                        dskdt = M.experiment(j).dsigmadtheta_fun{k}(theta);
                        dlkdt = M.experiment(j).dlambdadtheta_fun{k}(theta);
                        dxkdt = M.experiment(j).dxidtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dgkdt = sparse(dgkdt(options.grad_ind));
                        dskdt = sparse(dskdt(options.grad_ind));
                        dlkdt = sparse(dlkdt(options.grad_ind));
                        dxkdt = sparse(dxkdt(options.grad_ind));
                    end
                    
                    %% ASSIGNMENT OF AUXILIARY VARIABLES
                    zki   = (xi  - xk)/lk;
                    yki   = gk + sk*asinh(zki );
                    
                    %% DENSITY CALCULATION
                    phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                    fk  =     sk/lk*phik ./(sqrt(zki .^2 + 1));
                    Fk  = 0.5*erfc(-yki ./sqrt(2));
                    
                    %% GRADIENT CALCULATION
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + (  dgkdt*(wk*phik.*( - sk*yki )./(lk*sqrt(zki.^2+1))) ... % dgamma/dtheta
                            + dskdt*(wk*phik.*(1 - sk*yki.*asinh(zki))./(lk*sqrt(zki.^2+1))) ... % dsigma/dtheta
                            + dlkdt*(wk*phik.*(sk*zki.^2 + sk^2*zki.*yki.*sqrt(zki.^2+1) - sk*(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2))) ... % dlambda/dtheta
                            + dxkdt*(wk*phik.*(sk*zki + sk^2*yki.*sqrt(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2)))); % dxi/dtheta
                        dwkFkdt = dwkdt*Fk + (- dlkdt*(wk*sk/lk*(phik.*zki./sqrt(1+zki.^2))) ...
                            - dxkdt*(wk*sk/lk*(phik     ./sqrt(1+zki.^2))) ...
                            + dgkdt*(wk*       phik                      ) ...
                            + dskdt*(wk*       phik.*asinh(zki)          ));
                    end
                case 'gamma'
                    
                    %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES%% ASSIGN PARAMETERS
                    wk = M.experiment(j).w_fun{k}(theta);
                    ak = M.experiment(j).alpha_fun{k}(theta);
                    bk = M.experiment(j).beta_fun{k}(theta);
                    
                    if nargout >= 2 % Gradient computed
                        % Evaluation of gradients
                        dwkdt = M.experiment(j).dwdtheta_fun{k}(theta);
                        dakdt = M.experiment(j).dalphadtheta_fun{k}(theta);
                        dbkdt = M.experiment(j).dbetadtheta_fun{k}(theta);
                        % Elimination of unnecessary parameter directions
                        dwkdt = sparse(dwkdt(options.grad_ind));
                        dakdt = sparse(dakdt(options.grad_ind));
                        dbkdt = sparse(dbkdt(options.grad_ind));
                    end
                    
                    %% DENSITY CALCULATION
                    fk  =   gampdf(xi,ak,1/bk);
                    Fk  =   gamcdf(xi,ak,1/bk);
                    
                    if nargout >= 2 % Gradient computed
                        dwkfkdt = dwkdt*fk + wk*((dakdt*(fk.*(log(bk)+log(xi)-psi(ak)))+dbkdt*(fk.*(ak/bk-xi))));
                        dwkFkdt = dwkdt*Fk +wk*(dakdt*(1/gamma(ak)*(-(gamma(ak))^2*(bk*xi).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*xi)/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*xi,ak+zeros(length(xi),1)').*log(bk*xi))...
                            -Fk*psi(ak))+dbkdt*(1/gamma(ak)*(bk*xi).^(ak-1).*exp(-bk*xi).*xi));
                    end
                    
                otherwise
                    error('Unknown distribution type. Only ''normal'', ''log-normal'', ''gamma'' and ''Johnson SU'' are allowed!');
            end
            
            %% SUMMATION
            f  = f  + wk*fk;
            F  = F  + wk*Fk;
            if nargout >= 2 % Gradient computed
                dfdt = dfdt + dwkfkdt;
                dFdt = dFdt + dwkFkdt;
            end
            
        end
        
        if ~isempty(xc)
            % Loop: mixture elements censored data
            for k = 1:Mc.experiment(j).size
                switch Mc.mixture.type
                    case 'log-normal'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        mk = Mc.experiment(j).mu_fun{k}(theta);
                        sk = Mc.experiment(j).sigma_fun{k}(theta);
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dmkdt = Mc.experiment(j).dmudtheta_fun{k}(theta);
                            dskdt = Mc.experiment(j).dsigmadtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dmkdt = sparse(dmkdt(options.grad_ind));
                            dskdt = sparse(dskdt(options.grad_ind));
                        end
                        
                        %% ASSIGNMENT OF AUXILIARY VARIABLES
                        yki  = (log(xi )-mk)/sk;
                        
                        %% DENSITY CALCULATION
                        phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                        fk  =   phik ./(sk*xi);
                        Fk  = 0.5*erfc(-yki ./sqrt(2));
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + (  dmkdt*(wk/sk*yki.*fk) ...
                                + dskdt*(wk/sk*((yki.^2-1).*fk)));
                            dwkFkdt = dwkdt*Fk - (  dmkdt*(wk/sk*phik) ...
                                + dskdt*(wk/sk*(yki.*phik)));
                        end
                        
                    case 'Johnson SU'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        gk = Mc.experiment(j).gamma_fun{k}(theta);
                        sk = Mc.experiment(j).sigma_fun{k}(theta);
                        lk = Mc.experiment(j).lambda_fun{k}(theta);
                        xk = Mc.experiment(j).xi_fun{k}(theta);
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dgkdt = Mc.experiment(j).dgammadtheta_fun{k}(theta);
                            dskdt = Mc.experiment(j).dsigmadtheta_fun{k}(theta);
                            dlkdt = Mc.experiment(j).dlambdadtheta_fun{k}(theta);
                            dxkdt = Mc.experiment(j).dxidtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dgkdt = sparse(dgkdt(options.grad_ind));
                            dskdt = sparse(dskdt(options.grad_ind));
                            dlkdt = sparse(dlkdt(options.grad_ind));
                            dxkdt = sparse(dxkdt(options.grad_ind));
                        end
                        
                        %% ASSIGNMENT OF AUXILIARY VARIABLES
                        zki   = (xi  - xk)/lk;
                        yki   = gk + sk*asinh(zki );
                        
                        %% DENSITY CALCULATION
                        phik  = exp(-0.5*yki .^2)./sqrt(2*pi);
                        fk  = sk/lk*phik ./(sqrt(zki .^2 + 1));
                        Fk  = 0.5*erfc(-yki ./sqrt(2));
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + (  dgkdt*(wk*phik.*( - sk*yki )./(lk*sqrt(zki.^2+1))) ... % dgamma/dtheta
                                + dskdt*(wk*phik.*(1 - sk*yki.*asinh(zki))./(lk*sqrt(zki.^2+1))) ... % dsigma/dtheta
                                + dlkdt*(wk*phik.*(sk*zki.^2 + sk^2*zki.*yki.*sqrt(zki.^2+1) - sk*(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2))) ... % dlambda/dtheta
                                + dxkdt*(wk*phik.*(sk*zki + sk^2*yki.*sqrt(zki.^2+1))./(lk^2*(zki.^2+1).^(3/2)))); % dxi/dtheta
                            dwkFkdt = dwkdt*Fk + (- dlkdt*(wk*sk/lk*(phik.*zki./sqrt(1+zki.^2))) ...
                                - dxkdt*(wk*sk/lk*(phik     ./sqrt(1+zki.^2))) ...
                                + dgkdt*(wk*       phik                      ) ...
                                + dskdt*(wk*       phik.*asinh(zki)          ));
                        end
                        
                    case 'gamma'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES%% ASSIGN PARAMETERS
                        wk = Mc.experiment(j).w_fun{k}(theta);
                        ak = Mc.experiment(j).alpha_fun{k}(theta);
                        bk = Mc.experiment(j).beta_fun{k}(theta);
                        
                        if nargout >= 2 % Gradient computed
                            % Evaluation of gradients
                            dwkdt = Mc.experiment(j).dwdtheta_fun{k}(theta);
                            dakdt = Mc.experiment(j).dalphadtheta_fun{k}(theta);
                            dbkdt = Mc.experiment(j).dbetadtheta_fun{k}(theta);
                            % Elimination of unnecessary parameter directions
                            dwkdt = sparse(dwkdt(options.grad_ind));
                            dakdt = sparse(dakdt(options.grad_ind));
                            dbkdt = sparse(dbkdt(options.grad_ind));
                        end
                        
                        %% DENSITY CALCULATION
                        fk  =   gampdf(xi,ak,1/bk);
                        Fk  =   gamcdf(xi,ak,1/bk);
                        
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = dwkdt*fk + wk*((dakdt*(fk.*(log(bk)+log(xi)-psi(ak)))+dbkdt*(fk.*(ak/bk-xi))));
                            dwkFkdt = dwkdt*Fk +wk*(dakdt*(1/gamma(ak)*(-(gamma(ak))^2*(bk*xi).^ak.*(hypergeom([ak,ak],[ak+1, ak+1],-bk*xi)/gamma(ak+1)^2)+gamma(ak)*gammainc(bk*xi,ak+zeros(length(xi),1)').*log(bk*xi))...
                                -Fk*psi(ak))+dbkdt*(1/gamma(ak)*(bk*xi).^(ak-1).*exp(-bk*xi).*xi));
                        end
                        
                    case 'delta'
                        
                        %% ASSIGNMENT OF MIXTURE PARAMETERS AND DERIVATIVES
                        %                    wk = Mc.experiment(j).w_fun{k}(theta);
                        wk=1;
                        Tk = Mc.experiment(j).location;
                        
                        
                        %% DENSITY CALCULATION
                        
                        
                        fk(xi==Tk)=1;
                        fk(xi~=Tk)=0;
                        
                        
                        Fk(xi>=Tk)=1;
                        Fk(xi<Tk)=0;
                        
                        %% GRADIENT CALCULATION
                        if nargout >= 2 % Gradient computed
                            dwkfkdt = zeros(n_grad_theta,n_x );
                            dwkFkdt =zeros(n_grad_theta,n_x );
                        end
                        
                    otherwise
                        error('Unknown distribution type. Only ''normal'', ''log-normal'', ''gamma'', ''Johnson SU'' and ''delta'' are allowed!');
                end
                
                %% SUMMATION
                fc  = fc  + wk*fk;
                
                Fc  = Fc  + wk*Fk;
                
                if ~strcmp(Mc.mixture.type,'delta')
                    if nargout >= 2 % Gradient computed
                        dfcdt = dfcdt + dwkfkdt;
                        dFcdt = dFcdt + dwkFkdt;
                    end
                end
                
            end
        end
        %% LIKELIHOOD
        if  isfield(D{j},'data.cen_type')
            switch D{j}.data.cen_type
                case 'right'
                    % Event distribution
                    fA  = f .*(1-Fc);
                    fB  = fc.*(1-F );
                    
                    Px   = fA(ind );
                    Pxc  = fB(indc);
                case 'left'
                    % Event distribution
                    fA  = f .*(Fc);
                    fB  = fc.*(F);
                    
                    Px   = fA(ind );
                    Pxc  = fB(indc);
            end
        else
            % Event distribution
            fA  = f .*(1-Fc);
            fB  = fc.*(1-F );
            
            Px   = fA(ind );
            Pxc  = fB(indc);
        end
        
        
        % Likelihood
        if ~isempty(Px)
            logL = logL + sum(log(Px));
        end
        if ~isempty(Pxc)
            logL = logL + sum(log(Pxc));
        end
        
        
        
        % figure;
        %
        % subplot(2,2,1);
        % plot(xi,f,'b-'); hold on;
        % plot(xir,fr,'b-x'); hold on;
        % plot(xi,fc,'r-'); hold on;
        % plot(xir,fcr,'r-x'); hold on;
        % legend('f','f_r','fc','fc_r');
        %
        % subplot(2,2,2);
        % plot(xi,fA,'b-'); hold on;
        % plot(xir,fAr,'b-x'); hold on;
        % plot(xi,fB,'r-'); hold on;
        % plot(xir,fBr,'r-x'); hold on;
        % legend('fA','fA_r','fB','fB_r');
        %
        % subplot(2,2,3);
        % plot(x,Px,'bo'); hold on;
        % plot(x,Pxr,'bs'); hold on;
        % legend('Px','Px_r');
        %
        % subplot(2,2,4);
        % plot(xc,Pxc,'bo'); hold on;
        % plot(xc,Pxcr,'bs'); hold on;
        % legend('Pxc','Pxc_r');
        
        
        %% GRADIENT
        
        if nargout >= 2 % Gradient computed
            
            if  isfield(D{j},'data.cen_type')
                switch D{j}.data.cen_type
                    case 'right'
                        % Event distribution
                        dfAdt = [bsxfun(@times,dfdt(:,1:end),(1-Fc(1:end))) - bsxfun(@times,f(1:end),dFcdt(:,1:end))];
                        dfBdt = [bsxfun(@times,dfcdt(:,1:end),(1-F(1:end))) - bsxfun(@times,fc(1:end),dFdt(:,1:end))];
                    case 'left'
                        % Event distribution
                        dfAdt = [bsxfun(@times,dfdt(:,1:end),(Fc(1:end))) + bsxfun(@times,f(1:end),dFcdt(:,1:end))];
                        dfBdt = [bsxfun(@times,dfcdt(:,1:end),(F(1:end))) + bsxfun(@times,fc(1:end),dFdt(:,1:end))];
                end
            else
                % Event distribution
                dfAdt = [bsxfun(@times,dfdt(:,1:end),(1-Fc(1:end))) - bsxfun(@times,f(1:end),dFcdt(:,1:end))];
                dfBdt = [bsxfun(@times,dfcdt(:,1:end),(1-F(1:end))) - bsxfun(@times,fc(1:end),dFdt(:,1:end))];
            end
            
            
            dFAdt  = dfAdt;
            dFBdt  = dfBdt;
            
            % Probability of individual data points
            dPxdt  = dFAdt(:,ind );
            dPxcdt = dFBdt(:,indc);
            % Likelihood
            if ~isempty(Px)
                grad = grad + dPxdt*(Px.^-1)';
            end
            if ~isempty(Pxc)
                grad = grad + dPxcdt*(Pxc.^-1)';
            end
        end
        
    end
    
    
end

% Consideration of measurement error in input variable
if isfield(M.experiment(1),'cond')
    for j = 1:length(D)
        for k=1:length(M.experiment(j).cond.e)
            sig = M.experiment(j).cond.sigma{k} ;
            logL = logL - 0.5*(log(2*pi*sig^2)+(M.experiment(j).cond.e_fun{k}(theta)/sig)^2);
            if nargout >= 2 % Gradient computed
                grad = grad - 0.5*M.experiment(j).cond.e_fun{k}(theta)/sig^2*M.experiment(j).cond.dedtheta_fun{k}(theta);
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

