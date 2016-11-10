% getMixtureModel completes a given mixture model, by construction of function
%   handles to determine the parameter of the mixture from hyperparameters
%   theta. Merely, symbolic expression of the mixture parameters are
%   required to construct corresponding function handles and function
%   handles
%   which allow the evaluation of the derivative.
%
% USAGE:
% ======
% M = getMixtureModel(M,theta)
% [M,C] = getMixtureModel(M,theta)
%
% INPUTS:
% =======
% M ... mixture model containing symbolic experssions
%    .mixture.type ... type of mixture
%    .experiment{j} ... j-th experiment
%       .name ... name of experiment
%       .size ... number of mixture components
%       For mixture of normals:
%           f(x) = sum_k w_k*Phi((x-mu_k)/sigma_k)
%          .w{k}
%          .mu{k} 
%          .sigma{k} 
%             ... symbolic expression for mixture parameters w_k, mu_k and 
%                 sigma_k as function of the symbolic elements of theta.
%       For mixture of log-normals:
%           f(x) = sum_k w_k*Phi((log(x)-mu_k)/sigma_k)
%          .w{k}
%          .mu{k} 
%          .sigma{k} 
%             ... symbolic expression for mixture parameters w_k, mu_k and 
%                 sigma_k as function of the symbolic elements of theta.
%       For mixture of Johnson SU distributions:
%           f(x) = sum_k w_k*Phi(gamma_k + sigma_k*arcsinh((x-xi_k)/lambda_k))
%          .w{k}
%          .alpha{k}
%          .beta{k} 
%             ... symbolic expression for mixture parameters w_k, alpha_k, 
%                 and beta_k as function of symbolic elements of theta.
%       For mixture of gamma distribution:
%           f(x) = sum_k w_k*beta^alpha/Gamma(alpha)*x^(alpha-1)*exp(-beta*x)
%          .w_fun{k}
%          .alpha_fun{k} ... shape parameter
%          .beta_fun{k} ... rate parameter
% theta ... vector of parameters (all as symbolics)
%
% Outputs:
% ========
% M ... expanded mixture model containg function handles for the
%       computation of the parameters and their derivative with respect to
%       theta.
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
%
% 2012/05/16 Jan Hasenauer

% function M = getMixtureModel(M,theta)
function [M,C] = getMixtureModel(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 2
    M = varargin{1};
    theta = varargin{2};
else
    error('getMixtureModel requires two input arguments.');
end

%% CONSTRUCT MIXTURE MODEL

% Loop: Experiments
for j = 1:length(M.experiment)
    % Scaling
    if isfield(M.experiment(j),'scaling')
        M.experiment(j).scaling_fun = sym2fun(M.experiment(j).scaling,theta);
        M.experiment(j).dscalingdtheta_fun = sym2fun(transpose(simplify(jacobian(M.experiment(j).scaling,theta))),theta);
    end
    % Loop: Mixture components
    for k = 1:M.experiment(j).size
        switch M.mixture.type
            case {'normal','log-normal'}
                %% COMPUTATION OF FUNCTION
                M.experiment(j).w_fun{k}     = sym2fun(M.experiment(j).w{k}    ,theta);
                M.experiment(j).mu_fun{k}    = sym2fun(M.experiment(j).mu{k}   ,theta);
                M.experiment(j).sigma_fun{k} = sym2fun(M.experiment(j).sigma{k},theta);
                %% COMPUTATION OF DERIVATIVE
                M.experiment(j).dwdtheta_fun{k}     = sym2fun(transpose(simplify(jacobian(M.experiment(j).w{k}    ,theta))),theta);
                M.experiment(j).dmudtheta_fun{k}    = sym2fun(transpose(simplify(jacobian(M.experiment(j).mu{k}   ,theta))),theta);
                M.experiment(j).dsigmadtheta_fun{k} = sym2fun(transpose(simplify(jacobian(M.experiment(j).sigma{k},theta))),theta);

            case 'Johnson SU'
                %% COMPUTATION OF FUNCTION
                M.experiment(j).w_fun{k}      = sym2fun(M.experiment(j).w{k}     ,theta);
                M.experiment(j).gamma_fun{k}  = sym2fun(M.experiment(j).gamma{k} ,theta);
                M.experiment(j).sigma_fun{k}  = sym2fun(M.experiment(j).sigma{k} ,theta);
                M.experiment(j).lambda_fun{k} = sym2fun(M.experiment(j).lambda{k},theta);
                M.experiment(j).xi_fun{k}     = sym2fun(M.experiment(j).xi{k}    ,theta);
                %% COMPUTATION OF DERIVATIVE
                M.experiment(j).dwdtheta_fun{k}      = sym2fun(transpose(simplify(jacobian(M.experiment(j).w{k}     ,theta))),theta);
                M.experiment(j).dgammadtheta_fun{k}  = sym2fun(transpose(simplify(jacobian(M.experiment(j).gamma{k} ,theta))),theta);
                M.experiment(j).dsigmadtheta_fun{k}  = sym2fun(transpose(simplify(jacobian(M.experiment(j).sigma{k} ,theta))),theta);
                M.experiment(j).dlambdadtheta_fun{k} = sym2fun(transpose(simplify(jacobian(M.experiment(j).lambda{k},theta))),theta);
                M.experiment(j).dxidtheta_fun{k}     = sym2fun(transpose(simplify(jacobian(M.experiment(j).xi{k}    ,theta))),theta);

            case 'gamma'
                %% COMPUTATION OF FUNCTION
                M.experiment(j).w_fun{k}     = sym2fun(M.experiment(j).w{k}    ,theta);
                M.experiment(j).alpha_fun{k} = sym2fun(M.experiment(j).alpha{k},theta);
                M.experiment(j).beta_fun{k}  = sym2fun(M.experiment(j).beta{k} ,theta);
                %% COMPUTATION OF DERIVATIVE
                M.experiment(j).dwdtheta_fun{k}     = sym2fun(transpose(simplify(jacobian(M.experiment(j).w{k}    ,theta))),theta);
                M.experiment(j).dalphadtheta_fun{k} = sym2fun(transpose(simplify(jacobian(M.experiment(j).alpha{k},theta))),theta);
                M.experiment(j).dbetadtheta_fun{k}  = sym2fun(transpose(simplify(jacobian(M.experiment(j).beta{k} ,theta))),theta);
            
            case 'delta'
                %% COMPUTATION OF FUNCTION
%                 M.experiment(j).w_fun{k}        = sym2fun(M.experiment(j).w{k}       ,theta);
%                 M.experiment(j).location_fun{k} = sym2fun(M.experiment(j).location{k},theta);
%                 %% COMPUTATION OF DERIVATIVE
%                 M.experiment(j).dwdtheta_fun{k}        = sym2fun(transpose(simplify(jacobian(M.experiment(j).w{k}       ,theta))),theta);
%                 M.experiment(j).dlocationdtheta_fun{k} = sym2fun(transpose(simplify(jacobian(M.experiment(j).location{k},theta))),theta);
        end
    end
end

% %% MEASUREMENT ERRORS FOR CONDITION
% if isfield(M.experiment(1),'cond')
%     % Loop: Experiments
%     for j = 1:length(M.experiment)
%         %% COMPUTATION OF FUNCTION
%         M.experiment(j).cond.e_fun = sym2fun(M.experiment(j).cond.e,theta);
%         %% COMPUTATION OF DERIVATIVE
%         M.experiment(j).cond.dedtheta_fun = sym2fun(transpose(simplify(jacobian(M.experiment(j).cond.e,theta))),theta);
%     end
% end

%% MEASUREMENT ERRORS FOR CONDITION
if isfield(M.experiment(1),'cond')
    % Loop: Experiments
    for j = 1:length(M.experiment)
        for k=1:length(M.experiment(j).cond.e)
        %% COMPUTATION OF FUNCTION
        M.experiment(j).cond.e_fun{k} = sym2fun(M.experiment(j).cond.e{k},theta);
        %% COMPUTATION OF DERIVATIVE
        M.experiment(j).cond.dedtheta_fun{k} = sym2fun(transpose(simplify(jacobian(M.experiment(j).cond.e{k},theta))),theta);
        end
    end
end

%% CONSTRUCT CONSTRAINTS

if M.experiment(j).size == 1
    C.A = [];
    C.b = [];
    C.Aeq = [];
    C.beq = [];
    C.c = [];
    C.dcdtheta = [];
    C.ceq = [];
    C.dceqdtheta = [];
else
    % Initialization
    C.ceq = sym([]);
    C.dceqdtheta = sym([]);
    C.c = sym([]);
    C.dcdtheta = sym([]);

    % Loop: Experiments
    for j = 1:length(M.experiment)
        % Loop: Mixture components
        cj = sym(zeros(2*M.experiment(j).size,1));
        dcjdtheta = sym(zeros(2*M.experiment(j).size,length(theta)));
        ceqj = sym(-1);
        for k = 1:M.experiment(j).size
            ceqj = ceqj + M.experiment(j).w{k};
            cj(2*k-1) = -M.experiment(j).w{k};
            cj(2*k  ) =  M.experiment(j).w{k}-1;
            dcjdtheta(2*k-1,:) = transpose(simplify(jacobian(cj(2*k-1),theta)));
            dcjdtheta(2*k  ,:) = transpose(simplify(jacobian(cj(2*k  ),theta)));
        end
        ceqj = simplify(ceqj);
        dceqjdtheta(1,:) = transpose(simplify(jacobian(ceqj,theta)));
        % Assign
        C.c(end+1:end+2*M.experiment(j).size,1) = cj; 
        C.dcdtheta(end+1:end+2*M.experiment(j).size,:) = dcjdtheta;
        C.ceq(end+1,1) = ceqj; 
        C.dceqdtheta(end+1,:) = dceqjdtheta;
    end

    % Remove trivial equality constraints
    I = [];
    for i = 1:length(C.ceq)
        try
            ceqj = double(C.ceq(i));
            if ceqj ~= 0
                I(end+1) = i;
            end
        end
    end
    C.ceq = C.ceq(I);
    C.dceqdtheta = C.dceqdtheta(I,:);

    % Construct constraint matrix - inequality
    I = [];
    C.b = [];
    C.A = [];
    for i = 1:length(C.c)
        try 
            C.A(end+1,:) = double(C.dcdtheta(i,:));
            C.b(end+1,1) = -(C.c(i) - C.A(end,:)*theta);
        catch
            I(end+1) = i;
        end
    end
    C.c = C.c(I);
    C.dcdtheta = C.dcdtheta(I,:);

    % Construct constraint matrix - inequality
    I = [];
    C.beq = [];
    C.Aeq = [];
    for i = 1:length(C.ceq)
        try 
            C.Aeq(end+1,:) = double(C.dceqdtheta(i,:));
            C.beq(end+1,1) = -(C.ceq(i) - C.Aeq(end,:)*theta);
        catch
            I(end+1) = i;
        end
    end
    C.ceq = C.ceq(I);
    C.dceqdtheta = C.dceqdtheta(I,:);

    % Check
    if (length(C.ceq) >= 1) || (length(C.c) >= 1)
        error('Non-linear constraints not implemented yet!');
    end
end
