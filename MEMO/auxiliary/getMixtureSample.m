% getMixtureSample samples from the provided mixture model.
%
% USAGE:
% ======
% [Sample] = getMixtureSample(n,theta,M)
% [Sample] = getMixtureSample(n,theta,M,options)
%
% INPUTS:
% =======
% n ... sample size (either scalar or vector of number of experiments)
% theta ... parameter value for model
% M ... mixture model for given data
% options ... sampling options
%   .min ... lower bound for sampled values (default = -inf). 
%   .max ... upper bound for sampled values (default =  inf). 
%
% Outputs:
% ========
% Sample ... sample of event and censoring times
%    {j} ... j-th condition (experiment).
%
% 2012/07/02 Jan Hasenauer

% function Dens = getMixtureSample(n,theta,M,options)
function Sample = getMixtureSample(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 3
    n = varargin{1};
    theta = varargin{2};
    M = varargin{3};
else
    error('Not enought inputs.')
end
% Check dimension
if length(n) == 1
    n = n*ones(length(M.experiment),1);
elseif length(n) ~= length(M.experiment)
    error('Lenght of n and number of experiments must agree (or length(n) = 1).')
end

% Options
options.min = -inf;
options.max =  inf;
if nargin == 4
    options = setdefault(varargin{4},options);
end

% Number of sampled points
if (options.min == -inf) && (options.max == inf)
    N = n;
else
    N = 5*n;
end

%% SAMPLING OF EVENT
% Loop: Experiments
for j = 1:length(M.experiment)    
    X  = nan(N(j),1);
    w0 = 0;
    I0 = 0;
    switch M.mixture.type
        case 'log-normal'
            % Draw uniform random number to determine subpopulation
            R = rand(N(j),1);
            % Loop: Mixture components
            for k = 1:M.experiment(j).size
                % Assign parameters
                wk = M.experiment(j).w_fun{k}(theta); 
                mk = M.experiment(j).mu_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                % Determine number of cells in i-th sumpopulation
                Ik = sum((w0 <= R).*(R < (w0+wk)));
                % Sample event time of cells in this subpopulation
                X(I0+1:I0+Ik) = exp(sk*randn(Ik,1) + mk);
                % Update for next subpopulation
                w0 = w0 + wk;
                I0 = I0 + Ik;
            end            
            
        case 'Johnson SU'
            % Draw uniform random number to determine subpopulation
            R = rand(N(j),1);
            % Loop: Mixture components
            for k = 1:M.experiment(j).size
                % Assign parameters
                wk = M.experiment(j).w_fun{k}(theta);
                gk = M.experiment(j).gamma_fun{k}(theta);
                sk = M.experiment(j).sigma_fun{k}(theta);
                lk = M.experiment(j).lambda_fun{k}(theta);
                xk = M.experiment(j).xi_fun{k}(theta);
                % Determine number of cells in i-th sumpopulation
                Ik = sum((w0 <= R).*(R < (w0+wk)));
                % Sample event time of cells in this subpopulation
                X(I0+1:I0+Ik) = xk + lk*sinh((randn(Ik,1)-gk)/sk);
                % Update for next subpopulation
                w0 = w0 + wk;
                I0 = I0 + Ik;
            end            
        otherwise
            error('This option is not available.');
    end
    
    % Assignment
    if (options.min == -inf) && (options.max == inf)
        Sample{j} = X(randperm(n(j)));
    else
        X = X(randperm(N(j)));
        X = X(X>=options.min);
        X = X(X<=options.max);
        if length(X) >= n(j)
            Sample{j} = X(1:n(j));
        else
            warning('Probability distribution has large mass outside of support.')
        end
    end
end

