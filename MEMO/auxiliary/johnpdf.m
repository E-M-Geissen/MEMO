% johnpdf computes the probability of the values x under the Johnson SU
%   distribution with parameters gamma, sigma, lambda, and xi
%
% USAGE:
% ======
% f = johnpdf(x,gamma,sigma,lambda,xi)
%
% INPUTS:
% =======
% x ... vector containing the points at which the probability denisty 
%       is evaluated
% gamma ... parameter
% sigma (>0) ... parameter
% lambda (>0) ... parameter
% xi ... parameter
%
% Note: for gamma = 0, sigma = k, lambda = k*s and xi = mu, with k >> 1,
%       the Johnson SU distribution approaches the normal distribution with
%       mean mu and starndard deviation s.
%
% Outputs:
% ========
% f ... vector containing probability densities at x.
%
% 2012/05/16 Jan Hasenauer

% function f = johnpdf(x,gamma,sigma,lambda,xi)
function f = johnpdf(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin == 5
    x = varargin{1};
    gamma = varargin{2};
    sigma = varargin{3};
    lambda = varargin{4};
    xi = varargin{5};
else
    error('This routine requires five input variables.');
end

if (sigma <= 0) || (lambda <= 0)
    error('The parameters sigma and lambda must be positive.');
end

%% CALCULATION OF PROBABILITY DENSITY
z = (x-xi)./lambda;
f = sigma./(lambda*sqrt(2*pi)*sqrt(z.^2+1)) .* exp(-0.5*(gamma+sigma*asinh(z)).^2);
