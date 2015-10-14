% eval_performance calcualtes the AIC and BIC for a 
% maximum posterior estimate
%
% USAGE:
% ======
% 
% eval_performance(D,parameters)
%
% INPUTS:
% =======
%  D ... data structure
%  parameters ... struct containg information about parameters, at least:
%   .number number of parameters
%   .MS.MAP.logPost ... log-posterior at the maximum a posterior estimate
%
%  2015/09/23 Eva-Maria Geissen

function parameters = eval_performance(D, parameters)
%% EVALUATION OF THE NUMBER OF DATA POINTS
N = 0;
for j = 1:length(D)
    N = N + length(D{j}.data.censored) + length(D{j}.data.uncensored);
end

% evaluate performance
parameters.AIC  = -2*parameters.MS.MAP.logPost + 2*parameters.number;
parameters.BIC  = -2*parameters.MS.MAP.logPost + parameters.number*log(N);