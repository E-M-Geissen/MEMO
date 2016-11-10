% getSymbolicParameterVector determines the symbolic variables used to 
% 	define the model M. This is useful to convert the expressions to
% 	functions.
%
% USAGE:
% ======
% theta = getSymbolicParameterVector(M)
%
% INPUTS:
% =======
% theta ... parameter value for model
% M ... mixture model for given data
%
% Outputs:
% ========
% theta ... symbolic parameter vector
%
% 2012/05/25 Jan Hasenauer

function theta = getSymbolicParameterVector(M)

theta = [];
% Loop: Experiments
for j = 1:length(M.experiment)
    % Loop: Mixture components
    for k = 1:M.experiment(j).size
        % Select mixture type
        switch M.mixture.type
            case 'log-normal'
                theta = [theta; ...
                         symvar(M.experiment(j).w{k}).';...
                         symvar(M.experiment(j).mu{k}).';...
                         symvar(M.experiment(j).sigma{k}).'  ];
        case 'Johnson SU'
                theta = [theta; ...
                         symvar(M.experiment(j).w{k}).';...
                         symvar(M.experiment(j).gamma{k}).';...
                         symvar(M.experiment(j).sigma{k}).';...
                         symvar(M.experiment(j).lambda{k}).';...
                         symvar(M.experiment(j).xi{k}).'     ];
        end
    end
end
% Eliminate douple entries
theta = unique(theta);