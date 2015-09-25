function parameters = getConfidenceIntervals(parameters,alpha)

for i = 1:parameters.number
    % Confidence intervals computed using local approximation and a
    % threshold (-> similar to PL-based confidence intervals)
    parameters.CI.local_PL(i,1) = parameters.MS.MAP.par(i) - sqrt(icdf('chi2',alpha,1)/parameters.MS.MAP.hessian(i,i));
    parameters.CI.local_PL(i,2) = parameters.MS.MAP.par(i) + sqrt(icdf('chi2',alpha,1)/parameters.MS.MAP.hessian(i,i));

    % Confidence intervals computed using local approximation and the
    % probability mass (-> similar to Bayesian confidence intervals)
    if parameters.MS.MAP.hessian(i,i) > 1e-16
        parameters.CI.local_B(i,1)  = icdf('norm',  (1-alpha)/2,parameters.MS.MAP.par(i),inv(sqrt(parameters.MS.MAP.hessian(i,i))));
        parameters.CI.local_B(i,2)  = icdf('norm',1-(1-alpha)/2,parameters.MS.MAP.par(i),inv(sqrt(parameters.MS.MAP.hessian(i,i))));
    else
        parameters.CI.local_B(i,1) = -inf;
        parameters.CI.local_B(i,2) =  inf;
    end
    
    % Confidence intervals computed using profile likelihood
    if isfield(parameters,'P')
    if i <= length(parameters.P)
        % left bound
        ind  = find(parameters.P(i).par(i,:) <= parameters.MS.MAP.par(i));
        j = find(parameters.P(i).R(ind) <= exp(-icdf('chi2',alpha,1)/2),1,'last');
        if ~isempty(j)
            parameters.CI.PL(i,1) = interp1(parameters.P(i).R(ind([j,j+1])),...
                parameters.P(i).par(i,ind([j,j+1])),exp(-icdf('chi2',alpha,1)/2));
        else
            parameters.CI.PL(i,1) = -inf;
        end
        % right bound
        ind  = find(parameters.P(i).par(i,:) >= parameters.MS.MAP.par(i));
        j = find(parameters.P(i).R(ind) <= exp(-icdf('chi2',alpha,1)/2),1,'first');
        if ~isempty(j)
            parameters.CI.PL(i,2) = interp1(parameters.P(i).R(ind([j-1,j])),...
                parameters.P(i).par(i,ind([j-1,j])),exp(-icdf('chi2',alpha,1)/2));
        else
            parameters.CI.PL(i,2) = -inf;
        end
    end
    end
    
end
