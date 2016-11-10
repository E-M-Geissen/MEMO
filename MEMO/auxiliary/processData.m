% process data processes the experimental data for faster computation
%
% USAGE:
% ======
% processData(D)
% 
%
% INPUTS:
% =======
% D{i} ... information about i-th experiment
%     .name ... name of experiment
%     .data ... data
%         .uncensored ... column vector containing uncencored data
% 	      .censored ... column vector containing cencored data
%     .observation_interval ... inter-observation time
%
% Outputs:
% ========
% D{i} ... updated data object containing:
%    .unique_values ... unique values in data set
%    .position_uncensored ... positions of uncensored data points in
%                             unique_values
%    .position_censored ... positions of right censored data points in
%                            unique_values

function D = processData(D)

for i = 1:length(D)
    % Determine unique position
    [xi,~,I] = unique([D{i}.data.uncensored;D{i}.data.censored]);
    % Restore processed data
    D{i}.unique_values = xi;
    D{i}.position_uncensored = I(1:length(D{i}.data.uncensored));
    D{i}.position_censored = I(length(D{i}.data.uncensored)+1:end);
end