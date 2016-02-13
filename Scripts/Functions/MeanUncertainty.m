%% FUNCTION TO CALCULATE MEAN OF MULTIPLE VALUES WITH UNCERTAINTY

%% INPUTS
% y = data values
% sigma_y = uncertainty on data values

%% OUTPUTS
% mu = mean value
% sigma_mu = uncertainty in mean

function [mu, sigma_mu] = MeanUncertainty(y, sigma_y)

mu = sum(y./sigma_y.^2)/sum(1./sigma_y.^2); %Bevington and Robinson, Eq. 4.17
sigma_mu = sqrt(1/sum(1./sigma_y.^2)); %Bevington and Robinson, Eq. 4.19