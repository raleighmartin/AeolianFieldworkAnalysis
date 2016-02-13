%% FUNCTION TO CALCULATE CHI^2 VALUE FOR FIT

%% INPUTS
% y = data values
% sigma_y = uncertainty on data values
% y_fit = fit values
% nu = # degrees of freedom

%% OUTPUTS
% chi2 = chi^2 value
% p-value for chi2

function [chi2, p] = Chi2Calculation(y, sigma_y, y_fit, nu)

chi2 = sum((1./sigma_y.^2).*(y-y_fit).^2); %Bevington and Robinson, Eq. 8.4
p = chi2cdf(chi2,nu,'upper');