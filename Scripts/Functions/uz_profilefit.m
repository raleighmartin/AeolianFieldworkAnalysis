%% function to compute shear velocity and roughness length from wind profile

%%INPUTS
%z - observation heights (m)
%u - wind speed (m/s)
%sigma_u - wind uncertainty

%%OUTPUTS
%ust - shear velocity (m/s)
%zs - profile fit roughness height (m)
%sigma_ust - uncertainty in ust
%sigma_zs - uncertainty in zs

%% FITTING TO:
% u = ust/kappa*log(z/zs)
% u = a+b*log(z)
% ust = b*kappa
% zs = exp(-a*kappa/ust)

%use method of bevington and robinson for linear fitting and error estimation (p. 105)
function [ust,zs,sigma_ust,sigma_zs] = uz_profilefit(u, z, sigma_u)

%set kappa = 0.4
kappa = 0.4;

%compute log of height
logz = log(z);

%perform linear fit and estimate parameters from fit
if nargin == 2 %if no uncertainties given, use function linearfit without sigma_logqz
    [a, b, sigma_a, sigma_b] = linearfit(logz,u);
else %otherwise, use function with sigma_u
    [a, b, sigma_a, sigma_b] = linearfit(logz,u,sigma_u);
end
ust = b*kappa; %[m/s]
zs = exp(-a*kappa/ust); %[m]
sigma_ust = sigma_b*kappa; %[m/s]
sigma_zs = (kappa*zs/ust)*sqrt(sigma_a^2 + sigma_ust^2*a^2/ust^2); %[m]