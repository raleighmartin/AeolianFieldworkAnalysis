%% function to compute saltation flux profile assuming exponential

%%INPUTS
%z - observation heights (m)
%qz - vertical fluxes (g/m^2/s)
%sigma_z = height uncertainty
%sigma_qz = flux uncertainty

%%OUTPUTS
%q0 - profile fit (g/m^2/s)
%zq - profile fit e-folding height (m)
%Q - profile fit total flux (g/m/s)
%qz_fit - predicted q from profile fit
%sigma_q0 - uncertainty in q0
%sigma_Q - uncertainty in Q
%sigma_zq - uncertainty in zq
%sigma_qz_fit - uncertainty in predictions of qz (varies with qz)
%sigma_logqz_fit - uncertainty in predictions of log(qz)
%sigma2_q0zq - covariance in q0 and zq

%% FITTING TO:
% q = q0*exp(-z/zq)
% logq = a+b*z
% q0 = exp(a)
% zq = -1/b

%use method of bevington and robinson (p. 105)
function [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_fit,sigma_qz_fit,sigma_logqz_fit,sigma2_q0zq] = qz_profilefit(qz, z, sigma_qz, sigma_z)

%set input error values arbitrarily as one if input is not given
if nargin == 2
    sigma_qz = ones(size(qz)); 
    sigma_z = ones(size(z));
end

%compute log of flux and associated error
logqz = log(qz); %[log(q)]
sigma_logqz = sqrt((sigma_qz.^2./qz.^2)+(sigma_z.^2./z.^2)); %[log(q)] (combines height and flux error)

%perform linear fit and estimate parameters from fit
[a, b, sigma_a, sigma_b, ~, sigma_logqz_fit, sigma2_ab] = linearfit(z,logqz,sigma_logqz);
q0 = exp(a); %[g/m^2/s]
zq = -1/b; %[m]
sigma_q0 = sigma_a*q0; %[g/m^2/s]
sigma_zq = sigma_b*zq.^2; %[m]
qz_fit = q0*exp(-z/zq); %prediction of qz [g/m^2/s] from least squares fit
sigma_qz_fit = sigma_logqz_fit.*qz_fit; %confidence on prediction of qz [g/m^2/s] from least squares fit

%calculate Q and sigma_Q
Q = q0*zq; %get total flux [g/m/s]
sigma2_q0zq = (exp(a)/b^2)*sigma2_ab;
sigma_Q = sqrt((sigma_q0*zq)^2 + (sigma_zq*q0)^2 + 2*sigma2_q0zq*Q); %estimate uncertainty in total flux including covariance of parameters (Bevington and Robinson, Eq. 3.13)