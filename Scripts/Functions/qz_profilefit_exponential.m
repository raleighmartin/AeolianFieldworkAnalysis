%% function to compute saltation flux profile assuming exponential

%%INPUTS
%qz - vertical fluxes (g/m^2/s)
%z - observation heights (m)
%sigma_qz = flux uncertainty
%sigma_z = height uncertainty
%zq_estimated = value of zq for uncertainty estimation - set to 0.1 m if no value given

%%OUTPUTS
%q0 - profile fit (g/m^2/s)
%zq - profile fit e-folding height (m)
%Q - profile fit total flux (g/m/s)
%sigma_q0 - uncertainty in q0
%sigma_zq - uncertainty in zq
%sigma_Q - uncertainty in Q
%qz_fit - predicted q from profile fit
%sigma_qz_fit - uncertainty in predictions of qz (varies with qz)
%sigma_logqz_fit - uncertainty in predictions of log(qz)
%sigma2_q0zq - covariance in q0 and zq

%% FITTING TO:
% q = q0*exp(-z/zq)
% logq = a+b*z
% q0 = exp(a)
% zq = -1/b

%use method of bevington and robinson (p. 105)
function [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_fit,sigma_qz_fit,sigma_logqz_fit,sigma2_q0zq] = qz_profilefit_exponential(qz, z, sigma_qz, sigma_z, zq_estimated)

%set zq_estimated = 0.1 m if no value given
if nargin <= 4
    zq_estimated = 0.1; 
end

%compute log of flux
logqz = log(qz); %[log(q)]

%perform linear fit and estimate parameters from fit
if nargin == 2 %if no uncertainties given, use function linearfit without sigma_logqz
    [a, b, sigma_a, sigma_b, ~, sigma_logqz_fit, sigma2_ab, da_dy, db_dy] = linearfit(z,logqz);
else %otherwise, use function with sigma_logqz
    sigma_logqz = sqrt((sigma_qz.^2./qz.^2)+(sigma_z.^2./zq_estimated.^2)); %[log(q)] compute uncertainty in log(qz)
    %sigma_logqz = sqrt((sigma_qz.^2./qz.^2)+(sigma_z.^2./z.^2)); %[log(q)] compute uncertainty in log(qz)
    [a, b, sigma_a, sigma_b, ~, sigma_logqz_fit, sigma2_ab, da_dy, db_dy] = linearfit(z,logqz,sigma_logqz);
end
q0 = exp(a); %[g/m^2/s]
zq = -1/b; %[m]
sigma_q0 = sigma_a*q0; %[g/m^2/s]
sigma_zq = sigma_b*zq.^2; %[m]
qz_fit = q0*exp(-z/zq); %prediction of qz [g/m^2/s] from least squares fit
sigma_qz_fit = sigma_logqz_fit.*qz_fit; %confidence on prediction of qz [g/m^2/s] from least squares fit

%calculate Q and sigma_Q
Q = q0*zq; %get total flux [g/m/s]
sigma2_q0zq = q0*zq^2*sigma2_ab;
sigma_Q = sqrt((sigma_q0*zq)^2 + (sigma_zq*q0)^2 + 2*sigma2_q0zq*Q); %estimate uncertainty in total flux including covariance of parameters (Bevington and Robinson, Eq. 3.13)