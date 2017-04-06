%% function to compute saltation flux profile assuming exponential

%%INPUTS
%z - observation heights (m)
%qz - vertical fluxes (g/m^2/s)
%sigma_z = height uncertainty
%sigma_qz = flux uncertainty

%%OUTPUTS
%qp - profile fit (g/m^2/s)
%kz - profile fit power
%Q - profile fit total flux (g/m/s)
%qz_fit - predicted q from profile fit
%sigma_qp - uncertainty in qp
%sigma_Q - uncertainty in Q
%sigma_kz - uncertainty in kz
%sigma_qz_fit - uncertainty in predictions of qz (varies with qz)
%sigma_logqz_fit - uncertainty in predictions of log(qz)
%sigma2_qpkz - covariance in q0 and zq

%% FITTING TO:
% q = qp*z.^-kz
% logq = a+b*logz
% qp = exp(a)
% kz = -b

%use method of bevington and robinson (p. 105)
function [qp,kz,sigma_qp,sigma_kz,qz_fit,sigma_qz_fit,sigma_logqz_fit,sigma2_qpkz] = qz_profilefit_powerlaw(qz, z, sigma_qz, sigma_z)

% %set input error values arbitrarily as one if input is not given
% if nargin == 2
%     sigma_qz = ones(size(qz)); 
%     sigma_z = ones(size(z));
% end

%compute log of flux
logqz = log(qz); %[log(q)]

%compute log of height
logz = log(z); %[log(z)]

%perform linear fit and estimate parameters from fit
if nargin == 2 %if no uncertainties given, use function linearfit without sigma_logqz
    [a, b, sigma_a, sigma_b, ~, sigma_logqz_fit, sigma2_ab, da_dy, db_dy] = linearfit(logz,logqz);
else %otherwise, use function with sigma_logqz
    sigma_logqz = sqrt((sigma_qz.^2./qz.^2)+(sigma_z.^2./z.^2)); %[log(q)] compute uncertainty in log(qz)
    [a, b, sigma_a, sigma_b, ~, sigma_logqz_fit, sigma2_ab, da_dy, db_dy] = linearfit(logz,logqz,sigma_logqz);
end
qp = exp(a); %[g/m^2/s]
kz = -b;
sigma_qp = sigma_a*qp; %[g/m^2/s]
sigma_kz = sigma_b;
qz_fit = qp*z.^-kz; %prediction of qz [g/m^2/s] from least squares fit
sigma_qz_fit = sigma_logqz_fit.*qz_fit; %confidence on prediction of qz [g/m^2/s] from least squares fit

%calculate Q and sigma_Q
% Q = qp*kz; %get total flux [g/m/s]
sigma2_qpkz = sigma2_ab*qp;
% %sigma2_q0zq = (exp(a)/b^2)*sum((sigma_qz.^2/qz.^2).*da_dy.*db_dy); %%technically more correct, produces slightly larger covariance by accounting for z
% sigma_Q = sqrt((sigma_qp*kz)^2 + (sigma_kz*qp)^2 + 2*sigma2_qpkz*Q); %estimate uncertainty in total flux including covariance of parameters (Bevington and Robinson, Eq. 3.13)