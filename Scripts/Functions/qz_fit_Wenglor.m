%% function to compute saltation flux profile fit for Wenglors
% specifically addresses issues of:
% 1. repeated Wenglor heights
% 2. zero values for Wenglors

%%INPUTS
%qW = Wenglor partial fluxes
%zW = Wenglor heights
%sigma_qW = Wenglor partial flux uncertainties
%sigma_zW = Wenglor height uncertainties
%%%%%Q_min = detection limit for Q (g/m/s), set qW zero if below this value
%%%%zq_Q_min = assumed saltation height (m) for detection limit for exponential profile for detection limit for individual Wenglor
%zW_min = minimum number of Wenglors for profile fitting
%zq_estimated = estimated value of zq for uncertainty propagation

%%OUTPUTS
%q0 = profile fit scaling value
%zq = profile fit saltation layer height
%Q = profile fit total flux
%sigma_q0 = profile fit scaling value uncertainty
%sigma_zq = profile fit saltation layer height uncertainty
%sigma_Q = profile fit total flux uncertainty
%Chi2_Qfit = Chi2 value for profile fit
%df_Qfit = degrees of freedom for profile fit
%qW_unique = Wenglor partial fluxes for unique heights
%zW_unique = Wenglor heights for unique heights
%sigma_qW_unique = Wenglor partial flux uncertainties for unique heights
%sigma_zW_unique = Wenglor height uncertainties for unique heights
%N_z_fit = number of Wenglor heights for fit

%function [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,Chi2_Qfit,df_Qfit,qW_unique,zW_unique,sigma_qW_unique,sigma_zW_unique,N_z_fit] = qz_fit_Wenglor(qW, zW, sigma_qW, sigma_zW, Q_min, zq_Q_min, zW_min, zq_estimated)
function [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,Chi2_Qfit,df_Qfit,qW_unique,zW_unique,sigma_qW_unique,sigma_zW_unique,N_z_fit] = qz_fit_Wenglor(qW, zW, sigma_qW, sigma_zW, zW_min, zq_estimated)

%deal with repeated values
zW_unique = unique(zW);
N_zW_unique = length(zW_unique);
qW_unique = zeros(size(zW_unique));
sigma_qW_unique = zeros(size(zW_unique));
sigma_zW_unique = zeros(size(zW_unique));
for s = 1:N_zW_unique
    ind_zW = find(zW==zW_unique(s));
    if length(ind_zW)>1 %compute mean and uncertainty for repeated heights
%         q_min = (Q_min/zW_unique(s))*exp(-zW_unique(s)/zq_Q_min); %get min q detection limit for this height
%         qbar_z = mean(qW(ind_zW)); %get mean q for this height
%         if qbar_z < q_min %if mean q at height is below detection limit
%             qbar_z = 0; sigma_qbar_z = 0; %then just set values to zero
%         else %otherwise, use script "MeanUncertainty.m"
            [qbar_z, sigma_qbar_z] = MeanUncertainty(qW(ind_zW), sigma_qW(ind_zW)); %compute flux mean and uncertainy for repeated heights
%         end
        sigma_zW_unique(s) = mean(sigma_zW(ind_zW)); %no reduction in uncertainty for height, because all uncertainties are correlated, so just take mean of values
        qW_unique(s) = qbar_z;
        sigma_qW_unique(s) = sigma_qbar_z;
    else %otherwise, if only Wenglor at given height, just use existing values
        qW_unique(s) = qW(ind_zW);
        sigma_qW_unique(s) = sigma_qW(ind_zW);
        sigma_zW_unique(s) = sigma_zW(ind_zW);
    end
end

%Remove 0 values for fitting
ind_fit = find(qW_unique>0);
N_z_fit = length(ind_fit);
qW_fit = qW_unique(ind_fit);
zW_fit = zW_unique(ind_fit);
sigma_qW_fit = sigma_qW_unique(ind_fit);
sigma_zW_fit = zeros(1,N_z_fit); %neglect uncertainty in Wenglor height, which is already accounted for by calibration

%Perform profile fit to get q0, zq, and Q if sufficient points for fitting
if N_z_fit>=zW_min
    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qW_pred,sigma_qW_pred] = qz_profilefit(qW_fit,zW_fit,sigma_qW_fit,sigma_zW_fit,zq_estimated);
    qW_residuals = qW_pred - qW_fit; %residuals between observed and predicted q
    Chi2_Qfit = sum((qW_residuals./sigma_qW_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
    df_Qfit = N_z_fit-2; %compute degrees of freedom for Qfit
else %otherwise, set to NaN
    q0=NaN;
    zq=NaN;
    Q=NaN;
    sigma_q0=NaN;
    sigma_zq=NaN;
    sigma_Q = NaN;
    Chi2_Qfit=NaN;
    df_Qfit=NaN;
end