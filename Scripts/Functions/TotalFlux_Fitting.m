%% function to compute total saltation flux by fitting method

function [Q, sigma_Q, Chi2_Qfit, df_Qfit, zq, sigma_zq] = TotalFlux_Fitting(zW,qbar,sigma_qbar)               

%deal with repeated values
zW_unique = unique(zW);
N_zW_unique = length(zW_unique);
qbar_unique = zeros(size(zW_unique));
sigma_qbar_unique = zeros(size(zW_unique));
sigma_zW_unique = zeros(size(zW_unique));
for s = 1:N_zW_unique
    ind_zW = find(zW==zW_unique(s));
    if length(ind_zW)>1 %compute mean and uncertainty for repeated heights
        q_min = (Q_min/zW_unique(s))*exp(-zW_unique(s)/zq_Q_min); %get min q detection limit for this height
        qbar_z = mean(qbar(ind_zW)); %get mean q for this height
        if qbar_z < q_min %if mean q at height is below detection limit
            qbar_z = 0; sigma_qbar_z = 0; %then just set values to zero
        else %otherwise, use script "MeanUncertainty.m"
            [qbar_z, sigma_qbar_z] = MeanUncertainty(qbar(ind_zW), sigma_qbar(ind_zW)); %compute flux mean and uncertainy for repeated heights
        end
        sigma_zW_unique(s) = mean(sigma_zW(ind_zW)); %no reduction in uncertainty for height, because all uncertainties are correlated, so just take mean of values
        qbar_unique(s) = qbar_z;
        sigma_qbar_unique(s) = sigma_qbar_z;
    else %otherwise, if only Wenglor at given height, just use existing values
        qbar_unique(s) = qbar(ind_zW);
        sigma_qbar_unique(s) = sigma_qbar(ind_zW);
        sigma_zW_unique(s) = sigma_zW(ind_zW);
    end
end

%Remove 0 values for fitting
ind_fit = find(qbar_unique>0);
N_fit = length(ind_fit);
qbar_fit = qbar_unique(ind_fit);
zW_fit = zW_unique(ind_fit);
sigma_qbar_fit = sigma_qbar_unique(ind_fit);
sigma_zW_fit = zeros(1,N_fit); %neglect uncertainty in Wenglor height, which is already accounted for by calibration

%Perform profile fit to get q0, zq, and Q if sufficient points for fitting
if N_fit>=zW_limit
    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,q_pred,sigma_q_pred] = qz_profilefit(qbar_fit,zW_fit,sigma_qbar_fit,sigma_zW_fit);
    q_residuals = q_pred - qbar_fit; %residuals between observed and predicted q
    Chi2_Qfit = sum((q_residuals./sigma_q_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
    df_Qfit = N_fit-2; %compute degrees of freedom for Qfit
else %otherwise, set to NaN
    Q=NaN;
    sigma_Q = NaN;
    zq=NaN;
    sigma_q0=NaN;
    sigma_zq=NaN;
    Chi2_Qfit=NaN;
    df_Qfit=NaN;
end