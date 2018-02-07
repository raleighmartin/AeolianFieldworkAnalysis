%% function to compute total saltation flux by summation method

function [Qsum, sigma_Qsum] = TotalFlux_Summation(zW,qbar,sigma_qbar)               

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

%Perform summation to get Q
%compute delta_z weights for additive flux calculations
z1_Qsum = [0 sqrt(zW_unique(1:end-1).*zW_unique(2:end))]; %bottom of summation bins
z2_Qsum = [sqrt(zW_unique(1:end-1).*zW_unique(2:end)) Inf];
q0_Qsum = qbar_unique./(exp(-zW_unique/zq_BSNE));
sigma_q0_Qsum = sigma_qbar_unique./(exp(-zW_unique/zq_BSNE));
deltaQ = q0_Qsum.*zq_BSNE.*(exp(-z1_Qsum/zq_BSNE)-exp(-z2_Qsum/zq_BSNE));
sigma_deltaQ = sigma_q0_Qsum.*zq_BSNE.*(exp(-z1_Qsum/zq_BSNE)-exp(-z2_Qsum/zq_BSNE));
Qsum = sum(deltaQ);
sigma_Qsum = sqrt(sum(sigma_deltaQ.^2));