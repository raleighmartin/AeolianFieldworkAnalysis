%% function to compute saltation flux profile summation for Wenglors

%%INPUTS
%qW = Wenglor partial fluxes
%zW = Wenglor heights (m)
%sigma_qW = uncertainty in Wenglor partial fluxes
%zq_BSNE = saltation layer height from BSNE (m)
%sigma_zq_BSNE = uncertainty on saltation layer height from BSNE

%%OUTPUTS
%Qsum = total flux from summation
%sigma_Qsum = uncertainty in total flux calculation

function [Qsum,sigma_Qsum] = qz_summation_Wenglor(qW, zW, sigma_qW, zq_BSNE, sigma_zq_BSNE)

ind_summation = find(zW>0);


%compute delta_z weights for additive flux calculations
z1_Qsum = [0 sqrt(zW(1:end-1).*zW(2:end))]; %bottom of summation bins
z2_Qsum = [sqrt(zW(1:end-1).*zW(2:end)) Inf];
q0_Qsum = qW./(exp(-zW/zq_BSNE));
deltaQ = q0_Qsum.*zq_BSNE.*(exp(-z1_Qsum/zq_BSNE)-exp(-z2_Qsum/zq_BSNE));
sigma_deltaQ = deltaQ.*(sigma_qW./qW);
Qsum = sum(deltaQ);
sigma_Qsum = sqrt(sum(sigma_deltaQ(~isnan(sigma_deltaQ)).^2)+(Qsum*sigma_zq_BSNE/zq_BSNE)^2);