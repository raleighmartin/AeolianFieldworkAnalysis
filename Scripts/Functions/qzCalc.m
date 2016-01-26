%% function, given BSNE mass, height of opening, width of opening, and duration of observation
% to determine the mass flux
% Dependencies: NONE
% Used by: ProcessBSNEs

function [qz_gm2s, sigma_qz_gm2s] = qzCalc(MassBSNE_g, HeightBSNE_cm, WidthBSNE_cm, DurationBSNE, sigma_MassBSNE_g, sigma_HeightBSNE_cm, sigma_WidthBSNE_cm, sigma_DurationBSNE)

%compute mean value
qz_gm2s = (MassBSNE_g)/(1e-4*HeightBSNE_cm*WidthBSNE_cm*seconds(DurationBSNE));

%compute uncertainty
sigma_qz_gm2s = qz_gm2s*sqrt(...
    (sigma_MassBSNE_g/MassBSNE_g).^2+...
    (sigma_HeightBSNE_cm/HeightBSNE_cm).^2+...
    (sigma_WidthBSNE_cm/WidthBSNE_cm).^2+...
    (sigma_DurationBSNE/DurationBSNE).^2);