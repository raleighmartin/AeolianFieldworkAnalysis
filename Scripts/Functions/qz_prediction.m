%function to compute expected flux and associated uncertainty at a given
%height based on profile fit paramaters and their uncertainties

function [qz_pred, sigma_qz_pred] = qz_prediction(z, q0, zq, sigma_q0, sigma_zq, sigma2_q0zq)

%compute expected flux at given height
qz_pred = q0.*exp(-z/zq); %(g/m^2/s)

%compute uncertainty in expected flux at Wenglor height
sigma_qz_pred = ...
    qz_pred.*sqrt((sigma_q0/q0).^2 + ...
    (sigma_zq.*z/zq^2).^2 + ...
    (2*sigma2_q0zq.*z/(q0*zq^2)));
