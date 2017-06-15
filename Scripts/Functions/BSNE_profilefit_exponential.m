%% PURPOSE
% function to iteratively fit exponential to BSNE profile and determine
% characteristic heights of BSNEs from this profile

%%INPUT VALUES
% qz_profile = partial flux values (vector)
% z_bottom_profile = elevations of bottom of traps for partial flux values (vector)
% z_trapheight_profile = heights of trap openings (vector)
% sigma_qz_profile = uncertainties in qz values (vector)
% sigma_z_profile = uncertainties in elevations (vector)
% zq_estimated = estimated value of zq (saltation layer height) for uncertainty estimation (scalar)

%%OUTPUT VALUES
% z_profile = optimized trap heights (vector)
% q0 = calculated flux scaling parameter (scalar)
% zq = calculated saltation layer height (scalar)
% Q = calculated total flux (scalar)
% sigma_q0 = uncertainty in q0 (scalar)
% sigma_zq = uncertainty in zq (scalar)
% sigma_Q = uncertainty in Q (scalar)
% qz_pred = predicted values of partial flux for exponential fit (vector)
% sigma_qz_pred = uncertainty in qz_pred values (vector)
% sigma_logqz_pred = uncertainty in log(qz_pred) values (vector)
% sigma2_q0zq = covariance of q0 and zq uncertainties
% z_profile_geomean = z_profile vector, but calculated as geometric mean height of each trap
% q0_geomean = q0 vector, but calculated as geometric mean height of each trap
% zq_geomean = zq vector, but calculated as geometric mean height of each trap
% Q_geomean = Q vector, but calculated as geometric mean height of each trap

function [z_profile,q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_q0zq,...
    z_profile_geomean,q0_geomean,zq_geomean,Q_geomean] = ...
    BSNE_profilefit_exponential(qz_profile, z_bottom_profile, z_trapheight_profile, sigma_qz_profile, sigma_z_profile, zq_estimated)

%start with guess of BSNE heights as arithmetic mean of traps
z_profile = z_bottom_profile+z_trapheight_profile/2;

%exponential fit to flux profile
if nargin == 3 %if no uncertainties provided
    [~,zq] = qz_profilefit_exponential(qz_profile, z_profile);
else %if uncertainties provided
    [~,zq] = qz_profilefit_exponential(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile, zq_estimated);
end

%now that we have zq, redo calculation of BSNE heights
z_profile_old = z_profile; %document previous z-profile to see difference
z_profile = z_profile_calc_exponential(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights

%iterate until the z_profile_difference is minutely small
while(z_profile_difference>1e-8)
    if nargin == 3 %if no uncertainties provided
        [~,zq] = qz_profilefit_exponential(qz_profile, z_profile);
    else %if uncertainties provided
        [~,zq] = qz_profilefit_exponential(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile, zq_estimated);
    end
    z_profile_old = z_profile; %document previous z-profile to see difference
    z_profile = z_profile_calc_exponential(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
    z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
end

z_profile_geomean = sqrt(z_bottom_profile.*(z_bottom_profile+z_trapheight_profile)); %alternative z profile with BSNE geometric mean height
%final profile fit
if nargin == 3 %if no uncertainties provided
    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_q0zq] = qz_profilefit_exponential(qz_profile, z_profile); %height-integrated flux from exponential fit
    [q0_geomean,zq_geomean,Q_geomean] = qz_profilefit_exponential(qz_profile,z_profile_geomean); %same thing, but using geometric mean height
else %if uncertainties provided
    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_q0zq] = qz_profilefit_exponential(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile, zq_estimated); %height-integrated flux from exponential fit
    [q0_geomean,zq_geomean,Q_geomean] = qz_profilefit_exponential(qz_profile,z_profile_geomean, sigma_qz_profile, sigma_z_profile, zq_estimated); %same thing, but using geometric mean height
end