%function to iteratively fit exponential to BSNE profile and determine characteristic heights of BSNEs from this profile

function [z_profile,q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_q0zq] = ...
    BSNE_profilefit(qz_profile, z_bottom_profile, z_trapheight_profile, sigma_qz_profile, sigma_z_profile)

%start with guess of BSNE midpoint heights as arithmetic mean of traps
z_profile = z_bottom_profile+z_trapheight_profile/2;

%exponential fit to flux profile
[~,zq] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile);

%now that we have zq, redo calculation of BSNE heights
z_profile_old = z_profile; %document previous z-profile to see difference
z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights

%iterate until the z_profile is minutely small
while(z_profile_difference>1e-8)
    [~,zq] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile); %exponential fit to profile
    z_profile_old = z_profile; %document previous z-profile to see difference
    z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
    z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
end

%final profile fit
[q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_q0zq] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile); %height-integrated flux from exponential fit