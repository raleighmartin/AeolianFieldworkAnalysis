%function to iteratively fit power law to BSNE profile and determine characteristic heights of BSNEs from this profile

function [z_profile,qp,kz,sigma_qp,sigma_kz,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_qpkz,...
    z_profile_geomean,qp_geomean,kz_geomean] = ...
    BSNE_profilefit_powerlaw(qz_profile, z_bottom_profile, z_trapheight_profile, sigma_qz_profile, sigma_z_profile)

%start with guess of BSNE heights as arithmetic mean of traps
z_profile = z_bottom_profile+z_trapheight_profile/2;

%power law fit to flux profile
if nargin == 3 %if no uncertainties provided
    [~,kz] = qz_profilefit_powerlaw(qz_profile, z_profile);
else %if uncertainties provided
    [~,kz] = qz_profilefit_powerlaw(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile);
end

%now that we have kz, redo calculation of BSNE heights
z_profile_old = z_profile; %document previous z-profile to see difference
z_profile = z_profile_calc_powerlaw(z_bottom_profile,z_trapheight_profile,kz); %calculate new BSNE midpoint heights based on zq
z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights

%iterate until the z_profile_difference is minutely small
while(z_profile_difference>1e-8)
    if nargin == 3 %if no uncertainties provided
        [~,kz] = qz_profilefit_powerlaw(qz_profile, z_profile);
    else %if uncertainties provided
        [~,kz] = qz_profilefit_powerlaw(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile);
    end
    z_profile_old = z_profile; %document previous z-profile to see difference
    z_profile = z_profile_calc_powerlaw(z_bottom_profile,z_trapheight_profile,kz); %calculate new BSNE midpoint heights based on zq
    z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
end

z_profile_geomean = sqrt(z_bottom_profile.*(z_bottom_profile+z_trapheight_profile)); %alternative z profile with BSNE geometric mean height
%final profile fit
if nargin == 3 %if no uncertainties provided
    [qp,kz,sigma_qp,sigma_kz,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_qpkz] = qz_profilefit_powerlaw(qz_profile, z_profile); %height-integrated flux from exponential fit
    [qp_geomean,kz_geomean] = qz_profilefit_powerlaw(qz_profile,z_profile_geomean); %same thing, but using geometric mean height
else %if uncertainties provided
    [qp,kz,sigma_qp,sigma_kz,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_qpkz] = qz_profilefit_powerlaw(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile); %height-integrated flux from exponential fit
    [qp_geomean,kz_geomean] = qz_profilefit_powerlaw(qz_profile,z_profile_geomean,sigma_qz_profile,sigma_z_profile); %same thing, but using geometric mean height
end