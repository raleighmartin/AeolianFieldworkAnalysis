

z_Namikas03 = cell(N_Namikas03,1); %trap height (m)
q_Namikas03 = cell(N_Namikas03,1); %flux (g/m^2/s)

%% Iteratively compute BSNE heights, flux profiles, and uncertainties for each interval
for i = 1:N_Intervals
    %get flux and height profiles for fitting    
    qz_profile = FluxBSNE(i).qz.qz;
    z_bottom_profile = FluxBSNE(i).z.bottom; %bottom heights of BSNEs
    z_trapheight_profile = FluxBSNE(i).z.trapheight; %trap heights of BSNEs
    sigma_qz_profile = FluxBSNE(i).qz.sigma;
    sigma_z_profile = FluxBSNE(i).z.sigma_bottom;

    %start with guess of BSNE midpoint heights as arithmetic mean of traps
    z_profile = z_bottom_profile+z_trapheight_profile/2;
    
    %height-integrated flux from exponential fit
    [q0,zq,sigma_q0,sigma_zq,qz_pred,sigma_qz_pred,sigma_logqz_pred] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile);
    
    %now that we have zq, redo calculation of BSNE heights
    z_profile_old = z_profile; %document previous z-profile to see difference
    z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
    z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
    
    %iterate until the z_profile is minutely small
    while(z_profile_difference>1e-8)
        [q0,zq,sigma_q0,sigma_zq,qz_pred,sigma_qz_pred,sigma_logqz_pred] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile); %height-integrated flux from exponential fit
        z_profile_old = z_profile; %document previous z-profile to see difference
        z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
        z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
    end
    
    %% add to structured array
    FluxBSNE(i).qz.q0 = q0;
    FluxBSNE(i).qz.qz_pred = qz_pred;
    FluxBSNE(i).qz.sigma_q0 = sigma_q0;
    FluxBSNE(i).qz.sigma_qz_pred = sigma_qz_pred;
    FluxBSNE(i).qz.sigma_logqz_pred = sigma_logqz_pred;
    FluxBSNE(i).z.zq = zq;
    FluxBSNE(i).z.sigma_zq = sigma_zq;
    FluxBSNE(i).z.z = z_profile;
    FluxBSNE(i).z.sigma_z = sigma_z_profile;
end