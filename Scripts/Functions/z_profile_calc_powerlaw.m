%function to determine flux midpoint height assuming power law
%profile and:
%1. z_bottom_profile: profile of heights of bottoms of BSNEs
%2. z_trapheight_profile: profile of trap heights of BSNEs
%3. kz: assumed exponenent for power law fit

function z_profile = z_profile_calc_powerlaw(z_bottom_profile,z_trapheight_profile,kz)

%get number of heights, initialize z_profile
N_z = length(z_bottom_profile);
z_profile = zeros(N_z,1);

for i = 1:N_z
    z1 = z_bottom_profile(i);
    z2 = z_bottom_profile(i)+z_trapheight_profile(i);
    z_profile(i) = (1/((z2-z1)*(1-kz))*...
        (z2^(1-kz)-z1^(1-kz)))^(-1/kz); %new calc matching actual and expected z
%     figure(1); clf;
%     z_plot = linspace(z1,z2,100);
%     q_plot = z_plot.^(-kz);
%     plot(z_plot,q_plot);
end