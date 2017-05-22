%% PURPOSE
% function to determine flux midpoint height for each trap assuming power
% law flux profile and inputs below

%% INPUTS
% z_bottom_profile: profile of heights of bottoms of BSNEs
% z_trapheight_profile: profile of trap heights of BSNEs
% kz: assumed exponenent for power law fit

%% OUTPUTS
% z_profile = resulting values for flux midpoint heights 

function z_profile = z_profile_calc_powerlaw(z_bottom_profile,z_trapheight_profile,kz)

%get number of heights, initialize z_profile
N_z = length(z_bottom_profile);
z_profile = zeros(N_z,1);

for i = 1:N_z
    z1 = z_bottom_profile(i);
    z2 = z_bottom_profile(i)+z_trapheight_profile(i);
    z_profile(i) = (1/((z2-z1)*(1-kz))*...
        (z2^(1-kz)-z1^(1-kz)))^(-1/kz); %new calc matching actual and expected z
end