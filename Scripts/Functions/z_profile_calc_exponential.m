%% PURPOSE
% function to determine flux midpoint height for each trap assuming exponential flux
% profile and inputs below

%% INPUTS
% z_bottom_profile: profile of heights of bottoms of BSNEs
% z_trapheight_profile: profile of trap heights of BSNEs
% zq: assumed e-folding height of profile

%% OUTPUTS
% z_profile = resulting values for flux midpoint heights 

function z_profile = z_profile_calc_exponential(z_bottom_profile,z_trapheight_profile,zq)

%get number of heights, initialize z_profile
N_z = length(z_bottom_profile);
z_profile = zeros(N_z,1);

for i = 1:N_z
    z_profile(i) = -zq*log((-zq/z_trapheight_profile(i))*... %new calc matching actual and expected z
        exp(-z_bottom_profile(i)/zq)*(exp(-z_trapheight_profile(i)/zq)-1));
end