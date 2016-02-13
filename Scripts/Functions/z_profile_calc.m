%function to determine flux midpoint height assuming exponential flux
%profile and:
%1. z_bottom_profile: profile of heights of bottoms of BSNEs
%2. z_trapheight_profile: profile of trap heights of BSNEs
%3. zq: assumed e-folding height of profile

function z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq)

%get number of heights, initializq z_profile
N_z = length(z_bottom_profile);
z_profile = zeros(N_z,1);

for i = 1:N_z;
    z_profile(i) = -zq*log(0.5*(...
        exp(-(z_bottom_profile(i)+z_trapheight_profile(i))/zq)+...
        exp(-(z_bottom_profile(i)/zq))));
end
    