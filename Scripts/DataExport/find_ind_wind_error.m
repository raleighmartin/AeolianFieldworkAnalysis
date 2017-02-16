%% script for determining which indices in wind don't have errors
% u = raw streamwise wind speed
% v = raw lateral wind speed
% w = raw vertical wind speed
% diag = diagnostic flag (0 if no error)
% u_sigma_max = maximum # of std devs in total wind away from mean wind for error detection
% ind_wind_error = indices of wind speed time series with errors
% ind_wind_noerror = indices of wind speed time series with no errors

function [ind_wind_error, ind_wind_noerror] = find_ind_wind_error(u, v, w, diag, u_sigma_max)

%calculate total wind
u_total = sqrt(u.^2+v.^2+w.^2);

%calculate maximum total wind based on multiple of std dev
u_total_max = mean(u_total)+u_sigma_max*std(u_total);

%find error points
ind_diag = find(diag~=0); %find instances of nonzero diagnostic flag
ind_u_total = find(u_total > u_total_max); %find instances of u_total > u_total_max

%get indices of error and non-error points
ind_wind_error = union(ind_diag, ind_u_total); %indices of error points
[~, ind_wind_noerror] = setdiff((1:length(u))',ind_wind_error); %indices of non-error points