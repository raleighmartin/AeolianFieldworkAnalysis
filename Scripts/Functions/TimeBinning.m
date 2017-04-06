%% function to bin values by time

% OUTPUTS
% t_bins_mid = midpoint time in bin
% t_bins_lowerbound = earliest time in bin
% t_bins_upperbound = latest time in bin
% y_bin_avg = mean y-value in bin
% y_bin_SE = standard error y-value in bin
% y_bin_N = number of values in bin
% y_binned = all y-values in bin

% INPUTS
% t = times
% y = data values
% bin_delta_t = time increment for bins

function [t_bins_mid, t_bins_lowerbound, t_bins_upperbound, y_bin_avg, y_bin_SE, y_bin_N, y_binned] = TimeBinning(t,y,bin_delta_t)

%get maximum and minimum t
min_t = min(t);
max_t = max(t);

%get bin information - only allow bins that span full range of bin_delta_t
N_bins = floor(max_t/bin_delta_t); %determine number of bins
t_bins_lowerbound = min_t:bin_delta_t:min_t+(N_bins-1)*bin_delta_t; %lower time bound of bins
t_bins_upperbound = t_bins_lowerbound + bin_delta_t; %upper time bound of bins
t_bins_mid = (t_bins_upperbound+t_bins_lowerbound)/2; %get midpoint time in bin

%initialize binned values
y_binned = cell(N_bins,1); %list of all y's in bin
y_bin_avg = zeros(N_bins,1); %mean y in bin
y_bin_N = zeros(N_bins,1); %number of values in bin
y_bin_SE = zeros(N_bins,1); %standard error of y in bin

%group data into bins
for i = 1:N_bins
    ind_bin = find(t>=t_bins_lowerbound(i) & t<t_bins_upperbound(i)); %get indices of times in bin
    y_binned{i} = y(ind_bin); %get y-values in bin, add to cell array
    y_bin_N(i) = length(ind_bin); %get number of values in bin, add to vector
    y_bin_avg(i) = mean(y(ind_bin)); %get mean y-value in bin, add to vector
    y_bin_SE(i) = std(y(ind_bin))/sqrt(y_bin_N(i)); %get standard error y-value in bin, add to vector
end