%% FUNCTION TO GROUP VALUES TOGETHER INTO BINS
% This function allows bins to have a range of possible sizes

%% input values
% y = raw values for binning
% y_bin_minrange = minimum y range for grouping raw values together in one bin
% y_bin_maxrange = minimum y range for grouping raw values together in one bin
% bin_N_min = minimum number of values per bin

%% output values
% y_bin_values = values grouped by bin
% y_bin_N = number of values in each bin
% y_bin_min = minimum value in each bin
% y_bin_max = maximum value in each bin
% y_bin_avg = mean value in each bin
% y_bin_SE = standard error of values in each bin

function [y_bin_values, y_bin_N, y_bin_min, y_bin_max, y_bin_avg, y_bin_SE] = ...
    Binning(y, y_bin_minrange, y_bin_maxrange, bin_N_min)

y = sort(y); %sort y
N_y = length(y); %get number of y
y_bin_values = cell(1); %initialize list of y_bin_values
ind_bin = 1; %initialize index of bin
y_bin_list = []; %initialize list of y in bin
for i = 1:N_y
    if isempty(y_bin_list) %if no y's in bin, create first value for bin
        y_bin_list = y(i); %add y value to list
    elseif (y(i)-min(y_bin_list))<=y_bin_minrange; %if value of y is within min range, add to list
        y_bin_list = [y_bin_list; y(i)]; %add y value to list  
    elseif (y(i)-min(y_bin_list))>y_bin_maxrange; %if value of y is outside max range, create new bin
        y_bin_values{ind_bin,1} = y_bin_list; %add current y list to cell array
        ind_bin = ind_bin+1; %increment to next y bin
        y_bin_list = y(i); %add y value to new list
    elseif length(y_bin_list)<bin_N_min; %if next y is within min-max range, but have not reached minimum number of entries, add to list
        y_bin_list = [y_bin_list; y(i)]; %add y value to list
    elseif (N_y-i)<bin_N_min; %if new bin could be created, but there are few entries left at the end of the list, add to current list
        y_bin_list = [y_bin_list; y(i)]; %add y value to list
    else %create new bin
        y_bin_values{ind_bin,1} = y_bin_list; %add current y list to cell array
        ind_bin = ind_bin+1; %increment to next ybin
        y_bin_list = y(i); %add y value to new list
    end
    
    if i == N_y %if at the end of y's, add current y list to cell array
        y_bin_values{ind_bin,1} = y_bin_list; %add y value to list
    end
end

y_bin_N = cellfun(@length,y_bin_values);
y_bin_min = cellfun(@min,y_bin_values);
y_bin_max = cellfun(@max,y_bin_values);
y_bin_avg = cellfun(@mean,y_bin_values);
y_bin_std = cellfun(@std,y_bin_values);
y_bin_SE = y_bin_std./sqrt(y_bin_N);