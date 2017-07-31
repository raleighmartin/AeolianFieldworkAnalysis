%% compute moving average of timeseries y with timestamp t and averaging time T
% use MATLAB "filter" function, but correct for delay and bias at beginning
% and end of timeseries introduced by this method
% Dependencies: NONE
% Used by: interpolate_timeseries, WindFluxPlot

function [y_ma] = moving_average(y,t,T)

%compute dt based on most common increment in timeseries
dt = mode(diff(t));

%compute window size based on averaging time, T
windowSize = round(T/dt);

%compute filter shape based on window size
b = (1/windowSize)*ones(1,windowSize);

%apply filter to timeseries
y_filter = filter(b,1,y);

%initialize moving average
y_ma = zeros(size(y));

%number of timesteps to shift filtered timeseries to eliminate delay
shiftIndex = round(windowSize/2);
y_ma(1:shiftIndex) = y_filter(windowSize); %up to shift_index, use first full filtered value
y_ma((windowSize-shiftIndex):(end-shiftIndex)) = y_filter(windowSize:end); %starting from shift_index, shift filtered values
y_ma((end-shiftIndex):end) = y_filter(end); %for last few points, use the last filtered value

end