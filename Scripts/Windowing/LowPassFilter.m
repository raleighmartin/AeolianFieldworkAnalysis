%y_in is input timeseries
%dt is timestep of input timeseries (s)
%fc is cutoff frequency of low-pass filter (Hz)
%y_out is output timeseries

function [y_out] = LowPassFilter(y_in,dt,fc)

Wn = 2*pi*dt*fc; %get the cutoff frequency in radians/s
[b,a] = butter(6, Wn); %create 6th order Butterworth filter
y_out = filtfilt(b,a,y_in); %apply zero phase shift filtering to data