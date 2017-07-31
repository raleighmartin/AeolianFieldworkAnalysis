%% THIS SCRIPT TAKES BED ELEVATION TIMESERIES AND COMPUTES RIPPLE AMPLITUDES AND MIGRATION VELOCITIES TO OBTAIN RIPPLE FLUX

%initialize
clearvars;

%define starttimes and endtimes for analysis
starttimes_analysis = [datetime(2015,5,15,11,30,0),datetime(2015,5,23,12,0,0),datetime(2015,5,24,13,0,0)];
endtimes_analysis = [datetime(2015,5,15,15,30,0),datetime(2015,5,23,16,0,0),datetime(2015,5,24,16,0,0)];

%determine number of time intervals for analysis
N_analysis_windows = length(starttimes_analysis);

%initialize cell arrays of values
crest_times = cell(N_analysis_windows,1);

%go through each time interval to do analysis
for i = 1:N_analysis_windows

   %get times
   starttime_window = starttimes_analysis(i);
   endtime_window = endtimes_analysis(i);
  
   %something to load in the appropriate data for this analysis
   
   %run functions to do necessary steps for this window
   %pre-processing (smoothing, etc.)
   
   %getting crests and troughs and their times
   vector_crest_times_example = [1000, 1523, 1598, 21525]*i; 
      
   %add values to cell array
   crest_times{i} = vector_crest_times_example;
   
   %calculating values
end