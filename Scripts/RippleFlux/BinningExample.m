%% EXAMPLE SCRIPT FOR HOW TO BIN DATA INTO EQUAL TIME INTERVALS

%initialize
clearvars;

%time interval (s) for binning data
bin_delta_t = 1000;

%information for creating random values
min_t = 2000; %minimum t (s) for random values
max_t = 10500; %maximum t (s) for random values
N_rand = 100; %number of random times

%generate random times and values
t = sort(rand(N_rand,1)*(max_t-min_t)+min_t); %generate random times
y = rand(N_rand,1)+t/(max_t-min_t); %generate random values associated with times

%run time binning function - only return values needed for analysis
[t_bins_mid, ~, ~, y_bin_avg, y_bin_SE, ~, ~] = TimeBinning(t,y,bin_delta_t);

%plot all values then bin average values
figure(1); clf; hold on;
plot(t,y,'k.'); %plot all values
errorbar(t_bins_mid,y_bin_avg,y_bin_SE,'b+'); %plot bin-average values with error bars for standard deviation
legend('unbinned','binned','Location','NorthWest'); %create legend
xlabel('time (s)'); %x-label: times
ylabel('y'); %y-label: variable value