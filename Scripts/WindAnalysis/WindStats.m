%% ANALYZE STATISTICS OF WIND FLUCTUATIONS
clearvars; %initialize

%% choose values of dt for averaging
dt_min_sec = 0.02; %mininum dt, in seconds
dt_max_sec = 1000; %maximum dt, in seconds
N_dt_windowavg = 50; %intended number of dt's
dt_windowavg_sec = unique(round(logspace(0,log10(dt_max_sec/dt_min_sec),N_dt_windowavg)'))*dt_min_sec; %generate log-spaced dt's, round to nearest values, eliminate repeats
dt_windowavg = duration(0,0,dt_windowavg_sec); %convert to durations
N_dt_windowavg = length(dt_windowavg_sec); %get actual number of dt's

%% choose values of u for threshold crossing analysis
u_thr = 8:11; %m/s
N_thr = length(u_thr);

%% information about site for analysis
Site = 'Oceano';
Date = datetime(2015,5,18);
StartTime = datetime(2015,5,18,13,47,0);
EndTime = datetime(2015,5,18,17,47,0);
InstrumentType = 'Sonic';
Instrument = 'S1';

%% information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/WindStats/'; %folder for plots
SaveData_Path = strcat(folder_ProcessedData,'WindStats');
 
%% load processed data for site
ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Site);
Data = load(ProcessedData_Path); %load processed data

%% extract data for interval
[~, ~, IntervalN, IntervalInd] = ...
    ExtractVariableTimeInterval(ProcessedData.(InstrumentType).(Instrument),...
    StartTime,EndTime,'t','int','int'); %get indices of time interval
t = ProcessedData.(InstrumentType).(Instrument)(IntervalN).t.int(IntervalInd{1}); %get times
u_raw = ProcessedData.(InstrumentType).(Instrument)(IntervalN).u.int(IntervalInd{1}); %get raw u-values
v_raw = ProcessedData.(InstrumentType).(Instrument)(IntervalN).v.int(IntervalInd{1}); %get raw v-values
w_raw = ProcessedData.(InstrumentType).(Instrument)(IntervalN).w.int(IntervalInd{1}); %get raw w-values
[u_rot, v_rot, w_rot] = reorient_anemometers_vanboxel2004(u_raw, v_raw, w_raw); %rotate instrument

%% get timeseries with 10-second averages
dt_10s = duration(0,0,10);
[u_10s, t_10s] = window_average(u_rot, t, dt_10s);
[v_10s, ~] = window_average(v_rot, t, dt_10s);
[w_10s, ~] = window_average(w_rot, t, dt_10s);
theta_10s = atan(v_10s./u_10s)*180/pi;

%% go through each dt_avg, generate averaged timeseries, compute variances for different window averages
t_windowavg = cell(N_dt_windowavg,1);
u_windowavg = cell(N_dt_windowavg,1);
v_windowavg = cell(N_dt_windowavg,1);
w_windowavg = cell(N_dt_windowavg,1);
u_windowavg_var = zeros(N_dt_windowavg,1);
v_windowavg_var = zeros(N_dt_windowavg,1);
w_windowavg_var = zeros(N_dt_windowavg,1);
for i = 1:N_dt_windowavg
    [u_windowavg{i}, t_windowavg{i}] = window_average(u_rot, t, dt_windowavg(i));
    [v_windowavg{i}, ~] = window_average(v_rot, t, dt_windowavg(i));
    [w_windowavg{i}, ~] = window_average(w_rot, t, dt_windowavg(i));
    u_windowavg_var(i) = var(u_windowavg{i});
    v_windowavg_var(i) = var(v_windowavg{i});
    w_windowavg_var(i) = var(w_windowavg{i});
end

%% compute cumulative normalized variance (Eq. 6 in Furbish Time Averaging notes)
U = cumsum((u_windowavg_var(1:end)/u_windowavg_var(1)).*diff([0; dt_windowavg_sec]));
V = cumsum((v_windowavg_var(1:end)/v_windowavg_var(1)).*diff([0; dt_windowavg_sec]));
W = cumsum((w_windowavg_var(1:end)/w_windowavg_var(1)).*diff([0; dt_windowavg_sec]));

%% find excursion times for various thresholds

%get threshold up-crossing times for each threshold
u_thr_crossingtimes = cell(N_thr,1);
for i = 2:length(u_rot)
    for j = 1:N_thr
        if (u_rot(i)>=u_thr(j))&&(u_rot(i-1)<u_thr(j));
            u_thr_crossingtimes{j} = [u_thr_crossingtimes{j} t(i)];
        end
    end
end

%compute excursion durations based on up-crossing times
u_thr_excursiondurations = cell(N_thr,1);
u_thr_N_excursions = zeros(N_thr,1);
for j = 1:N_thr
    u_thr_excursiondurations{j} = sort(seconds(diff(u_thr_crossingtimes{j})));
    u_thr_N_excursions(j) = length(u_thr_excursiondurations{j});
end

%get autocorrelation of wind timeseries
[R_u, lags_u] = autocorrelation(u_rot,0.02,10000);
[R_v, lags_v] = autocorrelation(v_rot,0.02,10000);
[R_w, lags_w] = autocorrelation(w_rot,0.02,10000);

%% plot timeseries - using 10-second averages
figure(1); clf;
Markers = {'b-','r-','g-','k-'};
subplot(4,1,1);
plot(t_10s,u_10s,Markers{1},'LineWidth',2);
title([Site,', ',datestr(Date)]);
ylabel('u_{10s} (m/s)');
set(gca,'FontSize',16);
subplot(4,1,2);
plot(t_10s,v_10s,Markers{2},'LineWidth',2);
ylabel('v_{10s} (m/s)');
set(gca,'FontSize',16);
subplot(4,1,3);
plot(t_10s,w_10s,Markers{3},'LineWidth',2);
ylabel('w_{10s} (m/s)');
set(gca,'FontSize',16);
subplot(4,1,4);
plot(t_10s,theta_10s,Markers{4},'LineWidth',2);
ylabel('\theta_{10s} (\circ)');
set(gca,'FontSize',16);
print([folder_Plots,'WindTimeseries_10s_',Site,'_',datestr(Date),'.png'],'-dpng');

%% plot variance versus averaging time
figure(2); clf; hold on;
Markers = {'b-x','r-o','g-v'};
plot(dt_windowavg_sec,u_windowavg_var,Markers{1},'LineWidth',2);
plot(dt_windowavg_sec,v_windowavg_var,Markers{2},'LineWidth',2);
plot(dt_windowavg_sec,w_windowavg_var,Markers{3},'LineWidth',2);
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('T_{avg} (s)');
ylabel('Var(u_{i}(T)) (m^2/s^2)');
legend('u','v','w','Location','NorthEast');
title([Site,', ',datestr(Date)]);
set(gca,'FontSize',16);
print([folder_Plots,'WindVariance_Averaging_',Site,'_',datestr(Date),'.png'],'-dpng');

%% plot integrated normalized variance versus averaging time
figure(3); clf; hold on;
Markers = {'b-x','r-o','g-v'};
plot(dt_windowavg_sec,U,Markers{1},'LineWidth',2);
plot(dt_windowavg_sec,V,Markers{2},'LineWidth',2);
plot(dt_windowavg_sec,W,Markers{3},'LineWidth',2);
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('T_{avg} (s)');
ylabel('\int_{0}^{T} U_{i}(\tau) d\tau (s)');
legend('u','v','w','Location','NorthWest');
title([Site,', ',datestr(Date)]);
set(gca,'FontSize',16);
print([folder_Plots,'WindVariance_Cumulative_',Site,'_',datestr(Date),'.png'],'-dpng');

%% plot CDF of excursion durations
figure(4); clf; hold on;
Markers = {'b-','r-','g-','m-','c-'};
legend_values = cell(N_thr,1);
for i = 1:N_thr
    plot(u_thr_excursiondurations{i},...
        (1:u_thr_N_excursions(i))/u_thr_N_excursions(i),...
        Markers{i},'LineWidth',2);
    legend_values{i} = ['u_{thr} = ',int2str(u_thr(i)), ' m/s'];
end
legend(legend_values,'Location','SouthEast')
xlabel('excursion time (s)');
ylabel('CDF');
title([Site,', ',datestr(Date)]);
set(gca,'FontSize',16);
set(gca,'xscale','log');
print([folder_Plots,'WindThreshold_Excursions_',Site,'_',datestr(Date),'.png'],'-dpng');

%% plot autocorrelation of wind timeseries
figure(5); clf; hold on;
Markers = {'b-','r-'};
plot(lags_u,R_u,Markers{1},'LineWidth',2);
plot(lags_v,R_v,Markers{2},'LineWidth',2);
legend('u','v','Location','NorthEast');
xlabel('T (s)');
ylabel('R_{T}');
title([Site,', ',datestr(Date)]);
set(gca,'FontSize',16);
set(gca,'yscale','log');
print([folder_Plots,'WindAutocorrelation_',Site,'_',datestr(Date),'.png'],'-dpng');