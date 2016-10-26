%% SCRIPT TO GENERATE WINDOWS OF DISTANCE SENSOR DATA FOR RIPPLE
 
%% initialize
clearvars;
close all;

%% information about sensors
UpwindSensor = 'L2';
DownwindSensor = 'L1';

%% set times for analysis
StartTimes_Analysis = [... 
    datetime(2015,5,15,12,30,0);...
    datetime(2015,5,18,11,45,0);...
    datetime(2015,5,19,13,00,0);...
    datetime(2015,5,23,12,30,0);...
    datetime(2015,6,2,14,28,0);...
    datetime(2015,6,3,13,00,0)];

EndTimes_Analysis = [... 
    datetime(2015,5,15,18,03,0);...
    datetime(2015,5,18,17,00,0);...
    datetime(2015,5,19,17,45,0);...
    datetime(2015,5,23,18,10,0);...
    datetime(2015,6,2,18,16,0);...
    datetime(2015,6,3,18,34,0)];

N_Intervals = length(StartTimes_Analysis); %get number of intervals

%% load data - may need to change path
folder_Data = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
filename_Data = 'DistanceSensor_Oceano.mat'; %file name for retrieving data
load([folder_Data,filename_Data]); %load data

%% initialize cell arrays of distance sensor data for analysis
zU = cell(N_Intervals,1); %upwind sensor - elevations
tU = cell(N_Intervals,1); %upwind sensor - times
zD = cell(N_Intervals,1); %downwind sensor - elevations
tD = cell(N_Intervals,1); %downwind sensor - times

%% go through intervals
for i = 1:N_Intervals
    
    %% extract data for specific times for intervals - using function
    [zU_Interval, tU_Interval, Interval_N, ~] = ...
        ExtractVariableTimeInterval(Distance.(UpwindSensor),StartTimes_Analysis(i),EndTimes_Analysis(i),'z','int','int');
    [zD_Interval, tD_Interval, ~, ~] = ...
        ExtractVariableTimeInterval(Distance.(DownwindSensor),StartTimes_Analysis(i),EndTimes_Analysis(i),'z','int','int');

    %% add values to cell arrays
    zU{i} = zU_Interval;
    tU{i} = tU_Interval;
    zD{i} = zD_Interval;
    tD{i} = tD_Interval;
    
    %% get corresponding data for full intervals to plot
    date_Interval = Distance.(UpwindSensor)(Interval_N).Date;
    zU_full = Distance.(UpwindSensor)(Interval_N).z.int;
    tU_full = Distance.(UpwindSensor)(Interval_N).t.int;
    zD_full = Distance.(DownwindSensor)(Interval_N).z.int;
    tD_full = Distance.(DownwindSensor)(Interval_N).t.int;
    
    %% plot timeseries, with dashed lines showing data included in analysis
    figure(i); clf; hold on; %initialize figure
    plot(tU_full,zU_full,tD_full,zD_full); %plot full timeseries
    ylim_plot = ylim; %get vertical extent of plot
    
    %plot dashed lines for beginning and end of analysis periods
    plot([StartTimes_Analysis(i), StartTimes_Analysis(i)],ylim_plot,'k--',...
        [EndTimes_Analysis(i), EndTimes_Analysis(i)],ylim_plot,'k--');
     
    %annotate plot
    legend('Upwind','Downwind','Location','NorthOutside');
    ylabel('z (mm)');
    set(gca,'FontSize',16);
    
    %save plot
    print(['DistanceSensorData_',datestr(date_Interval),'.png'],'-dpng');
end

%% save data
save('DistanceSensor_Oceano_Analysis','StartTimes_Analysis','EndTimes_Analysis','zU','tU','zD','tD');