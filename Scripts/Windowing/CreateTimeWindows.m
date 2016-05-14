%% SCRIPT TO GENERATE TIME WINDOWS FOR SALTATION FLUX / STRESS ANALYSIS
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, CreateTimeBlocks.m, 

%% initialize
clearvars;

%% set time interval for computing wind/flux windows - full analysis
WindowTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
use_time_outside_BSNE = 0; %use times outside BSNE intervals? 0 for no, 1 for yes
SaveName = 'TimeWindows_30min';

%% information about sites for analysis - full analysis
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};
AnemometerType = {'Ultrasonic';'Ultrasonic';'Sonic'};
BaseAnemometer = {'U1';'U1';'S1'};
dt_wind_window = [0.04; 0.04; 0.02]; %time interval of sonic (s)
dt_flux_window = [0.04; 0.04; 0.04]; %time interval of Wenglor (s)

% %% set time interval for computing wind/flux windows - Kenyon
% WindowTimeInterval = duration(0,5,0); %duration of window for computations
% RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Kenyon';
% 
% %% set time interval for computing wind/flux windows - Yue
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,1,0); %offset to use for running averages, to enrich dataset
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Yue';
% 
% %% information about sites for analysis - Kenyon / Yue
% Sites = {'Oceano'};
% SiteNames = {'Oceano'};
% AnemometerType = {'Sonic'};
% BaseAnemometer = {'S1'};
% dt_u_s = [0.02]; %time interval of sonic (s)

N_RunningPerProfile = floor(WindowTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times
N_Sites = length(Sites);

%% information about where to load/save data, plots, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for storing outputs of this analysis
SaveData_Path = strcat(folder_SaveData,SaveName); %path for saving output data

%% load processed data for each site, get start times and end times for flux and wind observations
WindStartTimes = cell(N_Sites,1);
WindEndTimes = cell(N_Sites,1);
FluxStartTimes = cell(N_Sites,1);
FluxEndTimes = cell(N_Sites,1);
BSNEStartTimes = cell(N_Sites,1);
BSNEEndTimes = cell(N_Sites,1);

for i = 1:N_Sites
    LoadData_Path = strcat(folder_LoadData,'ProcessedData_',Sites{i});
    load(LoadData_Path); %load processed data
    WindData = ProcessedData.(AnemometerType{i}).(BaseAnemometer{i}); %data only for base anemometer
    FluxData = ProcessedData.FluxWenglor; %Wenglor flux data
    BSNEData = ProcessedData.FluxBSNE; %BSNE data
    WindStartTimes{i} = [WindData.StartTime]';
    WindEndTimes{i} = [WindData.EndTime]';
    FluxStartTimes{i} = [FluxData.StartTime]';
    FluxEndTimes{i} = [FluxData.EndTime]';
    BSNEStartTimes{i} = [BSNEData.StartTime]';
    BSNEEndTimes{i} = [BSNEData.EndTime]';
    
    clear WindData; %remove 'WindData' to clear up memory
    clear FluxData; %remove 'FluxData' to clear up memory
    clear BSNEData; %remove 'BSNEData' to clear up memory
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
end


%% initialize variable lists
Date_window = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_window = cell(N_Sites,1); %lists of window start times
EndTime_window = cell(N_Sites,1); %lists of window end times


%% GET TIME WINDOWS FOR EACH SITE
for i = 1:N_Sites
    
    if use_time_outside_BSNE==0
        %update FluxStartTimes and Flux EndTimes to account for intersecting flux and BSNE intervals
        [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(BSNEStartTimes{i},BSNEEndTimes{i},FluxStartTimes{i},FluxEndTimes{i});
        FluxStartTimes{i} = StartTimesIntersecting;
        FluxEndTimes{i} = EndTimesIntersecting;
    end
        
    %get start and end times for intersecting flux and wind intervals
    [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(WindStartTimes{i},WindEndTimes{i},FluxStartTimes{i},FluxEndTimes{i});
    
    %create time blocks based on intersecting time intervals
    [WindowStartTimes, WindowEndTimes] = ...
            CreateTimeBlocks(StartTimesIntersecting, EndTimesIntersecting, WindowTimeInterval);
    
    %Add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [WindowStartTimes_j, WindowEndTimes_j] = ...
            CreateTimeBlocks(StartTimesIntersecting+(j*RunningTimeInterval), EndTimesIntersecting, WindowTimeInterval);
        WindowStartTimes = [WindowStartTimes; WindowStartTimes_j];
        WindowEndTimes = [WindowEndTimes; WindowEndTimes_j];
    end
    WindowStartTimes = sort(WindowStartTimes);
    WindowEndTimes = sort(WindowEndTimes);
    WindowDates = datetime(WindowStartTimes.Year, WindowStartTimes.Month, WindowStartTimes.Day); %lists of dates corresponding to calculated values
    N_Windows = length(WindowStartTimes);
  
    %% populate lists of values
    Date_window{i} = WindowDates; %lists of dates corresponding to calculated values
    StartTime_window{i} = WindowStartTimes;
    EndTime_window{i} = WindowEndTimes;
end

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','N_Sites','AnemometerType','BaseAnemometer','*window'); %save file