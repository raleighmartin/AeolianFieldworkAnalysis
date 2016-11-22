%% SCRIPT TO GENERATE TIME WINDOWS FOR SALTATION FLUX / STRESS ANALYSIS
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, CreateTimeBlocks.m, 

%% initialize
clearvars;

% %% set time interval for computing wind/flux windows - full analysis - restricted by BSNE
WindowTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
use_time_outside_BSNE = 0; %use times outside BSNE intervals? 0 for no, 1 for yes
SaveName = 'TimeWindows_30min_Restricted';

% % set time interval for computing wind/flux windows - full analysis - unrestricted by BSNE
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_30min_Unrestricted';

%% information about sites for analysis
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};
AnemometerType = {'Ultrasonic';'Ultrasonic';'Sonic'};
AnemometerName_base = {'U1';'U1';'S1'};
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
% AnemometerName_base = {'S1'};
% dt_u_s = [0.02]; %time interval of sonic (s)

N_RunningPerProfile = floor(WindowTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times
N_Sites = length(Sites);

%% information about where to load/save data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for storing outputs of this analysis
SaveData_Path = strcat(folder_SaveData,SaveName); %path for saving output data

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Initialize start and end times for data intervals

%flux time information
FluxStartTimes = cell(N_Sites,1);
FluxEndTimes = cell(N_Sites,1);
BSNEStartTimes = cell(N_Sites,1);
BSNEEndTimes = cell(N_Sites,1);

%wind time information
WindStartTimes_base = cell(N_Sites,1);
WindEndTimes_base = cell(N_Sites,1);
WindStartTimes_profile = cell(N_Sites,1);
WindEndTimes_profile = cell(N_Sites,1);

%anemometer profile names and heights
AnemometerName_profile = cell(N_Sites,1); %anemometer names
N_Anemometers = zeros(N_Sites,1); %number of anemometers
zU_base = cell(N_Sites,1); %base anemometer height
zU_profile = cell(N_Sites,1); %profile anemometer heights

%% load processed data for each site, get start times and end times for flux and wind observations
for i = 1:N_Sites
    LoadData_Path = strcat(folder_LoadData,'ProcessedData_',Sites{i});
    load(LoadData_Path); %load processed data
    
    %get flux data information
    FluxData = ProcessedData.FluxWenglor; %Wenglor flux data
    BSNEData = ProcessedData.FluxBSNE; %BSNE data
    FluxStartTimes{i} = [FluxData.StartTime]';
    FluxEndTimes{i} = [FluxData.EndTime]';
    BSNEStartTimes{i} = [BSNEData.StartTime]';
    BSNEEndTimes{i} = [BSNEData.EndTime]';
    
    %get wind data information
    WindData = ProcessedData.(AnemometerType{i}); %data for all anemometers
    AnemometerName_profile{i} = fieldnames(WindData); %anemometer names
    N_Anemometers(i) = length(AnemometerName_profile{i}); %number of anemometers
    WindStartTimes_profile{i} = cell(N_Anemometers(i),1); %start times for wind data in profile
    WindEndTimes_profile{i} = cell(N_Anemometers(i),1); %end times for wind data in profile
    zU_profile{i} = cell(N_Anemometers(i),1); %heights of anemometers in profile
    for j = 1:N_Anemometers(i)
        WindStartTimes_profile{i}{j} = [WindData.(AnemometerName_profile{i}{j}).StartTime]';
        WindEndTimes_profile{i}{j} = [WindData.(AnemometerName_profile{i}{j}).EndTime]';
        zU_profile{i}{j} = zeros(length(WindStartTimes_profile{i}{j}),1);
        for k = 1:length(WindStartTimes_profile{i}{j})
            zU_profile{i}{j}(k) = [WindData.(AnemometerName_profile{i}{j})(k).z.z]';
        end
    end
    WindStartTimes_base{i} = [WindData.(AnemometerName_base{i}).StartTime]';
    WindEndTimes_base{i} = [WindData.(AnemometerName_base{i}).EndTime]';
    zU_base{i} = zeros(length(WindStartTimes_base{i}),1);
    for k = 1:length(WindStartTimes_base{i})
        zU_base{i}(k) = [WindData.(AnemometerName_base{i})(k).z.z]';
    end

    clear WindData; %remove 'WindData' to clear up memory
    clear FluxData; %remove 'FluxData' to clear up memory
    clear BSNEData; %remove 'BSNEData' to clear up memory
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
end


%% initialize window date time lists - windows with flux data
Date_window = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_window = cell(N_Sites,1); %lists of window start times
EndTime_window = cell(N_Sites,1); %lists of window end times

%% initialize window date time lists - windows without flux data
Date_window_noflux = cell(N_Sites,1); %lists of dates corresponding to calculated values - no flux windows
StartTime_window_noflux = cell(N_Sites,1); %lists of window start times - no flux windows
EndTime_window_noflux = cell(N_Sites,1); %lists of window end times - no flux windows

%% get heights of anemometers for windows with flux data
zU_base_window = cell(N_Sites,1); %height of base anemometer
zU_profile_window = cell(N_Sites,1); %height of all anemometers - set to NaN if data don't exist

%% get heights of anemometers for windows without flux data
zU_base_window_noflux = cell(N_Sites,1); %height of base anemometer
zU_profile_window_noflux = cell(N_Sites,1); %height of all anemometers - set to NaN if data don't exist

%% GET TIME WINDOWS FOR EACH SITE
for i = 1:N_Sites
    
    %% get time intervals for flux
    if use_time_outside_BSNE==0
        %update FluxStartTimes and Flux EndTimes to account for intersecting flux and BSNE intervals
        [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(BSNEStartTimes{i},BSNEEndTimes{i},FluxStartTimes{i},FluxEndTimes{i});
        FluxStartTimes{i} = StartTimesIntersecting;
        FluxEndTimes{i} = EndTimesIntersecting;
    end
    
    %% get flux/wind intersecting time intervals
    [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(WindStartTimes_base{i},WindEndTimes_base{i},FluxStartTimes{i},FluxEndTimes{i});
    
    
    %% Get time windows based on flux/wind intersecting time intervals
    
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
    
    %Get associated dates and number of windows
    WindowDates = datetime(WindowStartTimes.Year, WindowStartTimes.Month, WindowStartTimes.Day); %lists of dates corresponding to calculated values
    N_Windows = length(WindowStartTimes);

    
    %% Get additional time windows without flux data
    
    % first, specify wind time windows without flux
    WindStartTimes_noflux = [];
    WindEndTimes_noflux = [];
    for j = 1:length(WindStartTimes_base{i})
        %get indices of intersecting intervals within wind interval start/end time range
        IntervalInd = find(isbetween(StartTimesIntersecting,WindStartTimes_base{i}(j),WindEndTimes_base{i}(j)).* ...
            isbetween(EndTimesIntersecting,WindStartTimes_base{i}(j),WindEndTimes_base{i}(j)));
              
        %add start times to list (start time of wind interval, then intersecting interval end times within wind interval
        WindStartTimes_noflux = [WindStartTimes_noflux; WindStartTimes_base{i}(j); EndTimesIntersecting(IntervalInd)];

        %add end times to list (intersecting interval start times within wind interval, then end time of wind interval 
        WindEndTimes_noflux = [WindEndTimes_noflux; StartTimesIntersecting(IntervalInd); WindEndTimes_base{i}(j)];
    end
    
    %then, create time blocks based no flux data time intervals
    [WindowStartTimes_noflux, WindowEndTimes_noflux] = ...
            CreateTimeBlocks(WindStartTimes_noflux, WindEndTimes_noflux, WindowTimeInterval);
    
    %add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [WindowStartTimes_j, WindowEndTimes_j] = ...
            CreateTimeBlocks(WindStartTimes_noflux+(j*RunningTimeInterval), WindEndTimes_noflux, WindowTimeInterval);
        WindowStartTimes_noflux = [WindowStartTimes_noflux; WindowStartTimes_j];
        WindowEndTimes_noflux = [WindowEndTimes_noflux; WindowEndTimes_j];
    end
    WindowStartTimes_noflux = sort(WindowStartTimes_noflux);
    WindowEndTimes_noflux = sort(WindowEndTimes_noflux);
    
    %Get associated dates and number of windows
    WindowDates_noflux = datetime(WindowStartTimes_noflux.Year, WindowStartTimes_noflux.Month, WindowStartTimes_noflux.Day); %lists of dates corresponding to calculated values
    N_Windows_noflux = length(WindowStartTimes_noflux);
    
    %% populate lists of values
    Date_window{i} = WindowDates;
    StartTime_window{i} = WindowStartTimes;
    EndTime_window{i} = WindowEndTimes;
    
    Date_window_noflux{i} = WindowDates_noflux;
    StartTime_window_noflux{i} = WindowStartTimes_noflux;
    EndTime_window_noflux{i} = WindowEndTimes_noflux;  
    
    %% Now, go through windows to get data for anemometer profile
    zU_base_window{i} = zeros(N_Windows,1)*NaN; %height of base anemometer for windows with flux data
    zU_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %height of all anemometers for windows with flux data - set to NaN if data don't exist
    zU_base_window_noflux{i} = zeros(N_Windows_noflux,1)*NaN; %height of base anemometer for windows without flux data
    zU_profile_window_noflux{i} = zeros(N_Windows_noflux,N_Anemometers(i))*NaN; %height of all anemometers for windows without flux data - set to NaN if data don't exist
    
    %anemometer profile
    for j = 1:N_Anemometers(i)
        for k = 1:length(WindStartTimes_profile{i}{j})
            ind_window = find(isbetween(StartTime_window{i},WindStartTimes_profile{i}{j}(k),WindEndTimes_profile{i}{j}(k)).* ...
                isbetween(EndTime_window{i},WindStartTimes_profile{i}{j}(k),WindEndTimes_profile{i}{j}(k)));
            zU_profile_window{i}(ind_window,j) = zU_profile{i}{j}(k);
            ind_window_noflux = find(isbetween(StartTime_window_noflux{i},WindStartTimes_profile{i}{j}(k),WindEndTimes_profile{i}{j}(k)).* ...
                isbetween(EndTime_window_noflux{i},WindStartTimes_profile{i}{j}(k),WindEndTimes_profile{i}{j}(k)));
            zU_profile_window_noflux{i}(ind_window_noflux,j) = zU_profile{i}{j}(k);
        end
    end
    
    %base anemometer
    for k = 1:length(WindStartTimes_base{i})
        ind_window = find(isbetween(StartTime_window{i},WindStartTimes_base{i}(k),WindEndTimes_base{i}(k)).* ...
            isbetween(EndTime_window{i},WindStartTimes_base{i}(k),WindEndTimes_base{i}(k)));
        zU_base_window{i}(ind_window) = zU_base{i}(k);
        ind_window_noflux = find(isbetween(StartTime_window_noflux{i},WindStartTimes_base{i}(k),WindEndTimes_base{i}(k)).* ...
            isbetween(EndTime_window_noflux{i},WindStartTimes_base{i}(k),WindEndTimes_base{i}(k)));
        zU_base_window_noflux{i}(ind_window_noflux) = zU_base{i}(k);
    end
end

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','N_Sites','WindowTimeInterval','RunningTimeInterval','AnemometerType','AnemometerName_base','AnemometerName_profile','use_time_outside_BSNE','*window','*window_noflux'); %save file