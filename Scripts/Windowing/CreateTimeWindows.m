%% SCRIPT TO GENERATE TIME WINDOWS FOR SALTATION FLUX / STRESS ANALYSIS
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, CreateTimeBlocks.m, 

% %% initialize
% clearvars;

% %% information about sites for analysis
% Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
% SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};
% AnemometerType = {'Ultrasonic';'Ultrasonic';'Sonic'};
% AnemometerName_base = {'U1';'U1';'S1'};
% dt_wind_window = [0.04; 0.04; 0.02]; %time interval of sonic (s)
% dt_flux_window = [0.04; 0.04; 0.04]; %time interval of Wenglor (s)

% %% set time interval for computing wind/flux windows - full analysis - restricted by BSNE
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
% use_time_outside_BSNE = 0; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_30min_Restricted';

% %% set time interval for computing wind/flux windows - full analysis - unrestricted by BSNE
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_30min_Unrestricted';

%% information about sites for analysis - Kenyon / Yue
Sites = {'Oceano'};
SiteNames = {'Oceano'};
AnemometerType = {'Sonic'};
AnemometerName_base = {'S1'};
dt_wind_window = [0.02]; %time interval of sonic (s)
dt_flux_window = [0.04]; %time interval of Wenglor (s)

% %% set time interval for computing wind/flux windows - Kenyon
% WindowTimeInterval = duration(0,5,0); %duration of window for computations
% RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Kenyon';
% 
% %% set time interval for computing wind/flux windows - Yue - part 1
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,1,0); %offset to use for running averages, to enrich dataset
% min_datetime = [datetime(2015,5,26)];
% max_datetime = [datetime(2015,5,28)];
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Yue_1';

% %% set time interval for computing wind/flux windows - Yue - part 2
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,1,0); %offset to use for running averages, to enrich dataset
% min_datetime = [datetime(2015,5,29)];
% max_datetime = [datetime(2015,5,31)];
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Yue_2';

% %% set time interval for computing wind/flux windows - Yue - part 3
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,1,0); %offset to use for running averages, to enrich dataset
% min_datetime = [datetime(2015,6,1)];
% max_datetime = [datetime(2015,6,2)];
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Yue_3';
% 
%% set time interval for computing wind/flux windows - Yue - part 4
WindowTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,1,0); %offset to use for running averages, to enrich dataset
min_datetime = [datetime(2015,6,3)];
max_datetime = [datetime(2015,6,4)];
use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
SaveName = 'TimeWindows_Oceano_Yue_4';

%% number of offsets for starting times
N_RunningPerProfile = floor(WindowTimeInterval/RunningTimeInterval)-1;

%% number of sites
N_Sites = length(Sites);

%% information about where to load/save data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for storing outputs of this analysis
SaveData_Path = strcat(folder_SaveData,SaveName); %path for saving output data

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions


% %% Initialize start and end times for data intervals
% 
% %flux time information
% FluxStartTimes = cell(N_Sites,1);
% FluxEndTimes = cell(N_Sites,1);
% BSNEStartTimes = cell(N_Sites,1);
% BSNEEndTimes = cell(N_Sites,1);
% 
% %wind time information
% WindStartTimes_base = cell(N_Sites,1);
% WindEndTimes_base = cell(N_Sites,1);
% WindStartTimes_profile = cell(N_Sites,1);
% WindEndTimes_profile = cell(N_Sites,1);
% 
% %anemometer profile names and heights
% AnemometerName_profile = cell(N_Sites,1); %anemometer names
% N_Anemometers = zeros(N_Sites,1); %number of anemometers
% zU_base = cell(N_Sites,1); %base anemometer height
% zU_profile = cell(N_Sites,1); %profile anemometer heights
% 
% 
% %% load processed data for each site, get start times and end times for flux and wind observations
% for i = 1:N_Sites
%     LoadData_Path = strcat(folder_LoadData,'ProcessedData_',Sites{i});
%     load(LoadData_Path); %load processed data
%     
%     %get flux data information
%     FluxData = ProcessedData.FluxWenglor; %Wenglor flux data
%     BSNEData = ProcessedData.FluxBSNE; %BSNE data
%     FluxStartTimes{i} = [FluxData.StartTime]';
%     FluxEndTimes{i} = [FluxData.EndTime]';
%     BSNEStartTimes{i} = [BSNEData.StartTime]';
%     BSNEEndTimes{i} = [BSNEData.EndTime]';
%     
%     %get wind data information
%     WindData = ProcessedData.(AnemometerType{i}); %data for all anemometers
%     AnemometerName_profile{i} = fieldnames(WindData); %anemometer names
%     N_Anemometers(i) = length(AnemometerName_profile{i}); %number of anemometers
%     WindStartTimes_profile{i} = cell(N_Anemometers(i),1); %start times for wind data in profile
%     WindEndTimes_profile{i} = cell(N_Anemometers(i),1); %end times for wind data in profile
%     zU_profile{i} = cell(N_Anemometers(i),1); %heights of anemometers in profile
%     for j = 1:N_Anemometers(i)
%         WindStartTimes_profile{i}{j} = [WindData.(AnemometerName_profile{i}{j}).StartTime]';
%         WindEndTimes_profile{i}{j} = [WindData.(AnemometerName_profile{i}{j}).EndTime]';
%         zU_profile{i}{j} = zeros(length(WindStartTimes_profile{i}{j}),1);
%         for k = 1:length(WindStartTimes_profile{i}{j})
%             zU_profile{i}{j}(k) = [WindData.(AnemometerName_profile{i}{j})(k).z.z]';
%         end
%     end
%     WindStartTimes_base{i} = [WindData.(AnemometerName_base{i}).StartTime]';
%     WindEndTimes_base{i} = [WindData.(AnemometerName_base{i}).EndTime]';
%     zU_base{i} = zeros(length(WindStartTimes_base{i}),1);
%     for k = 1:length(WindStartTimes_base{i})
%         zU_base{i}(k) = [WindData.(AnemometerName_base{i})(k).z.z]';
%     end
% 
%     clear WindData; %remove 'WindData' to clear up memory
%     clear FluxData; %remove 'FluxData' to clear up memory
%     clear BSNEData; %remove 'BSNEData' to clear up memory
%     clear ProcessedData; %remove 'ProcessedData' to clear up memory
% end


%% initialize lists of dates, start times, and end times
Date_window = cell(N_Sites,1); %lists of dates corresponding to calculated values - all windows
StartTime_window = cell(N_Sites,1); %lists of window start times - all windows
EndTime_window = cell(N_Sites,1); %lists of window end times - all windows

%% initialize information about which windows have flux data
hasfluxdata_window = cell(N_Sites,1);

%% initialize anemometer heights
zU_base_window = cell(N_Sites,1); %height of base anemometer - all windows
zU_profile_window = cell(N_Sites,1); %height of all anemometers - all windows - will set to NaN if data don't exist


%% OBTAIN TIME WINDOWS FOR EACH SITE
for i = 1:N_Sites
       
    
    %% if not using times outside BSNE, update FluxStartTimes and FluxEndTimes to account for intersecting flux and BSNE intervals
    if use_time_outside_BSNE==0
        [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(BSNEStartTimes{i},BSNEEndTimes{i},FluxStartTimes{i},FluxEndTimes{i});
        FluxStartTimes{i} = StartTimesIntersecting;
        FluxEndTimes{i} = EndTimesIntersecting;
    end
    
    
    %% Get time windows based on flux/wind intersecting time intervals
    % get flux/wind intersecting time intervals
    [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(WindStartTimes_base{i},WindEndTimes_base{i},FluxStartTimes{i},FluxEndTimes{i});

    %create time blocks based on intersecting time intervals
    [WindowStartTimes_flux, WindowEndTimes_flux] = ...
            CreateTimeBlocks(StartTimesIntersecting, EndTimesIntersecting, WindowTimeInterval);
    
    %Add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [WindowStartTimes_j, WindowEndTimes_j] = ...
            CreateTimeBlocks(StartTimesIntersecting+(j*RunningTimeInterval), EndTimesIntersecting, WindowTimeInterval);
        WindowStartTimes_flux = [WindowStartTimes_flux; WindowStartTimes_j];
        WindowEndTimes_flux = [WindowEndTimes_flux; WindowEndTimes_j];
    end
    WindowStartTimes_flux = sort(WindowStartTimes_flux);
    WindowEndTimes_flux = sort(WindowEndTimes_flux);
    
    %Get associated dates and number of windows
    WindowDates_flux = datetime(WindowStartTimes_flux.Year, WindowStartTimes_flux.Month, WindowStartTimes_flux.Day); %lists of dates corresponding to calculated values
    
    %limit to only selected times, if necessary
    if exist('min_datetime','var')
        ind_window_flux = find(WindowDates_flux>=min_datetime & WindowDates_flux<=max_datetime);
        WindowStartTimes_flux = WindowStartTimes_flux(ind_window_flux);
        WindowEndTimes_flux = WindowEndTimes_flux(ind_window_flux);
        WindowDates_flux = WindowDates_flux(ind_window_flux);
    end
    
    %Get number of windows
    N_Windows_flux = length(WindowStartTimes_flux);
    
    
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
    
    %limit to only selected times, if necessary
    if exist('min_datetime','var')
        ind_window_noflux = find(WindowDates_noflux>=min_datetime & WindowDates_noflux<=max_datetime);
        WindowStartTimes_noflux = WindowStartTimes_noflux(ind_window_noflux);
        WindowEndTimes_noflux = WindowEndTimes_noflux(ind_window_noflux);
        WindowDates_noflux = WindowDates_noflux(ind_window_noflux);
    end
    
    %Get number of windows
    N_Windows_noflux = length(WindowStartTimes_noflux);
    
    
    %% combine start and end times for windows with and without flux data
    hasfluxdata = [ones(size(WindowStartTimes_flux)); zeros(size(WindowStartTimes_noflux))]; %list of zeros for windows with flux data, list of ones for windows with no flux data
    WindowStartTimes = [WindowStartTimes_flux; WindowStartTimes_noflux]; %combined window start times
    WindowEndTimes = [WindowEndTimes_flux; WindowEndTimes_noflux]; %combined window end times
    WindowDates = [WindowDates_flux; WindowDates_noflux]; %combined window dates
    N_Windows = length(WindowStartTimes); %total number of windows
    
    % sort combined values by start time
    [~, ind_StartTime_sort] = sort(WindowStartTimes);
    WindowStartTimes = WindowStartTimes(ind_StartTime_sort);
    WindowEndTimes = WindowEndTimes(ind_StartTime_sort);
    WindowDates = WindowDates(ind_StartTime_sort);
    hasfluxdata = hasfluxdata(ind_StartTime_sort);
    
    %% populate lists of values
    Date_window{i} = WindowDates;
    StartTime_window{i} = WindowStartTimes;
    EndTime_window{i} = WindowEndTimes;  
    hasfluxdata_window{i} = hasfluxdata;
    
    
    %% Now, go through windows to get heights for base anemometer and anemometer profile
    zU_base_window{i} = zeros(N_Windows,1)*NaN; %height of base anemometer for all windows
    zU_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %height of all anemometers for all windows
    
    %anemometer profile
    for j = 1:N_Anemometers(i) %go through each anemometer
        for k = 1:length(WindStartTimes_profile{i}{j}) %go through start times for data from each anemometer
            ind_window_flux = find(isbetween(WindowStartTimes,WindStartTimes_profile{i}{j}(k),WindEndTimes_profile{i}{j}(k)).* ...
                isbetween(WindowEndTimes,WindStartTimes_profile{i}{j}(k),WindEndTimes_profile{i}{j}(k))); %figure out which 30-minute windows correspond to anemometer window
            zU_profile_window{i}(ind_window_flux,j) = zU_profile{i}{j}(k); %add these heights to list
        end
    end
    
    %base anemometer
    for k = 1:length(WindStartTimes_base{i})
        ind_window_flux = find(isbetween(WindowStartTimes,WindStartTimes_base{i}(k),WindEndTimes_base{i}(k)).* ...
            isbetween(WindowEndTimes,WindStartTimes_base{i}(k),WindEndTimes_base{i}(k)));
        zU_base_window{i}(ind_window_flux) = zU_base{i}(k);
    end
end

% SAVE DATA
save(SaveData_Path,'Site*','N_Sites','AnemometerType','AnemometerName_base','AnemometerName_profile','*window'); %save file