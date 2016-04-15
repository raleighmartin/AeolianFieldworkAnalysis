%% SCRIPT TO GENERATE TIME WINDOWS FOR SALTATION FLUX / STRESS ANALYSIS
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, CreateTimeBlocks.m, 

%% initialize
clearvars;

% %% set time interval for computing wind/flux windows - full analysis
% WindowTimeInterval = duration(0,30,0); %duration of window for computations
% RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
% use_time_outside_BSNE = 0; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Analysis';
% 
% %% information about sites for analysis - full analysis
% Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
% SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};
% AnemometerType = {'Ultrasonic';'Ultrasonic';'Sonic'};
% BaseAnemometer = {'U1';'U1';'S1'};
% dt_u_s = [0.04; 0.04; 0.02]; %time interval of sonic (s)

% %% set time interval for computing wind/flux windows - Kenyon
% WindowTimeInterval = duration(0,5,0); %duration of window for computations
% RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
% use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
% SaveName = 'TimeWindows_Oceano_Kenyon';
% 
%% set time interval for computing wind/flux windows - Yue
WindowTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,1,0); %offset to use for running averages, to enrich dataset
use_time_outside_BSNE = 1; %use times outside BSNE intervals? 0 for no, 1 for yes
SaveName = 'TimeWindows_Oceano_Yue';

%% information about sites for analysis - Kenyon / Yue
Sites = {'Oceano'};
SiteNames = {'Oceano'};
AnemometerType = {'Sonic'};
BaseAnemometer = {'S1'};
dt_u_s = [0.02]; %time interval of sonic (s)

N_RunningPerProfile = floor(WindowTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times
N_Sites = length(Sites);

%% information about where to load/save data, plots, and functions
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for storing outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions
SaveData_Path = strcat(folder_AnalysisData,SaveName); %path for saving output data

%% load processed data for each site, get start times and end times for flux and wind observations
WindStartTimes = cell(N_Sites,1);
WindEndTimes = cell(N_Sites,1);
FluxStartTimes = cell(N_Sites,1);
FluxEndTimes = cell(N_Sites,1);
BSNEStartTimes = cell(N_Sites,1);
BSNEEndTimes = cell(N_Sites,1);

for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path); %load processed data
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
Date_all = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_all = cell(N_Sites,1); %lists of block start times
EndTime_all = cell(N_Sites,1); %lists of block end times


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
    [BlockStartTimes, BlockEndTimes] = ...
            CreateTimeBlocks(StartTimesIntersecting, EndTimesIntersecting, WindowTimeInterval);
    
    %Add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [BlockStartTimes_j, BlockEndTimes_j] = ...
            CreateTimeBlocks(StartTimesIntersecting+(j*RunningTimeInterval), EndTimesIntersecting, WindowTimeInterval);
        BlockStartTimes = [BlockStartTimes; BlockStartTimes_j];
        BlockEndTimes = [BlockEndTimes; BlockEndTimes_j];
    end
    BlockStartTimes = sort(BlockStartTimes);
    BlockEndTimes = sort(BlockEndTimes);
    BlockDates = datetime(BlockStartTimes.Year, BlockStartTimes.Month, BlockStartTimes.Day); %lists of dates corresponding to calculated values
    N_Blocks = length(BlockStartTimes);
  
    %% populate lists of values
    Date_all{i} = BlockDates; %lists of dates corresponding to calculated values
    StartTime_all{i} = BlockStartTimes;
    EndTime_all{i} = BlockEndTimes;
end

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','AnemometerType','BaseAnemometer','dt_u_s','N_Sites','Date_all','StartTime_all','EndTime_all'); %save file