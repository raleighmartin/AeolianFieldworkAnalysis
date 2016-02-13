%% SCRIPT TO GENERATE TIME WINDOWS FOR SALTATION FLUX / STRESS ANALYSIS
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, CreateTimeBlocks.m, 

%% initialize
clearvars;

%% set time interval for computing wind/flux windows
WindowTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,30,0); %do not use overlapping windows
%RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
N_RunningPerProfile = floor(WindowTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times

%% information about sites for analysis
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
AnemometerType = {'Ultrasonic';'Ultrasonic';'Sonic'};
BaseAnemometer = {'U1';'U1';'S1'};
dt_u_s = [0.04; 0.04; 0.02]; %time interval of sonic (s)
N_Sites = length(Sites);

%% information about where to load/save data, plots, and functions
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for storing outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
SaveData_Path = strcat(folder_AnalysisData,'TimeWindows'); %path for saving output data
addpath(folder_Functions); %point MATLAB to location of functions


% %% load processed and metadata for each site, add to structured arrays of all data and metadata
% Data = cell(N_Sites,1);
% for i = 1:N_Sites
%     ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
%     Data{i} = load(ProcessedData_Path); %load processed data
% end

%% load processed data for each site, get start times and end times for flux and wind observations
WindStartTimes = cell(N_Sites,1);
WindEndTimes = cell(N_Sites,1);
FluxStartTimes = cell(N_Sites,1);
FluxEndTimes = cell(N_Sites,1);

for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path); %load processed data
    WindData = ProcessedData.(AnemometerType{i}).(BaseAnemometer{i}); %data only for base anemometer
    FluxData = ProcessedData.FluxWenglor; %Wenglor flux data
    WindStartTimes{i} = [WindData.StartTime]';
    WindEndTimes{i} = [WindData.EndTime]';
    FluxStartTimes{i} = [FluxData.StartTime]';
    FluxEndTimes{i} = [FluxData.EndTime]';
    
    clear WindData; %remove 'WindData' to clear up memory
    clear FluxData; %remove 'FluxData' to clear up memory
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
end


%% initialize variable lists
Date_all = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_all = cell(N_Sites,1); %lists of block start times
EndTime_all = cell(N_Sites,1); %lists of block end times


%% GET TIME WINDOWS FOR EACH SITE
for i = 1:N_Sites
        
%     %% choose anemometer type based on site of interest
%     if strcmp(Sites{i},'Oceano')
%         AnemometerType = 'Sonic';
%         BaseAnemometer = 'S1';
%         dt_u_s = 0.02; %time interval of sonic (s)
%     elseif strcmp(Sites{i},'RanchoGuadalupe')
%         AnemometerType = 'Ultrasonic';
%         BaseAnemometer = 'U1';
%         dt_u_s = 0.04; %time interval of sonic (s)
%     elseif strcmp(Sites{i},'Jericoacoara')
%         AnemometerType = 'Ultrasonic';
%         BaseAnemometer = 'U1';
%         dt_u_s = 0.04; %time interval of sonic (s)
%     end
%     
%     
%     %% extract wind and flux data from overall processed data file
%     WindDataAll = Data{i}.ProcessedData.(AnemometerType); %all wind data
%     Anemometer_profile = fieldnames(WindDataAll); %get names of anemometers in profile
%     N_Anemometers = length(Anemometer_profile); %get number of anemometers in profile
%     WindDataBase = WindDataAll.(BaseAnemometer); %data only for base anemometer
%     FluxData = Data{i}.ProcessedData.FluxWenglor; %Wenglor flux data
%     
%     %% get start times and end times for flux and wind observations, using times for base anemometer
%     WindStartTimes = [WindData{i}.StartTime]';
%     WindEndTimes = [WindData{i}.EndTime]';
%     FluxStartTimes = [FluxData{i}.StartTime]';
%     FluxEndTimes = [FluxData{i}.EndTime]';
    
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
    N_Blocks = length(BlockStartTimes);
    
    
    %% populate lists of values
    Date_all{i} = datetime(BlockStartTimes.Year, BlockStartTimes.Month, BlockStartTimes.Day); %lists of dates corresponding to calculated values
    StartTime_all{i} = BlockStartTimes;
    EndTime_all{i} = BlockEndTimes;
end

% SAVE DATA
save(SaveData_Path,'Sites','Date_all','StartTime_all','EndTime_all'); %save file