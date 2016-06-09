%% initialize
clearvars;

%% parameter values
delta_t_lowpass_window = ... %time intervals for low-pass filtering
    [duration(0,0,5),...
    duration(0,0,25)...
    duration(0,0,100)];
N_lowpass_window = length(delta_t_lowpass_window); %number of time intervals for low-pass filtering

%% information about where to load/save data, plots, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for windowing data
DataWindow_Path = strcat(folder_LoadData,'DataWindows_30min'); %path for loading data windows
SaveData_Path = strcat(folder_SaveData,'LowPassWindows_30min'); %path for saving output data

%% load time windows
load(DataWindow_Path);

%% initialize variable lists
n_lowpass_window = cell(N_Sites,1); %low-pass count rates values
t_flux_lowpass_window = cell(N_Sites,1); %low-pass counts times
u_lowpass_window = cell(N_Sites,1); %low-pass u values
t_wind_lowpass_window = cell(N_Sites,1); %low-pass wind times

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    N_Windows = length(StartTime_window{i}); %get number of windows
        
    %% initialize lists of values
    n_lowpass_window{i} = cell(N_Windows,N_lowpass_window); %low-pass count rates values
    t_flux_lowpass_window{i} = cell(N_Windows,1); %low-pass counts times
    u_lowpass_window{i} = cell(N_Windows,N_lowpass_window); %low-pass u values
    t_wind_lowpass_window{i} = cell(N_Windows,1); %low-pass wind times
    
    %% go through time windows
    for j = 1:N_Windows

        %% display processing status
        processing_status = [SiteNames{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]
        
        %% get needed values for analysis
        n = n_int_window{i}{j}; %interpolated counts rate
        t_flux = t_flux_int_window{i}{j}; %times for interpolated flux
        N_zW = N_zW_window{i}(j); %number of Wenglors
        u = u_int_window{i}{j}; %rotated interpolated streamwise velocities
        t_wind = t_wind_int_window{i}{j}; %times for interpolated wind
        
        %% go through low-pass durations
        for k = 1:N_lowpass_window
        
            %% compute low-pass flux
            n_lowpass = zeros(length(n),N_zW); %initialize low-pass counts timeseries
            for l=1:N_zW %go through all Wenglors
                n_lowpass(:,l) = LowPassFilter(n(:,l), dt_flux_window(i), 1/seconds(delta_t_lowpass_window(k))); %compute low-pass counts timeseries
            end

            %% compute low-pass wind
            u_lowpass = LowPassFilter(u, dt_wind_window(i), 1/seconds(delta_t_lowpass_window(k))); %compute low-pass wind timeseries - using interpolated values

            %% add to lists
            n_lowpass_window{i}{j,k} = n_lowpass; %save low-pass counts rates
            u_lowpass_window{i}{j,k} = u_lowpass; %save low-pass wind speeds
        end
        t_flux_lowpass_window{i}{j} = t_flux; %save time for counts
        t_wind_lowpass_window{i}{j} = t_wind; %save time for wind
    end
end

% SAVE DATA
save(SaveData_Path,'SiteNames','N_Sites','zU_window','*lowpass_window'); %save reduced file in ProcessedData folder