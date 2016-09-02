%% SCRIPT TO GENERATE WINDOWS OF FLUX ACTIVITY, SALTATION FLUX, AND STRESS VALUES FOR ANALYSIS
% Delta_t = measurement interval duration (longer)
% delta_t = sampling interval duration (shorter)

%% initialize
clearvars;

%% parameter values
delta_t_avg_window = ... %time intervals for window averaging
    [duration(0,0,0.2),...
    duration(0,0,1)...
    duration(0,0,5)];
N_delta_t_avg = length(delta_t_avg_window); %number of time intervals for window averaging

%% information about where to load/save data, plots, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for windowing data
DataWindow_Path = strcat(folder_LoadData,'DataWindows_30min'); %path for loading data windows
SaveData_Path = strcat(folder_SaveData,'WindowAverageWindows_30min'); %path for saving output data

%% load time windows
load(DataWindow_Path);

%% initialize variable lists
t_avg_window = cell(N_Sites,1); %times for window average
n_avg_window = cell(N_Sites,1); %window-averaged count rates values
u_avg_window = cell(N_Sites,1); %window-averaged u values
ind_flux_err_avg_window = cell(N_Sites,1); %indices of error points in flux window-average
ind_wind_err_avg_window = cell(N_Sites,1); %indices of error points in wind window-average

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    N_Windows = length(StartTime_window{i}); %get number of windows
        
    %% initialize lists of values
    t_avg_window{i} = cell(N_Windows,N_delta_t_avg); %times for window average
    n_avg_window{i} = cell(N_Windows,N_delta_t_avg); %window-averaged count rates values
    u_avg_window{i} = cell(N_Windows,N_delta_t_avg); %window-averaged u values
    ind_flux_err_avg_window{i} = cell(N_Windows,N_delta_t_avg); %indices of error points in flux window-average
    ind_wind_err_avg_window{i} = cell(N_Windows,N_delta_t_avg); %indices of error points in window-average
    
    %% go through time windows
    for j = 1:N_Windows

        %% display processing status
        processing_status = [SiteNames{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]
        
        %% get needed values for analysis
        n = n_int_window{i}{j}; %interpolated counts rate
        t_flux = t_flux_int_window{i}{j}; %times for interpolated flux
        N_zW = N_zW_window{i}(j); %number of Wenglors
        ind_flux_err = ind_flux_err_window{i}{j}; %list of error time indices for flux
        u = u_int_window{i}{j}; %rotated interpolated streamwise velocities
        t_wind = t_wind_int_window{i}{j}; %times for interpolated wind
        ind_wind_err = ind_wind_err_window{i}{j}; %list of error time indices for wind
        
        %% go through window-average durations
        for k = 1:N_delta_t_avg
        
            %% compute window-averaged flux
            n_avg_1 = window_average(n(:,1), t_flux, delta_t_avg_window(k)); %compute sample window average counts timeseries to get number of values
            n_avg = zeros(length(n_avg_1),N_zW); %initialize full window averaged counts timeseries
            for l=1:N_zW %go through all Wenglors
                n_avg(:,l) = window_average(n(:,l), t_flux, delta_t_avg_window(k)); %compute window average counts timeseries
            end

            %% compute window-averaged wind
            [u_avg, t_avg] = window_average(u, t_wind, delta_t_avg_window(k)); %compute window average wind timeseries - using interpolated values

            %% indices of error points in window-averaged flux
            ind_err = zeros(length(n_avg),1); %initialize indices of error points in raw timeseries
            ind_err(ind_flux_err) = 1; %set error indices to 1
            ind_err_avg = window_average(ind_err, t_flux, delta_t_avg_window(k)); %compute window average 
            ind_flux_avg_err = find(ind_err_avg>0); %find window-averaged points > 0

            %% indices of error points in window-averaged wind
            ind_err = zeros(size(t_avg)); %initialize indices of error points in raw timeseries
            ind_err(ind_wind_err) = 1; %set error indices to 1
            ind_err_avg = window_average(ind_err, t_wind, delta_t_avg_window(k)); %compute window average 
            ind_wind_avg_err = find(ind_err_avg>0); %find window-averaged points > 0

            %% add to lists
            t_avg_window{i}{j,k} = t_avg; %save window-average times
            n_avg_window{i}{j,k} = n_avg; %save window-averaged counts rates
            u_avg_window{i}{j,k} = u_avg; %save window-average wind speeds
            ind_flux_err_avg_window{i}{j,k} = ind_flux_avg_err; %save list of error indices for window-averaged flux timeseries
            ind_wind_err_avg_window{i}{j,k} = ind_wind_avg_err; %save list of error indices for window-averaged wind timeseries
        end
    end
end

% SAVE DATA
save(SaveData_Path,'SiteNames','N_Sites','zU_window','*avg_window'); %save reduced file in ProcessedData folder