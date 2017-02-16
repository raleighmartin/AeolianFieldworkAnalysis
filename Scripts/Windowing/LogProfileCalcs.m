%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
kappa = 0.4; %von Karman parameter
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
g = 9.8; %gravity m/s^2
zU_max = 2.5; %maximum anemometer height for profile fit (m)

%% folders for loading and saving data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

%% paths for loading and saving data - restricted
LoadData_Path = strcat(folder_LoadData,'WindProfiles_30min_Restricted'); %path for 30 minute data - for flux law analysis
SaveData_Path = strcat(folder_SaveData,'LogProfiles_30min_Restricted'); %path for 30 minute data - for flux law analysis

% %% paths for loading and saving data - unrestricted
% create_profile = 1; %toggle whether to do calcs for full profile
% LoadData_Path = strcat(folder_LoadData,'WindProfiles_30min_Restricted'); %path for 30 minute data - for thresholds analysis
% SaveData_Path = strcat(folder_SaveData,'LogProfiles_30min_Restricted'); %path for 30 minute data - for thresholds analysis

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize variable lists
zsLog_all = cell(N_Sites,1); %roughness height from wind profile
sigma_zsLog_all = cell(N_Sites,1); %roughness height from wind profile - uncertainty
ustLog_all = cell(N_Sites,1); %shear velocity from wind profile
sigma_ustLog_all = cell(N_Sites,1); %shear velocity from wind profile - uncertainty
tauLog_all = cell(N_Sites,1); %shear stress from wind profile 
sigma_tauLog_all = cell(N_Sites,1); %shear stress from wind profile - uncertainty
   
%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% get number of windows
    WindowStartTimes = StartTime_window{i};
    WindowEndTimes = EndTime_window{i};
    N_Windows = length(WindowStartTimes);
    
    %% initialize lists of values 
    N_Anemometers = length(AnemometerName_profile{i}); %number of anemometers
    zsLog_all{i} = zeros(N_Windows,1)*NaN; %roughness height from wind profile
    sigma_zsLog_all{i} = zeros(N_Windows,1)*NaN; %roughness height from wind profile - uncertainty
    ustLog_all{i} = zeros(N_Windows,1)*NaN; %shear velocity from wind profile
    sigma_ustLog_all{i} = zeros(N_Windows,1)*NaN; %shear velocity from wind profile - uncertainty
    tauLog_all{i} = zeros(N_Windows,1)*NaN; %shear stress from wind profile 
    sigma_tauLog_all{i} = zeros(N_Windows,1)*NaN; %shear stress from wind profile - uncertainty

    %% go through time blocks
    for j = 1:N_Windows

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]

        %get specific start and end time
        StartTime = WindowStartTimes(j);
        EndTime = WindowEndTimes(j);

        %get duration of interval in seconds
        T_interval = seconds(EndTime-StartTime);
        
        %get anemometer heights and wind profile
        zU_profile = zU_profile_window{i}(j,:);
        u_profile = u_profile_window{i}(j,:); %uninterpolated wind

        %get indices of values in height range of interest
        ind_fit = find(zU_profile<=zU_max);

        %use only values in height range for fit
        [ustLog,zsLog,sigma_ustLog,sigma_zsLog] = uz_profilefit(u_profile(ind_fit), zU_profile(ind_fit));

        %add values to lists
        zsLog_all{i}(j) = zsLog; %roughness height from wind profile
        sigma_zsLog_all{i}(j) = sigma_zsLog; %roughness height from wind profile - uncertainty
        ustLog_all{i}(j) = ustLog; %shear velocity from wind profile
        sigma_ustLog_all{i}(j) = sigma_ustLog; %shear velocity from wind profile - uncertainty
        tauLog_all{i}(j) = rho_a(i)*ustLog.^2; %shear stress from wind profile 
        sigma_tauLog_all{i}(j) = 2*rho_a(i)*sigma_ustLog*ustLog; %shear stress from wind profile - uncertainty
    end
end

%RENAME ADDITIONAL VALUES NEEDED FOR FLUX LAW ANALYSIS WITH "ALL" ENDING
zU_profile_all = zU_profile_window; %anemometer profile heights

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','N_Sites','*_all'); %save reduced file in GitHub folder