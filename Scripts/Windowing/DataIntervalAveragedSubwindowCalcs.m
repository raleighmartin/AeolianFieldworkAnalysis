%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
kappa = 0.4; %von Karman parameter
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
z0 = [6.87e-5, 1.22e-4, 9.96e-5];
fD_min = 0.005; %minimum detection rate for zero flux

%% folders for loading and saving data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

% %% paths for loading and saving data - restricted
% LoadData_Path = strcat(folder_LoadData,'DataSubwindows_30min_Restricted'); %path for 30 minute data - for thresholds analysis
% SaveData_Path = strcat(folder_SaveData,'DataSubwindowCalcs_30min_Restricted'); %path for 30 minute data - for thresholds analysis
% 
%% paths for loading and saving data - unrestricted
LoadData_Path = strcat(folder_LoadData,'DataIntervalAveragedSubwindows_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis
SaveData_Path = strcat(folder_SaveData,'DataIntervalAveragedSubwindowCalcs_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%% get info about measurement and sampling intervals
N_Deltat = length(Deltat_all); %number of measurement intervals
N_deltat = length(deltat_all); %number of sampling intervals

%% initialize variable lists

%initialize lists of wind values for lowest anemometer
ubar_subwindow_all = cell(N_Sites,1); %mean wind velocity from lowest anemometer
vbar_subwindow_all = cell(N_Sites,1); %mean lateral wind velocity from lowest anemometer
wbar_subwindow_all = cell(N_Sites,1); %mean vertical wind velocity from lowest anemometer

%initialize lists of flux values
fD_subwindow_all = cell(N_Sites,1); %Wenglor detection frequency list
fQ_subwindow_all = cell(N_Sites,1); %Wenglor transport frequency matrix
lambda_subwindow_all = cell(N_Sites,1); %Wenglor arrival rate per sampling interval for flux frequency correction
uth_subwindow_all = cell(N_Sites,1); %threshold wind speed from flux frequency
tauth_subwindow_all = cell(N_Sites,1); %threshold wind stress from flux frequency
ustth_subwindow_all = cell(N_Sites,1); %threshold wind shear from flux frequency

%% GO THROUGH SITES
for i = 1:N_Sites
    
    %% initialize lists of values

    %wind values - lowest anemometer
    ubar_subwindow_all{i} = cell(N_Deltat,N_deltat); %observed wind velocity from lowest anemometer
    vbar_subwindow_all{i} = cell(N_Deltat,N_deltat); %mean transverse wind
    wbar_subwindow_all{i} = cell(N_Deltat,N_deltat); %mean vertical wind
    
    %flux values
    fD_subwindow_all{i} = cell(N_Deltat,N_deltat); %Wenglor detection frequency - 1 s
    fQ_subwindow_all{i} = cell(N_Deltat,N_deltat); %Wenglor total transport frequency - 1 s
    lambda_subwindow_all{i} = cell(N_Deltat,N_deltat); %Wenglor particle arrival rate - 1 s
    uth_subwindow_all{i} = cell(N_Deltat,N_deltat); %threshold wind speed from flux frequency
    tauth_subwindow_all{i} = cell(N_Deltat,N_deltat); %threshold wind stress from flux frequency
    ustth_subwindow_all{i} = cell(N_Deltat,N_deltat); %threshold wind shear from flux frequency
      
    %% go through measurement intervals
    for m = 1:N_Deltat
            
        %% go through sampling intervals
        for s = 1:N_deltat

            %% initialize lists of values
            
            %get number of subwindows
            N_subwindows = length(StartTime_subwindow{i}{m,s});
            
            %wind values - lowest anemometer
            ubar_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %observed wind velocity from lowest anemometer
            vbar_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %observed lateral wind velocity from lowest anemometer
            wbar_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %observed vertical wind velocity from lowest anemometer
            
            %flux values
            fD_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %Wenglor detection frequency - 1 s
            fQ_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %Wenglor total transport frequency - 1 s
            lambda_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %Wenglor particle arrival rate - 1 s
            uth_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %threshold wind speed from flux frequency
            tauth_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %threshold wind stress from flux frequency
            ustth_subwindow_all{i}{m,s} = zeros(N_subwindows,1)*NaN; %threshold wind shear from flux frequency
                        
            %% go through subwindows, perform calculations
            for j = 1:N_subwindows

                %get mean velocity value, add to list
                ubar_subwindow_all{i}{m,s}(j) = mean(u_subwindow{i}{m,s}{j});
                vbar_subwindow_all{i}{m,s}(j) = mean(v_subwindow{i}{m,s}{j});
                wbar_subwindow_all{i}{m,s}(j) = mean(w_subwindow{i}{m,s}{j});
               
                %only perform threshold calculations if there is flux data in subwindow
                if hasfluxdata_subwindow{i}{m,s}(j)==1
                    
                    %get non-error points, times, and counts rates
                    ind_noerr = ind_noerr_subwindow{i}{m,s}{j}; %indices without error
                    ntotal = ntotal_subwindow{i}{m,s}{j}(ind_noerr); %particle counts
                    u = u_subwindow{i}{m,s}{j}(ind_noerr); %horizontal wind speed
                    T = length(ntotal); %get number of time steps
                    
                    %calculate flux activity
                    [fD,fQ,lambda] = CalculateFluxActivity(ntotal,seconds(deltat_all(s)),fD_min);
                    
                    %determine wind corresponding to fQ (only non-error points)
                    u_sort = sort(u); %sort u's
                    ind_uth_fQ = round((1-fQ)*T); %get index in list of u's corresponding to threshold
                    if (ind_uth_fQ==0||ind_uth_fQ==T)||isnan(ind_uth_fQ) %points don't count if fQ = 0, 1, or is undefined
                        uth = NaN;
                        ustth = NaN;
                        tauth = NaN;
                    else
                        uth = u_sort(ind_uth_fQ); %threshold wind speed
                        ustth = (kappa*uth)/log(zU_subwindow{i}{m,s}(j)/z0(i)); %threshold shear velocity
                        tauth = rho_a(i)*ustth^2; %threshold shear stress
                    end

                    %add to arrays of values
                    fD_subwindow_all{i}{m,s}(j) = fD; %detection frequency
                    fQ_subwindow_all{i}{m,s}(j) = fQ; %flux frequency
                    lambda_subwindow_all{i}{m,s}(j) = lambda; %estimated particle detection rate
                    uth_subwindow_all{i}{m,s}(j) = uth; %threshold wind speed from flux frequency
                    tauth_subwindow_all{i}{m,s}(j) = tauth; %threshold wind stress from flux frequency
                    ustth_subwindow_all{i}{m,s}(j) = ustth; %threshold wind shear from flux frequency
                end
            end
        end
    end
end

%%RENAME ADDITIONAL VALUES NEEDED FOR FLUX LAW ANALYSIS
StartTime_subwindow_all = StartTime_subwindow; %start times
EndTime_subwindow_all = EndTime_subwindow; %start times
timeofday_subwindow_all = timeofday_subwindow; %time of day
zU_subwindow_all = zU_subwindow; %anemometer heights
theta_subwindow_all = theta_subwindow; %wind direction
ind_window_subwindow_all = ind_window_subwindow; %index of associated window

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','N_Sites','*all'); %save reduced file in GitHub folder