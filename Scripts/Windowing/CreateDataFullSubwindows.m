%% INITIALIZATION
%initialize
clearvars;
close all;

%% interval times
T_subwindow = [...
%     duration(0,0,0.04),...   %17 hrs
%     duration(0,0,0.08),...   %8 hrs
    duration(0,0,0.12);...   %6 hrs 
%     duration(0,0,0.2),...    %3.3 hrs
    duration(0,0,0.4);...    %1.7 hrs
%     duration(0,0,0.6),...    %1 hr
    duration(0,0,1);...    %40 minutes
    duration(0,0,2);...    %20 minutes
    duration(0,0,3);...    %14 minutes
    duration(0,0,5);...    %8 minutes
    duration(0,0,10);...   %4 minutes
    duration(0,0,15);...   %3 minutes
    duration(0,0,30);...   %2 minutes
    duration(0,1,0);...    %1 minute
    duration(0,2,0);...    %30 seconds
    duration(0,3,0);...    %20 seconds
    duration(0,5,0);...    %10 seconds
    duration(0,10,0);...   %5 seconds
    duration(0,15,0);...   %3 seconds
    duration(0,30,0)       %1 second
    ];
N_T_subwindow = length(T_subwindow); %number of subwindow durations
T_subwindow_s = seconds(T_subwindow); %duration of subwindows

%% other parameter values
T_window = duration(0,30,0);

%set info for plotting
Markers_Field = {'s','d','o'};
Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
PlotFont = 14;

%% LOAD DATA AND FUNCTIONS
%folders for loading data, saving data, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_DataBSNE = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder containing BSNE data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for data output
folder_Functions = '../Functions/'; %folder with functions

%paths for loading and saving data
LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Restricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'DataFullSubwindows_30min_Restricted'); %path for 30 minute data

%load data
load(LoadData_Path); %load window data

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize variable arrays
timeofday_subwindow = cell(N_Sites,1); %time of day array
StartTime_subwindow = cell(N_Sites,1); %start time array
EndTime_subwindow = cell(N_Sites,1); %end time array
ind_window_subwindow = cell(N_Sites,1); %index of corresponding window for subwindow

%% initialize wind arrays
zU_subwindow = cell(N_Sites,1); %anemometer height for subwindow
ubar_subwindow = cell(N_Sites,1); %horizontal winds for subwindow
vbar_subwindow = cell(N_Sites,1); %lateral winds for subwindow
wbar_subwindow = cell(N_Sites,1); %vertical winds for subwindow
theta_subwindow = cell(N_Sites,1); %wind angle for subwindow

%% initialize flux arrays
zq_BSNE_subwindow = cell(N_Sites,1); %associated BSNE zq's
zW_subwindow = cell(N_Sites,1); %flux heights for subwindow
sigma_zW_subwindow = cell(N_Sites,1); %uncertainty on flux heights for subwindow
qbar_subwindow = cell(N_Sites,1); %partial flux for subwindow
sigma_qbar_subwindow = cell(N_Sites,1); %partial flux uncertainty for subwindow
nbar_subwindow = cell(N_Sites,1); %particle counts for subwindow
sigma_nbar_subwindow = cell(N_Sites,1); %particle counts uncertainty for subwindow
Cqnbar_subwindow = cell(N_Sites,1); %calibration value for subwindow
sigma_Cqnbar_subwindow = cell(N_Sites,1); %calibration value uncertainty for subwindow

%% GO THROUGH SITES
for i = 1:N_Sites

    %% load BSNE data
    load(strcat(folder_DataBSNE,'FluxBSNE_',Sites{i})); 
    
    %% get information about windows and subwindows
    ind_windows = find(hasfluxdata_window{i}==1); %get indices of subwindows with flux
    N_windows = length(ind_windows); %get number of windows
    
    %% initialize variable arrays - times
    timeofday_subwindow{i} = cell(N_T_subwindow,1); %time of day array
    StartTime_subwindow{i} = cell(N_T_subwindow,1); %start time array
    EndTime_subwindow{i} = cell(N_T_subwindow,1); %end time array
    ind_window_subwindow{i} = cell(N_T_subwindow,1); %index of corresponding window for subwindow
    
    %% initialize variable arrays - wind
    zU_subwindow{i} = cell(N_T_subwindow,1); %anemometer height for subwindow
    ubar_subwindow{i} = cell(N_T_subwindow,1); %wind timeseries array
    vbar_subwindow{i} = cell(N_T_subwindow,1); %lateral winds for subwindow
    wbar_subwindow{i} = cell(N_T_subwindow,1); %vertical wind timeseries array
    theta_subwindow{i} = cell(N_T_subwindow,1); %wind angle for subwindow

    %% initialize variable arrays - flux
    zq_BSNE_subwindow{i} = cell(N_T_subwindow,1); %associated BSNE zq's
    zW_subwindow{i} = cell(N_T_subwindow,1); %flux heights for subwindow
    sigma_zW_subwindow{i} = cell(N_T_subwindow,1); %uncertainty on flux heights for subwindow
    qbar_subwindow{i} = cell(N_T_subwindow,1); %partial flux for subwindow
    sigma_qbar_subwindow{i} = cell(N_T_subwindow,1); %partial flux uncertainty for subwindow
    nbar_subwindow{i} = cell(N_T_subwindow,1); %particle counts for subwindow
    sigma_nbar_subwindow{i} = cell(N_T_subwindow,1); %particle counts uncertainty for subwindow
    Cqnbar_subwindow{i} = cell(N_T_subwindow,1); %calibration value for subwindow
    sigma_Cqnbar_subwindow{i} = cell(N_T_subwindow,1); %calibration value uncertainty for subwindow
    
    %% GO THROUGH MEASUREMENT INTERVALS
    for m = 1:N_T_subwindow

        %% display processing status
        processing_status = [Sites{i},', ',int2str(m),' of ',int2str(N_T_subwindow),', ',datestr(now)]
        
        %% get information about subwindows
        N_subwindows_per_window = T_window/T_subwindow(m);
        N_subwindows = N_windows*N_subwindows_per_window; %get number of subwindows
        
        %% initialize subwindow values - time info
        timeofday_subwindow{i}{m} = zeros(N_subwindows,1); %time of day array
        StartTime_subwindow{i}{m} = datetime(zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1)); %start time array
        EndTime_subwindow{i}{m} = datetime(zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1)); %end time array            
        ind_window_subwindow{i}{m} = zeros(N_subwindows,1); %array of associated window numbers
        
        %% initialize wind values
        zU_subwindow{i}{m} = zeros(N_subwindows,1); %anemometer height
        ubar_subwindow{i}{m} = zeros(N_subwindows,1); %wind timeseries array
        vbar_subwindow{i}{m} = zeros(N_subwindows,1); %transverse wind timeseries array
        wbar_subwindow{i}{m} = zeros(N_subwindows,1); %vertical wind timeseries array
        theta_subwindow{i}{m} = zeros(N_subwindows,1); %wind angle for subwindow

        %% initialize flux values
        zq_BSNE_subwindow{i}{m} = zeros(N_subwindows,1); %associated BSNE zq's
        zW_subwindow{i}{m} = cell(N_subwindows,1); %flux heights for subwindow
        sigma_zW_subwindow{i}{m} = cell(N_subwindows,1); %uncertainty on flux heights for subwindow
        qbar_subwindow{i}{m} = cell(N_subwindows,1); %partial flux for subwindow
        sigma_qbar_subwindow{i}{m} = cell(N_subwindows,1); %partial flux uncertainty for subwindow
        nbar_subwindow{i}{m} = cell(N_subwindows,1); %particle counts for subwindow
        sigma_nbar_subwindow{i}{m} = cell(N_T_subwindow,1); %particle counts uncertainty for subwindow
        Cqnbar_subwindow{i}{m} = cell(N_subwindows,1); %calibration value for subwindow
        sigma_Cqnbar_subwindow{i}{m} = cell(N_subwindows,1); %calibration value uncertainty for subwindow
        
        %% GO THROUGH 30-MINUTE WINDOWS
        for j = 1:N_windows

            %% display processing status
            processing_status = [int2str(j),' of ',int2str(N_windows),', ',datestr(now)]
            
            %% get information about BSNE
            ind_BSNE = find([FluxBSNE.StartTime] <= StartTime_window{i}(ind_windows(j)) & [FluxBSNE.EndTime] >= EndTime_window{i}(ind_windows(j)));
            zq_BSNE = FluxBSNE(ind_BSNE).z.zq; %get best fit saltation layer height for BSNE profile
            
            %% get wind values for window
            zU = zU_base_window{i}(ind_windows(j)); %anemometer height
            t_wind = t_wind_int_window{i}{ind_windows(j)}; %times for interpolated wind
            u = u_int_window{i}{ind_windows(j)}; %rotated interpolated streamwise velocities
            v = v_int_window{i}{ind_windows(j)}; %rotated interpolated lateral winds
            w = w_int_window{i}{ind_windows(j)}; %rotated interpolated vertical velocities

            %% get flux values for window
            t_flux = t_flux_int_window{i}{ind_windows(j)}; %times for interpolated flux
            zW = zW_window{i}{ind_windows(j)}; %get Wenglor heights
            sigma_zW = sigma_zW_window{i}{ind_windows(j)}; %get uncertainty in Wenglor heights
            q = q_int_window{i}{ind_windows(j)}; %get partial fluxes
            Cqn = Cqn_int_window{i}{ind_windows(j)}; %get calibration coefficients
            sigma_Cqn = sigma_Cqn_int_window{i}{ind_windows(j)}; %get uncertainties in calibration coefficients
            n = n_int_window{i}{ind_windows(j)}; %get counts rates
            
            %% get subwindow times within window
            window_StartTime_subwindow = StartTime_window{i}(ind_windows(j))+((1:N_subwindows_per_window)-1)*T_subwindow(m);
            window_EndTime_subwindow = EndTime_window{i}(ind_windows(j))+(1:N_subwindows_per_window)*T_subwindow(m);

            %% go through subwindows
            for k = 1:N_subwindows_per_window
                
                %% get index for adding values to array
                subwindow_array_ind = (j-1)*N_subwindows_per_window+k;

                %% get indices of all subwindow points within window
                subwindow_wind_ind = find(t_wind>=window_StartTime_subwindow(k)&t_wind<window_EndTime_subwindow(k));
                subwindow_flux_ind = find(t_flux>=window_StartTime_subwindow(k)&t_flux<window_EndTime_subwindow(k));

                %% get time of day, start time, end time, and other basic info about subwindow
                StartTime_subwindow{i}{m}(subwindow_array_ind) = window_StartTime_subwindow(k); %start time
                EndTime_subwindow{i}{m}(subwindow_array_ind) = window_EndTime_subwindow(k); %start time
                ind_window_subwindow{i}{m}(subwindow_array_ind) = ind_windows(j); %window corresponding to subwindow
                zU_subwindow{i}{m}(subwindow_array_ind) = zU; %anemometer height - from 30-minute window

                %% get wind data in subwindow
                ubar_subwindow{i}{m}(subwindow_array_ind) = mean(u(subwindow_wind_ind)); %winds for subwindow
                vbar_subwindow{i}{m}(subwindow_array_ind) = mean(v(subwindow_wind_ind)); %lateral winds for subwindow
                wbar_subwindow{i}{m}(subwindow_array_ind) = mean(w(subwindow_wind_ind)); %vertical winds for subwindow
                theta_subwindow{i}{m}(subwindow_array_ind) = 180/pi*atan(mean(v(subwindow_wind_ind))./mean(u((subwindow_wind_ind))));
                timeofday_subwindow{i}{m}(subwindow_array_ind) = hour(window_StartTime_subwindow(k))+minute(window_StartTime_subwindow(k))/60+second(window_StartTime_subwindow(k))/3600;        
                
                %% get flux data in subwindow
                q_subwindow = q(subwindow_flux_ind,:); %get partial fluxes
                n_subwindow = n(subwindow_flux_ind,:); %get counts rates
                Cqn_subwindow = Cqn(subwindow_flux_ind,:); %get calibration coefficient
                sigma_Cqn_subwindow = sigma_Cqn(subwindow_flux_ind,:); %get calibration coefficient uncertainty
                
                %% get mean q profile information for fitting
                qbar = mean(q_subwindow); %compute mean qz profile
                nbar = mean(n_subwindow)/dt_flux_window(i); %compute mean particle counts per second for each qz
                sigma_nbar = sqrt(nbar/T_subwindow_s(m)); %compute uncertainty in mean particle counts
                Cqnbar = mean(Cqn_subwindow); %compute mean calibration coefficient for each qz
                sigma_Cqnbar = mean(sigma_Cqn_subwindow); %compute mean calibration coefficient uncertainty for each qz
                sigma_qbar = sqrt((sigma_Cqnbar.*nbar).^2+(sigma_nbar.*Cqnbar).^2); %compute uncertainty in qz profile, from counts and calibration coefficient
                
                %% add flux values to arrays
                zq_BSNE_subwindow{i}{m}(subwindow_array_ind) = zq_BSNE; %zq for BSNE
                zW_subwindow{i}{m}{subwindow_array_ind} = zW; %flux heights for subwindow
                sigma_zW_subwindow{i}{m}{subwindow_array_ind} = sigma_zW; %uncertainty on flux heights for subwindow
                qbar_subwindow{i}{m}{subwindow_array_ind} = qbar; %partial fluxes
                sigma_qbar_subwindow{i}{m}{subwindow_array_ind} = sigma_qbar; %partial flux uncertainty for subwindow
                nbar_subwindow{i}{m}{subwindow_array_ind} = nbar; %counts rates
                sigma_nbar_subwindow{i}{m}{subwindow_array_ind} = sigma_nbar; %particle counts uncertainty for subwindow
                Cqnbar_subwindow{i}{m}{subwindow_array_ind} = Cqnbar; %calibration coefficient
                sigma_Cqnbar_subwindow{i}{m}{subwindow_array_ind} = sigma_Cqnbar; %calibration coefficient
            end           
        end
        
        %% get rid of extraneous variables
        clear q_subwindow;
        clear n_subwindow;
        clear Cqn_subwindow;
        clear sigma_Cqn_subwindow;
        
        %% save subset of data for measurement interval
        timeofday_subset = timeofday_subwindow{i}{m}; %time of day array
        StartTime_subset = StartTime_subwindow{i}{m}; %start time array
        EndTime_subset = EndTime_subwindow{i}{m}; %end time array            
        ind_window_subset = ind_window_subwindow{i}{m}; %array of associated window numbers
        zU_subset = zU_subwindow{i}{m}; %anemometer height
        ubar_subset = ubar_subwindow{i}{m}; %wind timeseries array
        vbar_subset = vbar_subwindow{i}{m}; %transverse wind timeseries array
        wbar_subset = wbar_subwindow{i}{m}; %vertical wind timeseries array
        theta_subset = theta_subwindow{i}{m}; %wind angle for subwindow
        zq_BSNE_subset = zq_BSNE_subwindow{i}{m}; %associated BSNE zq's
        zW_subset = zW_subwindow{i}{m}; %flux heights for subwindow
        sigma_zW_subset = sigma_zW_subwindow{i}{m}; %uncertainty on flux heights for subwindow
        qbar_subset = qbar_subwindow{i}{m}; %partial flux for subwindow
        sigma_qbar_subset = sigma_qbar_subwindow{i}{m}; %partial flux uncertainty for subwindow
        nbar_subset = nbar_subwindow{i}{m}; %particle counts for subwindow
        sigma_nbar_subset = sigma_nbar_subwindow{i}{m}; %particle counts uncertainty for subwindow
        Cqnbar_subset = Cqnbar_subwindow{i}{m}; %calibration value for subwindow
        sigma_Cqnbar_subset = sigma_Cqnbar_subwindow{i}{m}; %calibration value uncertainty for subwindow
        
        SaveData_Path_T = strcat(folder_SaveData,'ProcessingSteps/DataFullSubwindows_30min_Restricted_',Sites{i},'_',num2str(T_subwindow_s(m))); %path for saving data
        save(SaveData_Path_T,'*subset');
    end
    
    %% save subset of data for site
    timeofday_site = timeofday_subwindow{i}; %time of day array
    StartTime_site = StartTime_subwindow{i}; %start time array
    EndTime_site = EndTime_subwindow{i}; %end time array            
    ind_window_site = ind_window_subwindow{i}; %array of associated window numbers
    zU_site = zU_subwindow{i}; %anemometer height
    ubar_site = ubar_subwindow{i}; %wind timeseries array
    vbar_site = vbar_subwindow{i}; %transverse wind timeseries array
    wbar_site = wbar_subwindow{i}; %vertical wind timeseries array
    theta_site = theta_subwindow{i}; %wind angle for site
    zq_BSNE_site = zq_BSNE_subwindow{i}; %associated BSNE zq's
    zW_site = zW_subwindow{i}; %flux heights for subwindow
    sigma_zW_site = sigma_zW_subwindow{i}; %uncertainty on flux heights for subwindow
    qbar_site = qbar_subwindow{i}; %partial flux for subwindow
    sigma_qbar_site = sigma_qbar_subwindow{i}; %partial flux uncertainty for subwindow
    nbar_site = nbar_subwindow{i}; %particle counts for subwindow
    sigma_nbar_site = sigma_nbar_subwindow{i}; %particle counts uncertainty for subwindow
    Cqnbar_site = Cqnbar_subwindow{i}; %calibration value for subwindow
    sigma_Cqnbar_site = sigma_Cqnbar_subwindow{i}; %calibration value uncertainty for subwindow
    
    SaveData_Path_Site = strcat(folder_SaveData,'ProcessingSteps/DataFullSubwindows_30min_Restricted_',Sites{i}); %path for saving data
    save(SaveData_Path_Site,'*site');
end

%% save analysis data
save(SaveData_Path,'SiteNames','Sites','N_Sites','*subwindow'); 