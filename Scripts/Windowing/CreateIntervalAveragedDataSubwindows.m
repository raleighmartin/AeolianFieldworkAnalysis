%% SCRIPT TO CREATE DATA SUBWINDOWS WITH INTERVAL-AVERAGED TIMESERIES
% deltat = interval averaging duration (shorter)
% Deltat = measurement interval duration (longer)

%% measurement interval times
Deltat_all = ...
    [duration(0,0,30),...
    duration(0,1,0),...
    duration(0,2,0),...
    duration(0,5,0),...
    duration(0,10,0)];
N_Deltat = length(Deltat_all); %number of measurement intervals
subwindow_per_window = duration(0,30,0)./Deltat_all; %get number of measurement interval subwindows per window

%% calculate samples per subwindow
N_samples_per_subwindow = zeros(N_Deltat,N_deltat);
for m = 1:N_Deltat
    for s = 1:N_deltat
        N_samples_per_subwindow(m,s) = floor(Deltat_all(m)/deltat_all(s));
    end
end

%% folders for loading and saving data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

% %% paths for loading and saving data - restricted
% LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Restricted'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'DataIntervalAveragedSubwindows_30min_Restricted'); %path for 30 minute data

%% paths for loading and saving data - unrestricted
LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Unrestricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'DataIntervalAveragedSubwindows_30min_Unrestricted'); %path for 30 minute data

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize variable arrays
timeofday_subwindow = cell(N_Sites,1); %time of day array
StartTime_subwindow = cell(N_Sites,1); %start time array
EndTime_subwindow = cell(N_Sites,1); %end time array
ind_window_subwindow = cell(N_Sites,1); %index of corresponding window for subwindow
hasfluxdata_subwindow = cell(N_Sites,1); %array of which points have flux data
zU_subwindow = cell(N_Sites,1); %anemometer height for subwindow
t_wind_subwindow = cell(N_Sites,1); %times for subwindow with wind
ind_noerr_wind_subwindow = cell(N_Sites,1); %indices of times for subwindow with no wind error
u_subwindow = cell(N_Sites,1); %horizontal winds for subwindow
v_subwindow = cell(N_Sites,1); %lateral winds for subwindow
w_subwindow = cell(N_Sites,1); %vertical winds for subwindow
theta_subwindow = cell(N_Sites,1); %wind angle for subwindow

t_subwindow = cell(N_Sites,1); %times for subwindow with flux and wind
ind_noerr_subwindow = cell(N_Sites,1); %indices of times for subwindow with no flux or wind error
ntotal_subwindow = cell(N_Sites,1); %total counts rate per second for subwindow


%% GO THROUGH SITES
for i = 1:N_Sites
    
    %% get information about windows and subwindows
    N_windows = length(StartTime_window{i}); %get number of windows
    N_subwindows = N_windows*subwindow_per_window; %get number of subwindows for each measurement interval

    %% initialize variable arrays - by measurement interval and sampling interval
    timeofday_subwindow{i} = cell(N_Deltat,N_deltat); %time of day array
    StartTime_subwindow{i} = cell(N_Deltat,N_deltat); %start time array
    EndTime_subwindow{i} = cell(N_Deltat,N_deltat); %end time array
    ind_window_subwindow{i} = cell(N_Deltat,N_deltat); %index of corresponding window for subwindow
    hasfluxdata_subwindow{i} = cell(N_Deltat,N_deltat); %array of which points have flux data
    zU_subwindow{i} = cell(N_Deltat,N_deltat); %anemometer height for subwindow
    t_wind_subwindow{i} = cell(N_Deltat,N_deltat); %timeseries times array for wind
    ind_noerr_wind_subwindow{i} = cell(N_Deltat,N_deltat); %indices of times for no error in wind in subwindow
    u_subwindow{i} = cell(N_Deltat,N_deltat); %wind timeseries array
    v_subwindow{i} = cell(N_Deltat,N_deltat); %lateral winds for subwindow
    w_subwindow{i} = cell(N_Deltat,N_deltat); %vertical wind timeseries array
    theta_subwindow{i} = cell(N_Deltat,N_deltat); %wind angle for subwindow

    t_subwindow{i} = cell(N_Deltat,N_deltat); %timeseries times array for flux and wind
    ind_noerr_subwindow{i} = cell(N_Deltat,N_deltat); %indices of times for no error in flux or wind in subwindow
    ntotal_subwindow{i} = cell(N_Deltat,N_deltat); %counts rate timeseries array 
    
    %% initialize subwindow values
    for m = 1:N_Deltat
        for s = 1:N_deltat
            timeofday_subwindow{i}{m,s} = zeros(N_subwindows(m),1); %time of day array
            StartTime_subwindow{i}{m,s} = datetime(zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1)); %start time array
            EndTime_subwindow{i}{m,s} = datetime(zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1)); %end time array            
            ind_window_subwindow{i}{m,s} = zeros(N_subwindows(m),1); %array of which points have flux data (0=no flux data [default], 1=has flux data)
            hasfluxdata_subwindow{i}{m,s} = zeros(N_subwindows(m),1); %index of corresponding window for subwindow
            zU_subwindow{i}{m,s} = zeros(N_subwindows(m),1); %anemometer height
            t_wind_subwindow{i}{m,s} = cell(N_subwindows(m),1); %timeseries times array
            ind_noerr_wind_subwindow{i}{m,s} = cell(N_subwindows(m),1); %indices of times for no error
            u_subwindow{i}{m,s} = cell(N_subwindows(m),1); %wind timeseries array
            v_subwindow{i}{m,s} = cell(N_subwindows(m),1); %transverse wind timeseries array
            w_subwindow{i}{m,s} = cell(N_subwindows(m),1); %vertical wind timeseries array
            theta_subwindow{i}{m,s} = zeros(N_subwindows(m),1); %wind angle for subwindow
            
            t_subwindow{i}{m,s} = cell(N_subwindows(m),1); %timeseries times array
            ind_noerr_subwindow{i}{m,s} = cell(N_subwindows(m),1); %indices of times for no error
            ntotal_subwindow{i}{m,s} = cell(N_subwindows(m),1); %counts timeseries array
         end
    end
    
    %% GO THROUGH 30-MINUTE WINDOWS
    for j = 1:N_windows

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_windows),', ',datestr(now)]
        
        %% get wind values for window
        hasfluxdata = hasfluxdata_window{i}(j); %determine if window has flux data
        zU = zU_base_window{i}(j); %anemometer height
        t_wind = t_wind_int_window{i}{j}; %times for interpolated wind
        ind_wind_err = ind_wind_err_window{i}{j}; %list of error time indices for wind
        u = u_int_window{i}{j}; %rotated interpolated streamwise velocities
        v = v_int_window{i}{j}; %rotated interpolated lateral winds
        w = w_int_window{i}{j}; %rotated interpolated vertical velocities
        
        %% create timeseries with 0's and 1's for error points
        wind_err_binary = zeros(length(u),1); %initialize timeseries with 0's and 1's for error points
        wind_err_binary(ind_wind_err) = 1; %set error indices to 1

        %% GO THROUGH SAMPLING INTERVALS
        for s = 1:N_deltat

            %% compute window averages
            [u_deltat, t_deltat_u] = window_average(u, t_wind, deltat_all(s)); %wind timeseries
            [v_deltat, ~] = window_average(v, t_wind, deltat_all(s)); %vertical wind timeseries
            [w_deltat, ~] = window_average(w, t_wind, deltat_all(s)); %vertical wind timeseries
            wind_err_binary_deltat = window_average(wind_err_binary, t_wind, deltat_all(s)); %binary wind errors

            %% go through measurement interval durations
            for m = 1:N_Deltat

                %% get subwindow times for measurement interval
                StartTimes_Deltat = StartTime_window{i}(j)+((1:subwindow_per_window(m))-1)*Deltat_all(m);
                EndTimes_Deltat = StartTime_window{i}(j)+(1:subwindow_per_window(m))*Deltat_all(m);

                %% go through subwindows
                for k = 1:subwindow_per_window(m)

                    %get indices for adding values to array
                    subwindow_array_ind = (j-1)*subwindow_per_window(m)+k;

                    %get indices of all subwindow points within window
                    subwindow_window_ind = find(t_deltat_u>=StartTimes_Deltat(k)&t_deltat_u<EndTimes_Deltat(k));
                    
                    %get time of day, start time, end time, and other basic info about subwindow
                    StartTime_subwindow{i}{m,s}(subwindow_array_ind) = StartTimes_Deltat(k); %start time
                    EndTime_subwindow{i}{m,s}(subwindow_array_ind) = EndTimes_Deltat(k); %start time
                    ind_window_subwindow{i}{m,s}(subwindow_array_ind) = j; %window corresponding to subwindow
                    hasfluxdata_subwindow{i}{m,s}(subwindow_array_ind) = hasfluxdata; %indicate in array that this subwindow does contain flux data
                    zU_subwindow{i}{m,s}(subwindow_array_ind) = zU; %anemometer height - from 30-minute window
                     
                    %get all data in subwindow - including error points
                    t_wind_subwindow{i}{m,s}{subwindow_array_ind} = t_deltat_u(subwindow_window_ind); %times for subwindow
                    u_subwindow{i}{m,s}{subwindow_array_ind} = u_deltat(subwindow_window_ind); %winds for subwindow
                    v_subwindow{i}{m,s}{subwindow_array_ind} = v_deltat(subwindow_window_ind); %lateral winds for subwindow
                    w_subwindow{i}{m,s}{subwindow_array_ind} = w_deltat(subwindow_window_ind); %vertical winds for subwindow
                    theta_subwindow{i}{m,s}(subwindow_array_ind) = 180/pi*atan(mean(v_deltat(subwindow_window_ind))./mean(u_deltat((subwindow_window_ind))));
                    timeofday_subwindow{i}{m,s}(subwindow_array_ind) = hour(mean(t_wind_subwindow{i}{m,s}{subwindow_array_ind}))+minute(mean(t_wind_subwindow{i}{m,s}{subwindow_array_ind}))/60+second(mean(t_wind_subwindow{i}{m,s}{subwindow_array_ind}))/3600;
                    
                    %get indices of error points
                    ind_noerr_wind_subwindow{i}{m,s}{subwindow_array_ind} = find(wind_err_binary_deltat(subwindow_window_ind)==0);
                end
            end
        end
       
        %% get flux values for window - only perform calculation if there is flux data in window
        if hasfluxdata == 1
          
            %% get flux values for window
            ntotal = ntotal_int_window{i}{j}; %interpolated total counts rate (counts/s)
            t_flux = t_flux_int_window{i}{j}; %times for interpolated flux
            ind_flux_err = ind_flux_err_window{i}{j}; %list of error time indices for flux
            flux_err_binary = zeros(length(ntotal),1); %initialize timeseries with 0's and 1's for error points
            flux_err_binary(ind_flux_err) = 1; %set error indices to 1
        
            %% GO THROUGH SAMPLING INTERVALS
            for s = 1:N_deltat

                %% compute window averages
                [ntotal_deltat, t_deltat_n] = window_average(ntotal, t_flux, deltat_all(s)); %counts rate timeseries (counts/s)
                [u_deltat, t_deltat_u] = window_average(u, t_wind, deltat_all(s)); %wind timeseries
                flux_err_binary_deltat = window_average(flux_err_binary, t_flux, deltat_all(s)); %binary flux errors
                wind_err_binary_deltat = window_average(wind_err_binary, t_wind, deltat_all(s)); %binary wind errors

                %% get common t, reduce window average timeseries based on these values
                [t_deltat,ind_n,ind_u] = intersect(t_deltat_n,t_deltat_u);
                ntotal_deltat = ntotal_deltat(ind_n);
                flux_err_binary_deltat = flux_err_binary_deltat(ind_n);
                wind_err_binary_deltat = wind_err_binary_deltat(ind_u);
                err_binary_deltat = flux_err_binary_deltat+wind_err_binary_deltat; %combine window averages
                
                %% go through measurement interval durations
                for m = 1:N_Deltat

                    %% get subwindow times for measurement interval
                    StartTimes_Deltat = StartTime_window{i}(j)+((1:subwindow_per_window(m))-1)*Deltat_all(m);
                    EndTimes_Deltat = StartTime_window{i}(j)+(1:subwindow_per_window(m))*Deltat_all(m);

                    %% go through subwindows
                    for k = 1:subwindow_per_window(m)
                        
                        %get indices for adding values to array
                        subwindow_array_ind = (j-1)*subwindow_per_window(m)+k;
                        
                        %get all data in subwindow - including error points
                        subwindow_window_ind = find(t_deltat>=StartTimes_Deltat(k)&t_deltat<EndTimes_Deltat(k)); %get indices of all subwindow points within window
                        t_subwindow{i}{m,s}{subwindow_array_ind} = t_deltat(subwindow_window_ind); %times for subwindow
                        ntotal_subwindow{i}{m,s}{subwindow_array_ind} = ntotal_deltat(subwindow_window_ind); %counts for subwindow
                       
                        %get indices of error points
                        ind_noerr_subwindow{i}{m,s}{subwindow_array_ind} = find(err_binary_deltat(subwindow_window_ind)==0);
                    end
                end
            end
        end
    end
end

%% save analysis data
save(SaveData_Path,'SiteNames','Sites','N_Sites','*subwindow','*all'); 