%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% subwindow times
delta_t_subwindow_all = [duration(0,0,25),duration(0,1,0),duration(0,3,0),duration(0,5,0),duration(0,10,0),duration(0,30,0)];
delta_t_subwindow_all = [duration(0,0,25),duration(0,1,0),duration(0,3,0),duration(0,5,0),duration(0,10,0)];
delta_t_subwindow_all = fliplr(delta_t_subwindow_all);
N_delta_t_subwindow = length(delta_t_subwindow_all);

%% parameter values
delta_t_avg_analysis = duration(0,0,1); %time increment to use for primary window-average analysis
delta_t_lowpass_analysis = duration(0,0,25); %time increment to use for primary lowpass analysis
% u_ft = [9, 7.9, 7.3]; %estimated fluid threshold wind for each site
% u_it = [7.2, 6.3, 5.8]; %estimated impact threshold wind for each site
tau_ft = [0.188, 0.147, 0.129]; %estimated fluid threshold stress for each site (Pa)
tau_it = [0.128, 0.095, 0.084]; %estimated impact threshold stress for each site (Pa)
rho_a = 1.18; %air density kg/m^3
g = 9.8; %gravity m/s^2
kappa = 0.4; %von Karman parameter
z0 = 1e-4; %aerodynamic roughness length (m)
fD_min = 0.005; %minimum detection rate for zero flux

%% information about where to load/save data, plots, and functions
folder_LoadData_1 = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
folder_LoadData_2 = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/Thresholds/Windows/'; %folder for plots
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for loading window data
LoadData_1_Path = strcat(folder_LoadData_1,'TimeWindows_30min.mat'); %path for loading time window data
LoadData_2_Path = strcat(folder_LoadData_2,'WindowAverageWindows_30min.mat'); %path for loading window-averaged window data
LoadData_3_Path = strcat(folder_LoadData_2,'LowPassWindows_30min.mat'); %path for loading low-pass window data
%LoadData_4_Path = strcat(folder_LoadData_2,'DataWindows_30min.mat'); %path for loading full data windows

%% load window data
load(LoadData_1_Path);
load(LoadData_2_Path);
load(LoadData_3_Path);
%load(LoadData_4_Path);

%% record values so they don't get deleted in loop
save_t_avg_window_save = t_avg_window;
save_t_wind_lowpass_window_save = t_wind_lowpass_window;
save_t_flux_lowpass_window_save = t_flux_lowpass_window;
save_u_avg_window_save = u_avg_window;
save_u_lowpass_window_save = u_lowpass_window;
save_n_avg_window_save = n_avg_window;
save_n_lowpass_window_save = n_lowpass_window;

%% go through all subwindow times
for m = 1:N_delta_t_subwindow
    
    % determine specific subwindow time
    delta_t_subwindow = delta_t_subwindow_all(m)

    % path for saving output data
    SaveData_Path = strcat(folder_SaveData,'ActivityWindows_30min_',int2str(seconds(delta_t_subwindow)),'s.mat'); 

    % %info for demo plots
    % t_plot_s = 150; %duration of plot
    % N_subwindow_plot = t_plot_s/seconds(delta_t_subwindow); %number of subwindows for plotting
    % t_subwindow_plot = ((1:N_subwindow_plot)-0.5).*seconds(delta_t_subwindow); %times for subwindow plotting

    %% get time information for analysis
    ind_delta_t_avg_analysis = find(delta_t_avg_window==delta_t_avg_analysis); %index of window avg delta_t for primary analysis
    ind_delta_t_lowpass_analysis = find(delta_t_lowpass_window==delta_t_lowpass_analysis); %index of lowpass delta_t for primary analysis
    N_subwindows = WindowTimeInterval/delta_t_subwindow; %number of subwindows per window for analysis
    N_delta_t_avg = length(delta_t_avg_window);

    %% initialize variable list

    %transport detection / activity
    fD_avg_window = cell(N_Sites,1); %Wenglor detection frequency list for each window
    fQ_avg_window = cell(N_Sites,1); %Wenglor transport frequency matrix for each window
    lambda_avg_window = cell(N_Sites,1); %estimated particle arrival rate per averaging window

    %values based on window-averaged wind
    uth_fQ_avg_window = cell(N_Sites,1); %threshold wind estimated from histogram of winds and fQ
    ustth_fQ_avg_window = cell(N_Sites,1); %threshold shear velocity estimated uth and law-of-the-wall 
    tauth_fQ_avg_window = cell(N_Sites,1); %threshold wind stress estimated uth and law-of-the-wall 
    fplus_avg_window = cell(N_Sites,1); %fraction of time with u_avg above fluid threshold
    fminus_avg_window = cell(N_Sites,1); %fraction of time with u_avg below impact threshold
    fint_avg_window = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone
    fint_down_avg_window = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone from above
    fint_up_avg_window = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone from below
    t_ftdown_avg_window = cell(N_Sites,1); %times of fluid threshold downward crossings
    t_itup_avg_window = cell(N_Sites,1); %times of impact threshold upward crossings

    %values based on low-pass wind
    uth_fQ_lowpass_window = cell(N_Sites,1); %threshold wind estimated from histogram of winds and fQ
    ustth_fQ_lowpass_window = cell(N_Sites,1); %threshold shear velocity estimated uth and law-of-the-wall 
    tauth_fQ_lowpass_window = cell(N_Sites,1); %threshold wind stress estimated uth and law-of-the-wall 
    fplus_lowpass_window = cell(N_Sites,1); %fraction of time with u_lowpass above fluid threshold
    fminus_lowpass_window = cell(N_Sites,1); %fraction of time with u_lowpass below impact threshold
    fint_lowpass_window = cell(N_Sites,1); %fraction of time with u_lowpass in intermediate zone
    fint_down_lowpass_window = cell(N_Sites,1); %fraction of time with u_lowpass in intermediate zone from above
    fint_up_lowpass_window = cell(N_Sites,1); %fraction of time with u_lowpass in intermediate zone from below

    %transport detection / activity - subwindows
    fD_avg_subwindow = cell(N_Sites,1); %Wenglor detection frequency list for each window
    fQ_avg_subwindow = cell(N_Sites,1); %Wenglor transport frequency matrix for each window
    timeofday_subwindow = cell(N_Sites,1); %time associated with subwindow, in hour of day

    %values based on window-averaged wind - subwindows
    uth_fQ_avg_subwindow = cell(N_Sites,1); %threshold wind estimated from histogram of winds and fQ
    ustth_fQ_avg_subwindow = cell(N_Sites,1); %threshold shear velocity estimated uth and law-of-the-wall 
    tauth_fQ_avg_subwindow = cell(N_Sites,1); %threshold wind stress estimated uth and law-of-the-wall 
    fplus_avg_subwindow = cell(N_Sites,1); %fraction of time with u_avg above fluid threshold
    fminus_avg_subwindow = cell(N_Sites,1); %fraction of time with u_avg below impact threshold
    fint_avg_subwindow = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone
    fint_down_avg_subwindow = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone from above
    fint_up_avg_subwindow = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone from below

    %values based on low-pass wind - subwindows
    uth_fQ_lowpass_subwindow = cell(N_Sites,1); %threshold wind estimated from histogram of winds and fQ
    ustth_fQ_lowpass_subwindow = cell(N_Sites,1); %threshold shear velocity estimated uth and law-of-the-wall 
    tauth_fQ_lowpass_subwindow = cell(N_Sites,1); %threshold wind stress estimated uth and law-of-the-wall 
    fplus_lowpass_subwindow = cell(N_Sites,1); %fraction of time with u_lowpass above fluid threshold
    fminus_lowpass_subwindow = cell(N_Sites,1); %fraction of time with u_lowpass below impact threshold
    fint_lowpass_subwindow = cell(N_Sites,1); %fraction of time with u_lowpass in intermediate zone
    fint_down_lowpass_subwindow = cell(N_Sites,1); %fraction of time with u_lowpass in intermediate zone from above
    fint_up_lowpass_subwindow = cell(N_Sites,1); %fraction of time with u_lowpass in intermediate zone from below

    % %% INITIALIZE SAMPLE TIMESERIES PLOT
    % %Create sample timeseries plot
    % figure(1); clf;
    % subplot(2,3,1);
    % xlabel('t (s)');
    % ylabel('N (/s)');
    % title('raw');
    % 
    % subplot(2,3,2);
    % xlabel('t (s)');
    % ylabel('N (/s)');
    % title('{\delta}t window avg');
    % 
    % subplot(2,3,3);
    % xlabel('t (s)');
    % ylabel('f_Q');
    % title('{\Delta}t stats');
    % 
    % subplot(2,3,4);
    % xlabel('t (s)');
    % ylabel('u (m/s)');
    % 
    % subplot(2,3,5); hold on;
    % xlabel('t (s)');
    % ylabel('u (m/s)');
    % 
    % subplot(2,3,6);
    % xlabel('t (s)');
    % ylabel('u_{th} (m/s)');

    %% PERFORM ANALYSIS FOR EACH SITE
    for i = 1:N_Sites
        N_windows = length(n_avg_window{i}); %get number of windows

        %% initialize lists of values
        fD_avg_window{i} = zeros(N_windows,N_delta_t_avg); %Wenglor detection frequency list for each window
        fQ_avg_window{i} = zeros(N_windows,N_delta_t_avg); %Wenglor transport frequency matrix for each window
        lambda_avg_window{i} = zeros(N_windows,N_delta_t_avg); %estimated particle arrival rate per averaging interval

        %threshold values based on window-averaged wind
        uth_fQ_avg_window{i} = zeros(N_windows,N_delta_t_avg); %threshold wind estimated from histogram of winds and fQ
        ustth_fQ_avg_window{i} = zeros(N_windows,N_delta_t_avg); %threshold shear velocity estimated uth and law-of-the-wall 
        tauth_fQ_avg_window{i} = zeros(N_windows,N_delta_t_avg); %threshold wind stress estimated uth and law-of-the-wall 

        %wind hysteresis values based on window-averaged wind
        fplus_avg_window{i} = zeros(N_windows,N_delta_t_avg); %fraction of time with u_avg above fluid threshold
        fminus_avg_window{i} = zeros(N_windows,N_delta_t_avg); %fraction of time with u_avg below impact threshold
        fint_avg_window{i} = zeros(N_windows,N_delta_t_avg); %fraction of time with u_avg in intermediate zone
        fint_down_avg_window{i} = zeros(N_windows,N_delta_t_avg); %fraction of time with u_avg in intermediate zone from above
        fint_up_avg_window{i} = zeros(N_windows,N_delta_t_avg); %fraction of time with u_avg in intermediate zone from below
        t_ftdown_avg_window{i} = cell(N_windows,N_delta_t_avg); %times of fluid threshold downward crossings
        t_itup_avg_window{i} = cell(N_windows,N_delta_t_avg); %times of impact threshold upward crossings

        %threshold values based on low-pass wind
        uth_fQ_lowpass_window{i} = zeros(N_windows,N_delta_t_avg); %threshold wind estimated from histogram of winds and fQ
        ustth_fQ_lowpass_window{i} = zeros(N_windows,N_delta_t_avg); %threshold shear velocity estimated uth and law-of-the-wall 
        tauth_fQ_lowpass_window{i} = zeros(N_windows,N_delta_t_avg); %threshold wind stress estimated uth and law-of-the-wall 

        %wind hysteresis values based on low-pass wind
        fplus_lowpass_window{i} = zeros(N_windows,N_lowpass_window); %fraction of time with u_lowpass above fluid threshold
        fminus_lowpass_window{i} = zeros(N_windows,N_lowpass_window); %fraction of time with u_lowpass below impact threshold
        fint_lowpass_window{i} = zeros(N_windows,N_lowpass_window); %fraction of time with u_lowpass in intermediate zone
        fint_down_lowpass_window{i} = zeros(N_windows,N_lowpass_window); %fraction of time with u_lowpass in intermediate zone from above
        fint_up_lowpass_window{i} = zeros(N_windows,N_lowpass_window); %fraction of time with u_lowpass in intermediate zone from below

        %% initialize lists of values - subwindows
        fD_avg_subwindow{i} = zeros(N_windows,N_subwindows); %Wenglor detection frequency list for each window
        fQ_avg_subwindow{i} = zeros(N_windows,N_subwindows); %Wenglor transport frequency matrix for each window
        timeofday_subwindow{i} = zeros(N_windows,N_subwindows); %time associated with subwindow, in hour of day

        %threshold values based on window-averaged wind - subwindows
        uth_fQ_avg_subwindow{i} = zeros(N_windows,N_subwindows); %threshold wind estimated from histogram of winds and fQ
        ustth_fQ_avg_subwindow{i} = zeros(N_windows,N_subwindows); %threshold shear velocity estimated uth and law-of-the-wall 
        tauth_fQ_avg_subwindow{i} = zeros(N_windows,N_subwindows); %threshold wind stress estimated uth and law-of-the-wall 

        %wind hysteresis values based on window-averaged wind - subwindows
        fplus_avg_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_avg above fluid threshold
        fminus_avg_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_avg below impact threshold
        fint_avg_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_avg in intermediate zone
        fint_down_avg_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_avg in intermediate zone from above
        fint_up_avg_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_avg in intermediate zone from below

        %threshold values based on low-pass wind - subwindows
        uth_fQ_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %threshold wind estimated from histogram of winds and fQ
        ustth_fQ_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %threshold shear velocity estimated uth and law-of-the-wall 
        tauth_fQ_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %threshold wind stress estimated uth and law-of-the-wall 

        %wind hysteresis values based on low-pass wind - subwindows
        fplus_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_lowpass above fluid threshold
        fminus_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_lowpass below impact threshold
        fint_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_lowpass in intermediate zone
        fint_down_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_lowpass in intermediate zone from above
        fint_up_lowpass_subwindow{i} = zeros(N_windows,N_subwindows); %fraction of time with u_lowpass in intermediate zone from below

        %% go through time windows
        for j = 1:N_windows

            %display processing status
            processing_status = [SiteNames{i},', ',int2str(j),' of ',int2str(N_windows),', ',datestr(now)]

            %estimate thresholds for window
            u_ft = (sqrt(tau_ft(i)/rho_a)/kappa)*log(zU_window{i}(j)/z0);
            u_it = (sqrt(tau_it(i)/rho_a)/kappa)*log(zU_window{i}(j)/z0);

            %get lowpass wind for analysis - single lowpass time scale
            u_lowpass = u_lowpass_window{i}{j,ind_delta_t_lowpass_analysis}; %lowpass u's
            T_lowpass = length(u_lowpass); %lowpass number of timesteps
            u_lowpass_sort = sort(u_lowpass); %sort lowpass u's

            %% go through averaging times - total window analysis
            for k = 1:N_delta_t_avg

                %extract relevant info
                t_avg = t_avg_window{i}{j,k}; %get times
                n_avg = n_avg_window{i}{j,k}; %window-averaged count rates values
                u_avg = u_avg_window{i}{j,k}; %window-averaged u values
                ind_flux_err = ind_flux_err_avg_window{i}{j,k}; %indices of error points in flux window-average
                ind_wind_err = ind_wind_err_avg_window{i}{j,k}; %indices of error points in wind window-average

                %keep only nonerror points for calculation
                ind_err = union(ind_flux_err,ind_wind_err); %indices of all error points
                ind_noerr = setdiff((1:length(u_avg))',ind_err); %indices of nonerror points
                T_avg = length(ind_noerr); %number of timesteps (replace existing)
                n_avg = n_avg(ind_noerr,:); %counts for nonerror points (replace existing)
                u_avg = u_avg(ind_noerr); %winds for nonerror points (replace existing)

                %evaluate flux detections
                D = sum(n_avg')'>0; %particle detections (1 if detected, 0 if not)

                %calculate flux detection frequency and mean particle counts for no error points
                fD = sum(D)/T_avg; %calculate detection rate (only non-error points)
                N_s = sum(mean(n_avg)); %mean particle counts rate profile

                %estimate particle arrival rate per averaging window
                if fD<=fD_min %set to zero if below detection limit
                    lambda = 0;
                else %otherwise estimate arrival rate based on fD
                    lambda = N_s*seconds(delta_t_avg_window(k))/fD;
                end

                %estimate flux activity from flux detection rate and particle arrival rate
                if lambda==0 %if no (or negligible) flux, set frequencies to zero
                    fQ = 0;
                elseif fD==1 %if fD = 1, set fQ to 1
                    fQ = 1;
                else %otherwise, estimate fQ based on other parameters
                    fQ = fD/(1-exp(-lambda)); %calculate fQ
                end

                %add to lists of all values
                fD_avg_window{i}(j,k) = fD; %detection frequency
                fQ_avg_window{i}(j,k) = fQ; %flux frequency
                lambda_avg_window{i}(j,k) = lambda; %estimated particle detection rate

                %determine window-averaged wind corresponding to fQ
                u_avg_sort = sort(u_avg); %sort u's
                ind_uth_fQ = round((1-fQ)*T_avg); %get index in list of u's corresponding to threshold
                if ind_uth_fQ==0||isnan(ind_uth_fQ)
                    uth_avg_fQ = NaN;
                    ustth_avg_fQ = NaN;
                    tauth_avg_fQ = NaN;
                else
                    uth_avg_fQ = u_avg_sort(ind_uth_fQ); %threshold wind speed
                    ustth_avg_fQ = (kappa*uth_avg_fQ)/log(zU_window{i}(j)/z0); %threshold shear velocity
                    tauth_avg_fQ = rho_a*ustth_avg_fQ^2; %threshold shear stress
                end

                %add to list of values from window-averaged winds
                uth_fQ_avg_window{i}(j,k) = uth_avg_fQ; %threshold wind from flux frequency
                tauth_fQ_avg_window{i}(j,k) = tauth_avg_fQ; %threshold wind from flux frequency
                ustth_fQ_avg_window{i}(j,k) = ustth_avg_fQ; %threshold wind from flux frequency

                %determine fractions of u in different regions - window-averaged u
                [fplus_avg,fminus_avg,fint_avg,fint_up_avg,fint_down_avg] = CalculateWindHysteresis(u_avg,u_ft,u_it);

                %determine times of threshold crossings
                t_ftdown_avg_window{i}{j,k} = t_avg(diff(u_avg>=u_ft)==-1); %times of fluid threshold downward crossings
                t_itup_avg_window{i}{j,k} = t_avg(diff(u_avg<=u_it)==-1); %times of impact threshold upward crossings

                %add to list of values from window-averaged winds
                fplus_avg_window{i}(j,k) = fplus_avg; %fraction of time with u_avg above fluid threshold
                fminus_avg_window{i}(j,k) = fminus_avg; %fraction of time with u_avg below impact threshold
                fint_avg_window{i}(j,k) = fint_avg; %fraction of time with u_avg in intermediate zone
                fint_down_avg_window{i}(j,k) = fint_down_avg; %fraction of time with u_avg in intermediate zone with downward crossing
                fint_up_avg_window{i}(j,k) = fint_up_avg; %fraction of time with u_avg in intermediate zone with upward crossing

                %determine low-pass wind corresponding to fQ
                ind_uth_lowpass_fQ = round((1-fQ)*T_lowpass); %get index in list of u's corresponding to threshold
                if ind_uth_lowpass_fQ==0||isnan(ind_uth_lowpass_fQ)
                    uth_lowpass_fQ = NaN;
                    ustth_lowpass_fQ = NaN;
                    tauth_lowpass_fQ = NaN;
                else
                    uth_lowpass_fQ = u_lowpass_sort(ind_uth_lowpass_fQ); %threshold wind speed
                    ustth_lowpass_fQ = (kappa*uth_lowpass_fQ)/log(zU_window{i}(j)/z0); %threshold shear velocity
                    tauth_lowpass_fQ = rho_a*ustth_lowpass_fQ^2; %threshold shear stress
                end

                %add to list of values from low-pass winds
                uth_fQ_lowpass_window{i}(j,k) = uth_lowpass_fQ; %threshold wind from flux frequency
                tauth_fQ_lowpass_window{i}(j,k) = tauth_lowpass_fQ; %threshold wind from flux frequency
                ustth_fQ_lowpass_window{i}(j,k) = ustth_lowpass_fQ; %threshold wind from flux frequency  
            end

            %% go through low-pass filtering times to investigate their effects on whole window wind partitioning
            for k=1:N_lowpass_window

                %determine fractions of u in different regions - low-pass u
                [fplus_lowpass,fminus_lowpass,fint_lowpass,fint_up_lowpass,fint_down_lowpass] = CalculateWindHysteresis(u_lowpass_window{i}{j,k},u_ft,u_it);

                %add to list of values from low-pass winds
                fplus_lowpass_window{i}(j,k) = fplus_lowpass; %fraction of time with u_lowpass above fluid threshold
                fminus_lowpass_window{i}(j,k) = fminus_lowpass; %fraction of time with u_lowpass below impact threshold
                fint_lowpass_window{i}(j,k) = fint_lowpass; %fraction of time with u_lowpass in intermediate zone
                fint_down_lowpass_window{i}(j,k) = fint_down_lowpass; %fraction of time with u_lowpass in intermediate zone with downward crossing
                fint_up_lowpass_window{i}(j,k) = fint_up_lowpass; %fraction of time with u_lowpass in intermediate zone with upward crossing

            end

            %% go through subwindows

            %create subwindows
            BoundingTimes_subwindow = (StartTime_window{i}(j):delta_t_subwindow:EndTime_window{i}(j))';
            StartTime_subwindow = BoundingTimes_subwindow(1:end-1);
            EndTime_subwindow = BoundingTimes_subwindow(2:end);

            %get window-averaged times, wind, particle counts, and error points for analysis
            t_avg = t_avg_window{i}{j,ind_delta_t_avg_analysis}; %times for window average
            n_avg = n_avg_window{i}{j,ind_delta_t_avg_analysis}; %window-averaged count rates values
            u_avg = u_avg_window{i}{j,ind_delta_t_avg_analysis}; %window-averaged u values
            ind_flux_err = ind_flux_err_avg_window{i}{j,ind_delta_t_avg_analysis}; %indices of error points in flux window-average
            ind_wind_err = ind_wind_err_avg_window{i}{j,ind_delta_t_avg_analysis}; %indices of error points in wind window-average
            ind_err = union(ind_flux_err,ind_wind_err); %indices of all error points
            ind_noerr = setdiff((1:length(t_avg))',ind_err); %indices of nonerror points

            %times of threshold crossings
            t_ftdown_avg = t_ftdown_avg_window{i}{j,ind_delta_t_avg_analysis}; %times of fluid threshold downward crossings
            t_itup_avg = t_itup_avg_window{i}{j,ind_delta_t_avg_analysis}; %times of impact threshold upward crossings

            %get lowpass times and wind
            t_lowpass = t_wind_lowpass_window{i}{j}; %times for lowpass values
            u_lowpass = u_lowpass_window{i}{j,ind_delta_t_lowpass_analysis}; %low-pass u values

            %% go through subwindows
            for k = 1:N_subwindows

                %get window-averaged points in subwindow
                ind_avg_subwindow = find(t_avg>=StartTime_subwindow(k)&t_avg<EndTime_subwindow(k)); %get indices of all subwindow points
                ind_avg_subwindow = intersect(ind_noerr,ind_avg_subwindow); %get indices only of subwindow points with no error (replace existing)
                T_avg_subwindow = length(ind_avg_subwindow); %number of timesteps
                n_avg_subwindow = n_avg(ind_avg_subwindow,:); %counts for subwindow
                u_avg_subwindow = u_avg(ind_avg_subwindow); %winds for subwindow
                timeofday_subwindow{i}(j,k) = hour(mean(t_avg(ind_avg_subwindow))); %time associated with subwindow, in hour of day

                %evaluate flux detections - subwindows
                D = sum(n_avg_subwindow')'>0; %particle detections (1 if detected, 0 if not)

                %calculate flux detection frequency and mean particle counts for no error points  - subwindows
                fD = sum(D)/T_avg_subwindow; %calculate detection rate (only non-error points)
                N_s = sum(mean(n_avg_subwindow)); %mean particle counts rate profile

                %get particle arrival rate per averaging window from whole window value
                lambda = lambda_avg_window{i}(j,ind_delta_t_avg_analysis);

                %estimate flux activity from flux detection rate and particle arrival rate
                if lambda==0 %if no (or negligible) flux, set frequencies to zero
                    fQ = 0;
                elseif fD==1 %if fD = 1, set fQ to 1
                    fQ = 1;
                else %otherwise, estimate fQ based on other parameters
                    fQ = fD/(1-exp(-lambda)); %calculate fQ
                end

                %add to lists of all values - subwindows
                fD_avg_subwindow{i}(j,k) = fD; %detection frequency
                fQ_avg_subwindow{i}(j,k) = fQ; %flux frequency

                %determine window-averaged wind corresponding to fQ - subwindows
                u_avg_sort = sort(u_avg_subwindow); %sort u's
                ind_uth_fQ = round((1-fQ)*T_avg_subwindow); %get index in list of u's corresponding to threshold
                if ind_uth_fQ==0||isnan(ind_uth_fQ)
                    uth_avg_fQ = NaN;
                    ustth_avg_fQ = NaN;
                    tauth_avg_fQ = NaN;
                else
                    uth_avg_fQ = u_avg_sort(ind_uth_fQ); %threshold wind speed
                    ustth_avg_fQ = (kappa*uth_avg_fQ)/log(zU_window{i}(j)/z0); %threshold shear velocity
                    tauth_avg_fQ = rho_a*ustth_avg_fQ^2; %threshold shear stress
                end

                %add to list of values from window-averaged winds - subwindows
                uth_fQ_avg_subwindow{i}(j,k) = uth_avg_fQ; %threshold wind from flux frequency
                tauth_fQ_avg_subwindow{i}(j,k) = tauth_avg_fQ; %threshold wind from flux frequency
                ustth_fQ_avg_subwindow{i}(j,k) = ustth_avg_fQ; %threshold wind from flux frequency

                %get state of last threshold crossing
                t_subwindow_start = t_avg(ind_avg_subwindow(1)); %time for start of subwindow

                %index of last downward crossing
                t_ftdown_last = t_ftdown_avg(find(t_ftdown_avg<=t_subwindow_start, 1, 'last' ));
                if isempty(t_ftdown_last)
                    t_ftdown_last = datetime(0,0,0); %set to zero if none found
                end

                %index of last upward crossing
                ind_itup_last = t_itup_avg(find(t_itup_avg<=t_subwindow_start, 1, 'last' ));
                if isempty(ind_itup_last)
                    ind_itup_last = datetime(0,0,0); %set to zero if none found
                end

                %set init_state based on whether downward crossing or upward crossing was more recent
                if t_ftdown_last>ind_itup_last
                    init_state=1;
                elseif ind_itup_last>t_ftdown_last
                    init_state=-1;
                else
                    init_state=0;
                end

                %determine fractions of u in different regions - window-averaged u - subwindows
                [fplus_avg,fminus_avg,fint_avg,fint_up_avg,fint_down_avg] = CalculateWindHysteresis(u_avg_subwindow,u_ft,u_it,init_state);

                %add to list of values from window-averaged winds - subwindows
                fplus_avg_subwindow{i}(j,k) = fplus_avg; %fraction of time with u_avg above fluid threshold
                fminus_avg_subwindow{i}(j,k) = fminus_avg; %fraction of time with u_avg below impact threshold
                fint_avg_subwindow{i}(j,k) = fint_avg; %fraction of time with u_avg in intermediate zone
                fint_down_avg_subwindow{i}(j,k) = fint_down_avg; %fraction of time with u_avg in intermediate zone with downward crossing
                fint_up_avg_subwindow{i}(j,k) = fint_up_avg; %fraction of time with u_avg in intermediate zone with upward crossing

                %extract data for subwindow
                ind_lowpass_subwindow = find(t_lowpass>=StartTime_subwindow(k)&t_lowpass<EndTime_subwindow(k)); %indices for lowpass points in subwindow
                u_lowpass_subwindow = u_lowpass(ind_lowpass_subwindow); %winds for subwindow
                T_lowpass_subwindow = length(ind_lowpass_subwindow); %number of timesteps in subwindow

                %determine low-pass wind corresponding to fQ - subwindows
                u_lowpass_sort = sort(u_lowpass_subwindow); %sort u's
                ind_uth_lowpass_fQ = round((1-fQ)*T_lowpass_subwindow); %get index in list of u's corresponding to threshold
                if ind_uth_lowpass_fQ==0||isnan(ind_uth_lowpass_fQ)
                    uth_lowpass_fQ = NaN;
                    ustth_lowpass_fQ = NaN;
                    tauth_lowpass_fQ = NaN;
                else
                    uth_lowpass_fQ = u_lowpass_sort(ind_uth_lowpass_fQ); %threshold wind speed
                    ustth_lowpass_fQ = (kappa*uth_lowpass_fQ)/log(zU_window{i}(j)/z0); %threshold shear velocity
                    tauth_lowpass_fQ = rho_a*ustth_lowpass_fQ^2; %threshold shear stress
                end

                %add to list of values from low-pass winds - subwindows
                uth_fQ_lowpass_subwindow{i}(j,k) = uth_lowpass_fQ; %threshold wind from flux frequency
                tauth_fQ_lowpass_subwindow{i}(j,k) = tauth_lowpass_fQ; %threshold wind from flux frequency
                ustth_fQ_lowpass_subwindow{i}(j,k) = ustth_lowpass_fQ; %threshold wind from flux frequency

                %determine fractions of u in different regions - low-pass u
                [fplus_lowpass,fminus_lowpass,fint_lowpass,fint_up_lowpass,fint_down_lowpass] = CalculateWindHysteresis(u_lowpass_subwindow,u_ft,u_it);

                %add to list of values from low-pass winds
                fplus_lowpass_subwindow{i}(j,k) = fplus_lowpass; %fraction of time with u_lowpass above fluid threshold
                fminus_lowpass_subwindow{i}(j,k) = fminus_lowpass; %fraction of time with u_lowpass below impact threshold
                fint_lowpass_subwindow{i}(j,k) = fint_lowpass; %fraction of time with u_lowpass in intermediate zone
                fint_down_lowpass_subwindow{i}(j,k) = fint_down_lowpass; %fraction of time with u_lowpass in intermediate zone with downward crossing
                fint_up_lowpass_subwindow{i}(j,k) = fint_up_lowpass; %fraction of time with u_lowpass in intermediate zone with upward crossing
            end

    %         %% Create sample timeseries plot
    %         figure(1);
    %                 
    %         subplot(2,3,1); cla;
    %         ind_t_flux = find(seconds(t_flux_window{i}{j}-min(t_flux_window{i}{j}))<=t_plot_s);
    %         t_flux_plot = linspace(0,t_plot_s,length(ind_t_flux));
    %         N_flux_plot = sum(n_window{i}{j}(ind_t_flux,:),2)/dt_flux_window(i);
    %         plot(t_flux_plot,N_flux_plot);
    %         xlim([0 t_plot_s]);
    %         ylim_n = ylim;
    %         set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    %         xlabel('t (s)');
    %         ylabel('N (s^{-1})');
    %         title('raw');
    %                 
    %         subplot(2,3,2); cla; hold on;
    %         plot(seconds(t_avg-min(t_avg)),sum(n_avg')');
    %         for k = 1:N_subwindow_plot-1
    %             plot((seconds(delta_t_subwindow)*k)*[1,1],ylim_n,'k--');
    %         end
    %         xlim([0 t_plot_s]);
    %         ylim(ylim_n);
    %         set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    %         xlabel('t (s)');
    %         ylabel('N (s^{-1})');
    %         title('{\delta}t window avg');
    %         
    %         subplot(2,3,3); cla;
    %         bar(t_subwindow_plot,fQ_avg_subwindow{i}(j,1:N_subwindow_plot));
    %         xlim([0 t_plot_s]);
    %         ylim([0 1]);
    %         xlabel('t (s)');
    %         ylabel('f_Q');
    %         title('{\Delta}t stats');
    %         
    %         subplot(2,3,4); cla;
    %         ind_t_wind = find(seconds(t_wind_window{i}{j}-min(t_wind_window{i}{j}))<=t_plot_s);
    %         t_wind_plot = linspace(0,t_plot_s,length(ind_t_wind));
    %         u_plot = u_window{i}{j}(ind_t_wind);
    %         plot(t_wind_plot,u_plot);
    %         xlim([0 t_plot_s]);
    %         ylim_u = ylim;
    %         set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    %         xlabel('t (s)');
    %         ylabel('u (m/s)'); 
    %                 
    %         subplot(2,3,5); cla; hold on;
    %         plot(seconds(t_avg-min(t_avg)),u_avg);
    %         for k = 1:N_subwindow_plot-1
    %             plot((seconds(delta_t_subwindow)*k)*[1,1],ylim_u,'k--');
    %         end
    %         xlim([0 t_plot_s]);
    %         ylim(ylim_u);
    %         set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    %         xlabel('t (s)');
    %         ylabel('u (m/s)');
    %         
    %         subplot(2,3,6); cla;
    %         bar(t_subwindow_plot,uth_fQ_avg_subwindow{i}(j,1:N_subwindow_plot));
    %         xlim([0 t_plot_s]);
    %         ylim(ylim_u);
    %         xlabel('t (s)');
    %         ylabel('u_{th} (m/s)');
    %         
    %         %print plot
    %         set(gca, 'LooseInset', get(gca,'TightInset'));
    %         set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 8]);
    %         print([folder_Plots,'TFEM_demo_',Sites{i},'_',int2str(j),'.png'],'-dpng');

        end
    end

    % SAVE DATA (first clear extraneous values)
    clear('n*','t_*','u_avg*','u_lowpass*');
    save(SaveData_Path,'SiteNames','N_Sites','*avg_window','*lowpass_window','*avg_subwindow','*lowpass_subwindow',...
        'timeofday_subwindow','delta_t_subwindow','zU_window','z0','kappa','rho_a'); %save reduced file in GitHub folder
    
    %% restore deleted values
    t_avg_window = save_t_avg_window_save;
    t_wind_lowpass_window = save_t_wind_lowpass_window_save;
    t_flux_lowpass_window = save_t_flux_lowpass_window_save;
    u_avg_window = save_u_avg_window_save;
    u_lowpass_window = save_u_lowpass_window_save;
    n_avg_window = save_n_avg_window_save;
    n_lowpass_window = save_n_lowpass_window_save;
end