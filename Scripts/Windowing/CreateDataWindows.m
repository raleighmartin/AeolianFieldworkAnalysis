%% SCRIPT TO GENERATE WINDOWS OF FLUX ACTIVITY, SALTATION FLUX, AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
zW_min = 0.018; %minimum Wenglor height (m) = 1.5*height of instrument (to allow one full instrument height between bottom of instrument and bed)
u_sigma_max = 5; %maximum standard deviation in total wind for error detection
AnemometerProfile = {{'U1';'U2';'U3'};{'U1';'U2';'U3'};{'S1';'S2';'S3'}};

%% information about where to load/save data, plots, and functions
folder_TimeData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for windowing data
LoadTimeData_Path = strcat(folder_TimeData,'TimeWindows_30min'); %path for loading time windows
SaveDataWindows_Path = strcat(folder_ProcessedData,'DataWindows_30min'); %path for saving output data
SaveWindProfiles_Path = strcat(folder_ProcessedData,'WindProfiles_30min'); %path for saving output data

%% load time windows
load(LoadTimeData_Path);
SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};

%% load processed data for each site, add to cell arrays of all sites 
WindData = cell(N_Sites,1);
WindProfileData = cell(N_Sites,1);
FluxData = cell(N_Sites,1);
WeatherData = cell(N_Sites,1);
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path); %load processed data
    WindData{i} = ProcessedData.(AnemometerType{i}).(BaseAnemometer{i}); %data only for base anemometer
    WindProfileData{i} = cell(length(AnemometerProfile{i}),1);
    for k = 1:length(AnemometerProfile{i})
        WindProfileData{i}{k} = ProcessedData.(AnemometerType{i}).(AnemometerProfile{i}{k}); %data for anemometer profile
    end
    FluxData{i} = ProcessedData.FluxWenglor; %Wenglor flux data
    WeatherData{i} = ProcessedData.Weather.WS; %all weather station data
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
end

%% initialize variable lists

%initialize lists of flux values (excluding error points)
t_flux_window = cell(N_Sites,1); %times for flux
q_window = cell(N_Sites,1); %partial flux
n_window = cell(N_Sites,1); %counts rate
W_ID_window = cell(N_Sites,1); %get IDs of Wenglors in profiles
zW_window = cell(N_Sites,1); %Wenglor flux height list
sigma_zW_window = cell(N_Sites,1); %Wenglor flux height uncertainty list
N_zW_window = cell(N_Sites,1); %number of Wenglors
zW_min_window = cell(N_Sites,1); %Wenglor min height list
zW_max_window = cell(N_Sites,1); %Wenglor min height list
Cqn_window = cell(N_Sites,1); %calibration factor list
sigma_Cqn_window = cell(N_Sites,1); %calibration factor uncertainty list

%flux values (including interpolated points)
t_flux_int_window = cell(N_Sites,1); %times for interpolated flux
q_int_window = cell(N_Sites,1); %interpolated partial flux
n_int_window = cell(N_Sites,1); %interpolated counts rate
Cqn_int_window = cell(N_Sites,1); %interpolated calibration factor list
sigma_Cqn_int_window = cell(N_Sites,1); %interpolated calibration factor uncertainty list
ind_flux_err_window = cell(N_Sites,1); %times for flux error points
N_flux_err_window = cell(N_Sites,1); %number of flux error points

%initialize lists of wind values for lowest anemometer (excluding error points)
t_wind_window = cell(N_Sites,1); %times for wind
zU_window = cell(N_Sites,1); %height of anemometer
theta_window = cell(N_Sites,1); %wind angle for time window
u_window = cell(N_Sites,1); %rotated streamwise wind for time window
w_window = cell(N_Sites,1); %rotated vertical wind for time window
Temp_window = cell(N_Sites,1); %temperature for time window
RH_window = cell(N_Sites,1); %relative humidity for window

%wind values (including interpolated points)
t_wind_int_window = cell(N_Sites,1); %times for interpolated wind
u_int_window = cell(N_Sites,1); %rotated interpolated streamwise wind for time window
w_int_window = cell(N_Sites,1); %rotated interpolated vertical wind for time window
ind_wind_err_window = cell(N_Sites,1); %times for wind error points
N_wind_err_window = cell(N_Sites,1); %number of wind error points

%initialize mean values for anemometer profile
zU_profile_window = cell(N_Sites,1); %heights of anemometer profile
theta_profile_window = cell(N_Sites,1); %wind angles of anemometer profile
u_profile_window = cell(N_Sites,1); %rotated streamwise wind speed of anemometer profile
w_profile_window = cell(N_Sites,1); %rotated vertical wind speed of anemometer profile
Temp_profile_window = cell(N_Sites,1); %temperature for anemometer profile

% %full profile values
% t_wind_profile = cell(N_Sites,1); %times for wind
% zU_profile = cell(N_Sites,1); %height of anemometer
% theta_profile = cell(N_Sites,1); %wind angle for time window
% u_profile = cell(N_Sites,1); %rotated streamwise wind for time window
% w_profile = cell(N_Sites,1); %rotated vertical wind for time window
% Temp_profile = cell(N_Sites,1); %temperature for time window
% t_wind_int_profile = cell(N_Sites,1); %times for interpolated wind
% u_int_profile = cell(N_Sites,1); %rotated interpolated streamwise wind for time window
% w_int_profile = cell(N_Sites,1); %rotated interpolated vertical wind for time window
% ind_wind_err_profile = cell(N_Sites,1); %times for wind error points
% N_wind_err_profile = cell(N_Sites,1); %number of wind error points

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% get start times and end times from file
    WindowStartTimes = StartTime_window{i};
    WindowEndTimes = EndTime_window{i};
    N_Windows = length(WindowStartTimes);
    
    %% get number of anemometers in profile
    N_Profile = length(AnemometerProfile{i});
    
    %% initialize lists of values
    %flux values (excluding error points)
    t_flux_window{i} = cell(N_Windows,1); %times for flux
    q_window{i} = cell(N_Windows,1); %partial flux
    n_window{i} = cell(N_Windows,1); %counts rate
    W_ID_window{i} = cell(N_Windows,1); %get names of Wenglors in profiles
    zW_window{i} = cell(N_Windows,1); %Wenglor height
    sigma_zW_window{i} = cell(N_Windows,1); %Wenglor height uncertainty
    N_zW_window{i} = zeros(N_Windows,1); %number of Wenglors
    zW_min_window{i} = zeros(N_Windows,1); %minimum Wenglor height
    zW_max_window{i} = zeros(N_Windows,1); %maximum Wenglor height
    Cqn_window{i} = cell(N_Windows,1); %calibration factor
    sigma_Cqn_window{i} = cell(N_Windows,1); %uncertainty in calibration factor
    
    %flux values (including interpolated points)
    t_flux_int_window{i} = cell(N_Windows,1); %times for interpolated flux
    q_int_window{i} = cell(N_Windows,1); %interpolated partial flux
    n_int_window{i} = cell(N_Windows,1); %interpolated counts rate
    Cqn_int_window{i} = cell(N_Windows,1); %interpolated calibration factor list
    sigma_Cqn_int_window{i} = cell(N_Windows,1); %interpolated calibration factor uncertainty list
    ind_flux_err_window{i} = cell(N_Windows,1); %times for flux error points
    N_flux_err_window{i} = zeros(N_Windows,1); %number of flux error points
    
    %wind values (excluding error points)
    t_wind_window{i} = cell(N_Windows,1); %times for wind
    zU_window{i} = zeros(N_Windows,1)*NaN; %height of anemometer
    theta_window{i} = zeros(N_Windows,1)*NaN; %wind angles for anemometer
    u_window{i} = cell(N_Windows,1); %rotated streamwise wind for time window
    w_window{i} = cell(N_Windows,1); %rotated vertical wind for time window
    Temp_window{i} = cell(N_Windows,1); %temperature for time window
    RH_window{i} = zeros(N_Windows,1)*NaN; %relative humidity for window

    %wind values (including interpolated points)
    t_wind_int_window{i} = cell(N_Windows,1); %times for interpolated wind
    u_int_window{i} = cell(N_Windows,1); %rotated interpolated streamwise wind for time window
    w_int_window{i} = cell(N_Windows,1); %rotated interpolated vertical wind for time window
    ind_wind_err_window{i} = cell(N_Windows,1); %times for wind error points
    N_wind_err_window{i} = zeros(N_Windows,1); %number of wind error points
    
    %initialize mean values for anemometer profile
    zU_profile_window{i} = zeros(N_Windows,N_Profile)*NaN; %height of anemometer
    theta_profile_window{i} = zeros(N_Windows,N_Profile)*NaN; %wind angles for anemometer
    u_profile_window{i} = zeros(N_Windows,N_Profile)*NaN; %rotated streamwise wind speed of anemometer profile
    w_profile_window{i} = zeros(N_Windows,N_Profile)*NaN; %rotated vertical wind speed of anemometer profile
    Temp_profile_window{i} = zeros(N_Windows,N_Profile)*NaN; %temperature for anemometer profile

%     %full profile values
%     t_wind_profile{i} = cell(N_Windows,N_Profile); %times for wind
%     zU_profile{i} = zeros(N_Windows,N_Profile)*NaN; %height of anemometer
%     theta_profile{i} = zeros(N_Windows,N_Profile)*NaN; %wind angles for anemometer
%     u_profile{i} = cell(N_Windows,N_Profile); %rotated streamwise wind for time window
%     w_profile{i} = cell(N_Windows,N_Profile); %rotated vertical wind for time window
%     Temp_profile{i} = cell(N_Windows,N_Profile); %temperature for time window
%     t_wind_int_profile{i} = cell(N_Windows,N_Profile); %times for interpolated wind
%     u_int_profile{i} = cell(N_Windows,N_Profile); %rotated interpolated streamwise wind for time window
%     w_int_profile{i} = cell(N_Windows,N_Profile); %rotated interpolated vertical wind for time window
%     ind_wind_err_profile{i} = cell(N_Windows,N_Profile); %times for wind error points
%     N_wind_err_profile{i} = zeros(N_Windows,N_Profile); %number of wind error points
        
    %% go through time windows
    for j = 1:N_Windows

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]

        %get specific start and end time
        StartTime = WindowStartTimes(j);
        EndTime = WindowEndTimes(j);

        %get duration of interval in seconds
        T_interval = seconds(EndTime-StartTime);
        
        %% FLUX DATA FOR INTERVAL

        %extract time interval - get interval number and indices within interval for analysis
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(FluxData{i},StartTime,EndTime,'t','t','t');

        %use only longest interval
        ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
        IntervalN = IntervalN(ind_longest); %number of interval
        IntervalInd = IntervalInd{ind_longest}; %indices for interval
        t_flux = FluxData{i}(IntervalN).t.t(IntervalInd); %times for interval
        
        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(FluxData{i}(IntervalN).t.t,FluxData{i}(IntervalN).t.err); %indices of error times
        ind_flux_err = intersect(ErrInd,IntervalInd)-IntervalInd(1); %get adjusted indices of error times
        IntervalInd_noerr = setdiff(IntervalInd,ErrInd); %indices for interval excluding error times
        t_flux_noerr = FluxData{i}(IntervalN).t.t(IntervalInd_noerr); %times for interval excluding error times
        
        %get Wenglor heights
        zW_profile = mean(FluxData{i}(IntervalN).qz.z(IntervalInd_noerr,:)); %compute mean Wenglor heights
        ind_zW_positive = find(zW_profile>zW_min); %find indices of points definitely above zero (based on zW_min)
        zW_profile = mean(FluxData{i}(IntervalN).qz.z(IntervalInd_noerr,ind_zW_positive)); %compute mean Wenglor heights (but now only positive values)
        sigma_zW_profile = mean(FluxData{i}(IntervalN).qz.sigma_z(IntervalInd_noerr,ind_zW_positive)); %compute uncertainty in Wenglor heights (but now only positive values)
        N_zW = length(zW_profile); %get number of Wenglors
        
        %get raw qz, Cqn, sigma_Cqn and n values
        qz = FluxData{i}(IntervalN).qz.qz(IntervalInd_noerr,ind_zW_positive); %partial fluxes
        n = FluxData{i}(IntervalN).qz.n(IntervalInd_noerr,ind_zW_positive); %particle counts per increment (not per second)
        Cqn = FluxData{i}(IntervalN).qz.Cqn(IntervalInd_noerr,ind_zW_positive); %calibration factors
        sigma_Cqn = FluxData{i}(IntervalN).qz.sigma_Cqn(IntervalInd_noerr,ind_zW_positive); %uncertainty in calibration factors
        W_ID = FluxData{i}(IntervalN).qz.WenglorID(ind_zW_positive); %get associated Wenglor IDs only for heights included in profile
        
        %add to list
        t_flux_window{i}{j} = t_flux_noerr; %times
        q_window{i}{j} = qz; %partial flux
        n_window{i}{j} = n; %counts rate per increment
        W_ID_window{i}{j} = W_ID; %get names of Wenglors in profiles
        zW_window{i}{j} = zW_profile; %Wenglor flux height
        sigma_zW_window{i}{j} = sigma_zW_profile; %uncertainty in Wenglor flux height
        N_zW_window{i}(j) = N_zW; %number of Wenglors
        zW_min_window{i}(j) = min(zW_profile); %minimum Wenglor height
        zW_max_window{i}(j) = max(zW_profile); %maximum Wenglor height
        Cqn_window{i}{j} = Cqn; %calibration factor
        sigma_Cqn_window{i}{j} = sigma_Cqn; %uncertainty in calibration factor
        
        %% get flux values including error points
        qz = FluxData{i}(IntervalN).qz.qz(IntervalInd,ind_zW_positive); %partial fluxes
        n = FluxData{i}(IntervalN).qz.n(IntervalInd,ind_zW_positive)/dt_flux_window(i); %particle counts rate (per second) - account for sampling frequency
        Cqn = FluxData{i}(IntervalN).qz.Cqn(IntervalInd,ind_zW_positive); %calibration factors
        sigma_Cqn = FluxData{i}(IntervalN).qz.sigma_Cqn(IntervalInd,ind_zW_positive); %uncertainty in calibration factors
        
        t_flux_int_window{i}{j} = t_flux; %times for interpolated flux
        q_int_window{i}{j} = qz; %interpolated partial flux
        n_int_window{i}{j} = n; %interpolated counts rate
        Cqn_int_window{i}{j} = Cqn; %interpolated calibration factor list
        sigma_Cqn_int_window{i}{j} = sigma_Cqn; %interpolated calibration factor uncertainty list
      
        %add number of errors to list
        ind_flux_err_window{i}{j} = ind_flux_err; %record list of error time indices
        N_flux_err_window{i}(j) = length(ind_flux_err); %record number of error times
        
        
        %% WIND DATA FOR INTERVAL - BASE ANEMOMETER
        %extract time interval
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData{i},StartTime,EndTime,'u','int','int');
        
        %use only longest interval
        ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
        IntervalN = IntervalN(ind_longest); %number of interval
        IntervalInd = IntervalInd{ind_longest}; %indices for interval
        t_wind = WindData{i}(IntervalN).t.int(IntervalInd); %times for longest interval
        
        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(WindData{i}(IntervalN).t.int,WindData{i}(IntervalN).t.err);
        ind_wind_err = intersect(ErrInd,IntervalInd)-IntervalInd(1); %get adjusted indices of error times
        IntervalInd_noerr = setdiff(IntervalInd,ErrInd); %indices for interval excluding error times
        t_wind_noerr = WindData{i}(IntervalN).t.int(IntervalInd_noerr); %times for interval excluding error times
        
        %get anemometer height
        zU = WindData{i}(IntervalN).z.z;

        %get velocity values using no error points
        u = WindData{i}(IntervalN).u.int(IntervalInd_noerr);
        v = WindData{i}(IntervalN).v.int(IntervalInd_noerr);
        w = WindData{i}(IntervalN).w.int(IntervalInd_noerr);
        Temp = WindData{i}(IntervalN).T.int(IntervalInd_noerr);

        %remove additional error points based on large deviations in wind velocity
        u_total = sqrt(u.^2+v.^2+w.^2); %total wind
        u_total_max = mean(u_total)+u_sigma_max*std(u_total); %maximum total wind based on multiple of std dev
        ind_good_wind = find(u_total<=u_total_max); %get indices of points with total wind below upper limit
        t_wind_noerr = t_wind_noerr(ind_good_wind); %keep only good t
        u = u(ind_good_wind); %keep only good u
        v = v(ind_good_wind); %keep only good v
        w = w(ind_good_wind); %keep only good w
        Temp = Temp(ind_good_wind); %keep only good temperatures
        
        %compute wind angle
        theta = atan(mean(v./u))*180/pi; %wind angle

        %rotate instrument, call these 'rot' values
        [u_rot, ~, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

        %extract humidity values
        [H_Interval, ~, ~, ~] = ExtractVariableTimeInterval(WeatherData{i},StartTime,EndTime,'H','int','int');
        
        %add to window lists
        t_wind_window{i}{j} = t_wind_noerr; %wind times
        theta_window{i}(j) = theta; %wind angle
        zU_window{i}(j) = zU; %height of anemometer
        u_window{i}{j} = u_rot; %rotated streamwise wind
        w_window{i}{j} = w_rot; %rotated vertical wind
        Temp_window{i}{j} = Temp; %temperature
        RH_window{i}(j) = mean(H_Interval); %relative humidity

        %% get wind values including interpolated points
        u_int = WindData{i}(IntervalN).u.int(IntervalInd);
        v_int = WindData{i}(IntervalN).v.int(IntervalInd);
        w_int = WindData{i}(IntervalN).w.int(IntervalInd);
        [u_rot_int, ~, ~] = reorient_anemometers_vanboxel2004(u_int, v_int, w_int); %rotate instrument

        %recompute ind_wind_err to include points with wind outside 5 sigma
        u_total = sqrt(u_int.^2+v_int.^2+w_int.^2); %total wind
        ind_bad_wind = find(u_total>u_total_max); %get indices of points with total wind above upper limit
        ind_wind_err = union(ind_wind_err,ind_bad_wind); %combine these indices with ErrInd points
                
        %add to window lists
        t_wind_int_window{i}{j} = t_wind; %interpolated wind times
        u_int_window{i}{j} = u_rot_int; %rotated interpolated streamwise velocities
        w_int_window{i}{j} = u_rot_int; %rotated interpolated vertical velocities
        
        %add number of errors to list
        ind_wind_err_window{i}{j} = ind_wind_err; %record list of error time indices
        N_wind_err_window{i}(j) = length(ind_wind_err); %record number of error times
        
        %% WIND DATA FOR INTERVAL - PROFILE
        for k = 1:N_Profile
            
            %extract time interval
            [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindProfileData{i}{k},StartTime,EndTime,'u','int','int');

            %use only longest interval
            ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
            IntervalN = IntervalN(ind_longest); %number of interval
            IntervalInd = IntervalInd{ind_longest}; %indices for interval
            t_wind = WindProfileData{i}{k}(IntervalN).t.int(IntervalInd); %times for longest interval

            %further reduce IntervalInd based on eliminating error times
            [~, ErrInd, ~] = intersect(WindProfileData{i}{k}(IntervalN).t.int,WindProfileData{i}{k}(IntervalN).t.err);
            ind_wind_err = intersect(ErrInd,IntervalInd)-IntervalInd(1); %get adjusted indices of error times
            IntervalInd_noerr = setdiff(IntervalInd,ErrInd); %indices for interval excluding error times
            t_wind_noerr = WindProfileData{i}{k}(IntervalN).t.int(IntervalInd_noerr); %times for interval excluding error times

            %get anemometer height
            zU = WindProfileData{i}{k}(IntervalN).z.z;

            %get velocity values using no error points
            u = WindProfileData{i}{k}(IntervalN).u.int(IntervalInd_noerr);
            v = WindProfileData{i}{k}(IntervalN).v.int(IntervalInd_noerr);
            w = WindProfileData{i}{k}(IntervalN).w.int(IntervalInd_noerr);
            Temp = WindProfileData{i}{k}(IntervalN).T.int(IntervalInd_noerr);

            %remove additional error points based on large deviations in wind velocity
            u_total = sqrt(u.^2+v.^2+w.^2); %total wind
            u_total_max = mean(u_total)+u_sigma_max*std(u_total); %maximum total wind based on multiple of std dev
            ind_good_wind = find(u_total<=u_total_max); %get indices of points with total wind below upper limit
            t_wind_noerr = t_wind_noerr(ind_good_wind); %keep only good t
            u = u(ind_good_wind); %keep only good u
            v = v(ind_good_wind); %keep only good v
            w = w(ind_good_wind); %keep only good w
            Temp = Temp(ind_good_wind); %keep only good temperatures

            %compute wind angle
            theta = atan(mean(v./u))*180/pi; %wind angle

            %rotate instrument, call these 'rot' values
            [u_rot, ~, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

            %initialize mean values for anemometer profile
            zU_profile_window{i}(j,k) = zU; %height of anemometer
            theta_profile_window{i}(j,k) = theta; %wind angle
            u_profile_window{i}(j,k) = mean(u_rot); %mean rotated streamwise wind
            w_profile_window{i}(j,k) = mean(w_rot); %mean rotated vertical wind
            Temp_profile_window{i}(j,k) = mean(Temp); %temperature
            
%             %add to window lists - all values
%             t_wind_profile{i}{j,k} = t_wind_noerr; %wind times
%             theta_profile{i}(j,k) = theta; %wind angle
%             zU_profile{i}(j,k) = zU; %height of anemometer
%             u_profile{i}{j,k} = u_rot; %rotated streamwise wind
%             w_profile{i}{j,k} = w_rot; %rotated vertical wind
%             Temp_profile{i}{j,k} = Temp; %temperature
%             
%             %% get wind values including interpolated points
%             u_int = WindProfileData{i}{k}(IntervalN).u.int(IntervalInd);
%             v_int = WindProfileData{i}{k}(IntervalN).v.int(IntervalInd);
%             w_int = WindProfileData{i}{k}(IntervalN).w.int(IntervalInd);
%             [u_rot_int, ~, ~] = reorient_anemometers_vanboxel2004(u_int, v_int, w_int); %rotate instrument
% 
%             %recompute ind_wind_err to include points with wind outside 5 sigma
%             u_total = sqrt(u_int.^2+v_int.^2+w_int.^2); %total wind
%             ind_bad_wind = find(u_total>u_total_max); %get indices of points with total wind above upper limit
%             ind_wind_err = union(ind_wind_err,ind_bad_wind); %combine these indices with ErrInd points
% 
%             %add to window lists
%             t_wind_int_profile{i}{j,k} = t_wind; %interpolated wind times
%             u_int_profile{i}{j,k} = u_rot_int; %rotated interpolated streamwise velocities
%             w_int_profile{i}{j,k} = w_rot_int; %rotated interpolated streamwise velocities
%             
%             %add number of errors to list
%             ind_wind_err_profile{i}{j,k} = ind_wind_err; %record list of error time indices
%             N_wind_err_profile{i}(j,k) = length(ind_wind_err); %record number of error times
        end
    end
end

% SAVE DATA
save(SaveDataWindows_Path,'Sites','SiteNames','N_Sites','BaseAnemometer','*window'); %save primary data reduced file in ProcessedData folder
%save(SaveWindProfiles_Path,'Sites','SiteNames','N_Sites','AnemometerProfile','*profile'); %save wind profile reduced file in ProcessedData folder