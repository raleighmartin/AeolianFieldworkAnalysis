%% SCRIPT TO GENERATE WINDOWS OF FLUX ACTIVITY, SALTATION FLUX, AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;
close all;

%% parameter values
zW_min = 0.018; %minimum Wenglor height (m) = 1.5*height of instrument (to allow one full instrument height between bottom of instrument and bed)
u_sigma_max = 5; %maximum standard deviation in total wind for error detection

%% information about where to load/save data, plots, and functions
folder_TimeData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

% information for sample plots
folder_TimeseriesPlot = '../../PlotOutput/Timeseries/'; %folder for sample timeseries plots
dt_plot = duration(0,0,5); %time for window-average for plot
h1 = figure(1);
set(h1,'visible','off');

%% Specific information for windowing data - restricted windows
LoadTimeData_Path = strcat(folder_TimeData,'TimeWindows_30min_Restricted'); %path for loading time windows
SaveDataWindows_Path = strcat(folder_ProcessedData,'DataWindows_30min_Restricted'); %path for saving output data

% %% Specific information for windowing data - unrestricted windows
% LoadTimeData_Path = strcat(folder_TimeData,'TimeWindows_30min_Unrestricted'); %path for loading time windows
% SaveDataWindows_Path = strcat(folder_ProcessedData,'DataWindows_30min_Unrestricted'); %path for saving output data

%% load time windows
load(LoadTimeData_Path);

%% load processed data for each site, add to cell arrays of all sites 
FluxData = cell(N_Sites,1);
WeatherData = cell(N_Sites,1);
WindData_base = cell(N_Sites,1);
WindData_profile = cell(N_Sites,1);
N_Anemometers = zeros(N_Sites,1); %number of anemometers
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path); %load processed data
    
    %get flux and weather data
    FluxData{i} = ProcessedData.FluxWenglor; %Wenglor flux data
    WeatherData{i} = ProcessedData.Weather.WS; %all weather station data
    
    %get wind data
    WindData = ProcessedData.(AnemometerType{i}); %data only for base anemometer
    N_Anemometers(i) = length(AnemometerName_profile{i}); %number of anemometers
    WindData_base{i} = WindData.(AnemometerName_base{i}); %data only for base anemometer
    WindData_profile{i} = cell(length(AnemometerName_profile{i}),1);
    for k = 1:length(AnemometerName_profile{i})
        WindData_profile{i}{k} = WindData.(AnemometerName_profile{i}{k}); %data for anemometer profile
    end
        
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
end

%% initialize variable lists

%initialize lists of flux values (excluding error points)
t_flux_window = cell(N_Sites,1); %times for flux
q_window = cell(N_Sites,1); %partial flux
n_window = cell(N_Sites,1); %counts rate (per Wenglor per increment)
ntotal_window = cell(N_Sites,1); %total counts rate (all Wenglors per second)
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
ntotal_int_window = cell(N_Sites,1); %interpolated total counts rate
Cqn_int_window = cell(N_Sites,1); %interpolated calibration factor list
sigma_Cqn_int_window = cell(N_Sites,1); %interpolated calibration factor uncertainty list
ind_flux_err_window = cell(N_Sites,1); %times for flux error points
N_flux_err_window = cell(N_Sites,1); %number of flux error points

%initialize lists of wind values for lowest anemometer (excluding error points)
t_wind_window = cell(N_Sites,1); %times for wind
theta_window = cell(N_Sites,1); %wind angle for time window
u_window = cell(N_Sites,1); %rotated streamwise wind for time window
v_window = cell(N_Sites,1); %rotated transverse wind for time window
w_window = cell(N_Sites,1); %rotated vertical wind for time window
Temp_window = cell(N_Sites,1); %temperature for time window
RH_window = cell(N_Sites,1); %relative humidity for window

%wind values (including interpolated points)
t_wind_int_window = cell(N_Sites,1); %times for interpolated wind
u_int_window = cell(N_Sites,1); %rotated interpolated streamwise wind for time window
v_int_window = cell(N_Sites,1); %rotated interpolated transverse wind for time window
w_int_window = cell(N_Sites,1); %rotated interpolated vertical wind for time window
ind_wind_err_window = cell(N_Sites,1); %times for wind error points
N_wind_err_window = cell(N_Sites,1); %number of wind error points

%full profile values (not including interpolated points)
theta_profile_window = cell(N_Sites,1); %wind angle for time window
u_profile_window = cell(N_Sites,1); %rotated streamwise wind for time window
v_profile_window = cell(N_Sites,1); %rotated transverse wind for time window
w_profile_window = cell(N_Sites,1); %rotated vertical wind for time window
Temp_profile_window = cell(N_Sites,1); %temperature for time window

%% combine start and end times for windows with and without flux data
hasfluxdata_window = cell(N_Sites,1);

for i = 1:N_Sites
    %list of zeros for windows with flux data, list of ones for windows with no flux data
    hasfluxdata_window{i} = [ones(size(StartTime_window{i})); zeros(size(StartTime_window_noflux{i}))];
   
    % get combined start times and end times - replace existing arrays
    StartTime_window{i} = [StartTime_window{i}; StartTime_window_noflux{i}]; %include windows with no flux data
    EndTime_window{i} = [EndTime_window{i}; EndTime_window_noflux{i}]; %include windows with no flux data
    Date_window{i} = [Date_window{i}; Date_window_noflux{i}]; %include windows with no flux data
    zU_base_window{i} = [zU_base_window{i}; zU_base_window_noflux{i}]; %include windows with no flux data
    zU_profile_window{i} = [zU_profile_window{i}; zU_profile_window_noflux{i}]; %include windows with no flux data
end

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% get start times and end times
    WindowStartTimes = StartTime_window{i};
    WindowEndTimes = EndTime_window{i};
    N_Windows = length(WindowStartTimes);
  
    %% initialize lists of values
    %flux values (excluding error points)
    t_flux_window{i} = cell(N_Windows,1); %times for flux
    q_window{i} = cell(N_Windows,1); %partial flux
    n_window{i} = cell(N_Windows,1); %counts rate
    ntotal_window{i} = cell(N_Windows,1); %total counts rate
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
    ntotal_int_window{i} = cell(N_Windows,1); %interpolated counts rate
    Cqn_int_window{i} = cell(N_Windows,1); %interpolated calibration factor list
    sigma_Cqn_int_window{i} = cell(N_Windows,1); %interpolated calibration factor uncertainty list
    ind_flux_err_window{i} = cell(N_Windows,1); %times for flux error points
    N_flux_err_window{i} = zeros(N_Windows,1); %number of flux error points
    
    %wind values (excluding error points)
    t_wind_window{i} = cell(N_Windows,1); %times for wind
    theta_window{i} = zeros(N_Windows,1)*NaN; %wind angles for anemometer
    u_window{i} = cell(N_Windows,1); %rotated streamwise wind for time window
    v_window{i} = cell(N_Windows,1); %rotated streamwise wind for time window
    w_window{i} = cell(N_Windows,1); %rotated vertical wind for time window
    Temp_window{i} = cell(N_Windows,1); %temperature for time window
    RH_window{i} = zeros(N_Windows,1)*NaN; %relative humidity for window

    %wind values (including interpolated points)
    t_wind_int_window{i} = cell(N_Windows,1); %times for interpolated wind
    u_int_window{i} = cell(N_Windows,1); %rotated interpolated streamwise wind for time window
    v_int_window{i} = cell(N_Windows,1); %rotated interpolated transverse wind for time window
    w_int_window{i} = cell(N_Windows,1); %rotated interpolated vertical wind for time window
    ind_wind_err_window{i} = cell(N_Windows,1); %times for wind error points
    N_wind_err_window{i} = zeros(N_Windows,1); %number of wind error points
    
    %mean profile values (not including interpolated values)
    theta_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %wind angles for anemometer
    u_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %rotated streamwise wind for time window
    v_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %rotated transverse wind for time window
    w_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %rotated vertical wind for time window
    Temp_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %temperature for time window

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

        %go through computations if window has flux data
        if hasfluxdata_window{i}(j)==1
            
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
            zW = mean(FluxData{i}(IntervalN).qz.z(IntervalInd_noerr,:)); %compute mean Wenglor heights
            ind_zW_positive = find(zW>zW_min); %find indices of points definitely above zero (based on zW_min)
            zW = mean(FluxData{i}(IntervalN).qz.z(IntervalInd_noerr,ind_zW_positive)); %compute mean Wenglor heights (but now only positive values)
            sigma_zW = mean(FluxData{i}(IntervalN).qz.sigma_z(IntervalInd_noerr,ind_zW_positive)); %compute uncertainty in Wenglor heights (but now only positive values)
            N_zW = length(zW); %get number of Wenglors

            %get raw qz, Cqn, sigma_Cqn and n values
            qz = FluxData{i}(IntervalN).qz.qz(IntervalInd_noerr,ind_zW_positive); %partial fluxes
            n = FluxData{i}(IntervalN).qz.n(IntervalInd_noerr,ind_zW_positive); %particle counts per increment (not per second)
            ntotal = sum(n')'/dt_flux_window(i); %total counts rate per second
            Cqn = FluxData{i}(IntervalN).qz.Cqn(IntervalInd_noerr,ind_zW_positive); %calibration factors
            sigma_Cqn = FluxData{i}(IntervalN).qz.sigma_Cqn(IntervalInd_noerr,ind_zW_positive); %uncertainty in calibration factors
            W_ID = FluxData{i}(IntervalN).qz.WenglorID(ind_zW_positive); %get associated Wenglor IDs only for heights included in profile
            
            %add to list
            t_flux_window{i}{j} = t_flux_noerr; %times
            q_window{i}{j} = qz; %partial flux
            n_window{i}{j} = n; %counts rate per increment
            ntotal_window{i}{j} = ntotal; %total counts rate per increment
            W_ID_window{i}{j} = W_ID; %get names of Wenglors in profiles
            zW_window{i}{j} = zW; %Wenglor flux height
            sigma_zW_window{i}{j} = sigma_zW; %uncertainty in Wenglor flux height
            N_zW_window{i}(j) = N_zW; %number of Wenglors
            zW_min_window{i}(j) = min(zW); %minimum Wenglor height
            zW_max_window{i}(j) = max(zW); %maximum Wenglor height
            Cqn_window{i}{j} = Cqn; %calibration factor
            sigma_Cqn_window{i}{j} = sigma_Cqn; %uncertainty in calibration factor

            %% get flux values including error points
            qz = FluxData{i}(IntervalN).qz.qz(IntervalInd,ind_zW_positive); %partial fluxes
            n = FluxData{i}(IntervalN).qz.n(IntervalInd,ind_zW_positive); %interpolated particle counts rate
            ntotal = sum(n')'/dt_flux_window(i); %total counts rate per second
            Cqn = FluxData{i}(IntervalN).qz.Cqn(IntervalInd,ind_zW_positive); %calibration factors
            sigma_Cqn = FluxData{i}(IntervalN).qz.sigma_Cqn(IntervalInd,ind_zW_positive); %uncertainty in calibration factors

            %add to list
            t_flux_int_window{i}{j} = t_flux; %times for interpolated flux
            q_int_window{i}{j} = qz; %interpolated partial flux
            n_int_window{i}{j} = n; %interpolated counts rate
            ntotal_int_window{i}{j} = ntotal; %interpolated total counts rate
            Cqn_int_window{i}{j} = Cqn; %interpolated calibration factor list
            sigma_Cqn_int_window{i}{j} = sigma_Cqn; %interpolated calibration factor uncertainty list
            ind_flux_err_window{i}{j} = ind_flux_err; %record list of error time indices
            N_flux_err_window{i}(j) = length(ind_flux_err); %record number of error times
        
        %otherwise set all values as NaN if no flux data exist
        else
            %non-error values
            t_flux_window{i}{j} = NaN; %times
            q_window{i}{j} = NaN; %partial flux
            n_window{i}{j} = NaN; %counts rate per increment
            ntotal_window{i}{j} = NaN; %total counts rate
            W_ID_window{i}{j} = NaN; %get names of Wenglors in profiles
            zW_window{i}{j} = NaN; %Wenglor flux height
            sigma_zW_window{i}{j} = NaN; %uncertainty in Wenglor flux height
            N_zW_window{i}(j) = NaN; %number of Wenglors
            zW_min_window{i}(j) = NaN; %minimum Wenglor height
            zW_max_window{i}(j) = NaN; %maximum Wenglor height
            Cqn_window{i}{j} = NaN; %calibration factor
            sigma_Cqn_window{i}{j} = NaN; %uncertainty in calibration factor
            %interpolated values
            t_flux_int_window{i}{j} = NaN; %times for interpolated flux
            q_int_window{i}{j} = NaN; %interpolated partial flux
            n_int_window{i}{j} = NaN; %interpolated counts rate
            ntotal_int_window{i}{j} = NaN; %total counts rate
            Cqn_int_window{i}{j} = NaN; %interpolated calibration factor list
            sigma_Cqn_int_window{i}{j} = NaN; %interpolated calibration factor uncertainty list
            ind_flux_err_window{i}{j} = NaN; %record list of error time indices
            N_flux_err_window{i}(j) = NaN; %record number of error times
        end
            
        %% WIND DATA FOR INTERVAL - BASE ANEMOMETER
        %extract time interval
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData_base{i},StartTime,EndTime,'u','int','int');
        
        %use only longest interval
        ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
        IntervalN = IntervalN(ind_longest); %number of interval
        IntervalInd = IntervalInd{ind_longest}; %indices for interval
        t_wind = WindData_base{i}(IntervalN).t.int(IntervalInd); %times for longest interval
        
        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(WindData_base{i}(IntervalN).t.int,WindData_base{i}(IntervalN).t.err);
        ind_wind_err = intersect(ErrInd,IntervalInd)-IntervalInd(1); %get adjusted indices of error times
        IntervalInd_noerr = setdiff(IntervalInd,ErrInd); %indices for interval excluding error times
        t_wind_noerr = WindData_base{i}(IntervalN).t.int(IntervalInd_noerr); %times for interval excluding error times

        %get velocity values using no error points
        u = WindData_base{i}(IntervalN).u.int(IntervalInd_noerr);
        v = WindData_base{i}(IntervalN).v.int(IntervalInd_noerr);
        w = WindData_base{i}(IntervalN).w.int(IntervalInd_noerr);
        Temp = WindData_base{i}(IntervalN).T.int(IntervalInd_noerr);

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
        [u_rot, v_rot, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

        %extract humidity values
        [H_Interval, ~, ~, ~] = ExtractVariableTimeInterval(WeatherData{i},StartTime,EndTime,'H','int','int');
        
        %add to window lists
        t_wind_window{i}{j} = t_wind_noerr; %wind times
        theta_window{i}(j) = theta; %wind angle
        u_window{i}{j} = u_rot; %rotated streamwise wind
        v_window{i}{j} = v_rot; %rotated transverse wind
        w_window{i}{j} = w_rot; %rotated vertical wind
        Temp_window{i}{j} = Temp; %temperature
        RH_window{i}(j) = mean(H_Interval); %relative humidity

        %% get wind values including interpolated points
        u_int = WindData_base{i}(IntervalN).u.int(IntervalInd);
        v_int = WindData_base{i}(IntervalN).v.int(IntervalInd);
        w_int = WindData_base{i}(IntervalN).w.int(IntervalInd);
        [u_rot_int, v_rot_int, w_rot_int] = reorient_anemometers_vanboxel2004(u_int, v_int, w_int); %rotate instrument

        %recompute ind_wind_err to include points with wind outside 5 sigma
        u_total = sqrt(u_int.^2+v_int.^2+w_int.^2); %total wind
        ind_bad_wind = find(u_total>u_total_max); %get indices of points with total wind above upper limit
        ind_wind_err = union(ind_wind_err,ind_bad_wind); %combine these indices with ErrInd points
                
        %add to window lists
        t_wind_int_window{i}{j} = t_wind; %interpolated wind times
        u_int_window{i}{j} = u_rot_int; %rotated interpolated streamwise velocities
        v_int_window{i}{j} = v_rot_int; %rotated interpolated transverse velocities
        w_int_window{i}{j} = w_rot_int; %rotated interpolated vertical velocities
        
        %add number of errors to list
        ind_wind_err_window{i}{j} = ind_wind_err; %record list of error time indices
        N_wind_err_window{i}(j) = length(ind_wind_err); %record number of error times

        %% create sample timeseries plot
        
        %get window averages
        t_plot = t_wind(1):dt_plot:t_wind(end-1); %times
        u_plot = window_average(u_rot_int, t_wind, dt_plot); %compute window average
        w_plot = window_average(w_rot_int, t_wind, dt_plot); %compute window average
        qz_plot = zeros(length(t_plot),N_zW);
        for k=1:N_zW
            qz_plot(:,k) = window_average(qz(:,k), t_flux, dt_plot); %compute window average
        end
        
        %create plot
        subplot(6,1,1:2);
        plot(t_plot,u_plot);
        ylabel('streamwise wind speed, $$u$$ (ms$$^{-1}$$)','Interpreter','Latex');
        title([SiteNames{i},', ',datestr(StartTime, 'yyyy-mm-dd HH:MM'),' - ',datestr(EndTime, 'HH:MM')]);
        set(gca,'FontSize',13,'XMinorTick','On','YMinorTick','On','Box','On');
        ylims = [floor(min(u_plot)) ceil(max(u_plot))];
        ylim(ylims);
        text(min(t_plot)+duration(0,0,30),ylims(2)*0.97,'(a)','FontSize',13);

        subplot(6,1,3);
        plot(t_plot,w_plot);
        ylabel('vert. wind speed, $$w$$ (ms$$^{-1}$$)','Interpreter','Latex');
        set(gca,'FontSize',13,'XMinorTick','On','YMinorTick','On','Box','On');
        ylims = [floor(min(w_plot)*10) ceil(max(w_plot)*10)]/10;
        ylim(ylims);
        text(min(t_plot)+duration(0,0,30),ylims(2)*0.65,'(b)','FontSize',13);

        subplot(6,1,4:6); cla; hold on;
        legend_items = cell(N_zW,1);
        for k=1:N_zW
            plot(t_plot,qz_plot(:,k));
            legend_items{k} = ['$$z_',int2str(k),'$$ = ',num2str(round(zW(k)*100,1)),' cm'];
        end
        legend(legend_items,'Location','NorthWest','Interpreter','Latex');
        ylabel('partial saltation flux, $$q_i$$ (gm$$^{-2}$$s$$^{-1}$$)','Interpreter','Latex');
        set(gca,'FontSize',13,'XMinorTick','On','YMinorTick','On','Box','On');
        ylims = [0 ceil(max(max(qz_plot))/10)]*10;
        if range(ylims)==0; ylims = [0 10]; end; %adjust ylims if no saltation
        ylim(ylims);
        text(min(t_plot)+duration(0,7,0),ylims(2)*0.95,'(c)','FontSize',13);
        set(gcf, 'PaperPosition',[0 0 10.5 13]);
        print([folder_TimeseriesPlot,'SampleTimeseries_',Sites{i},'_',int2str(j),'.png'],'-dpng');

        %% WIND DATA FOR INTERVAL - PROFILE
        for k = 1:N_Anemometers(i)
            
            %only perform extraction if there are data in window
            if ~isnan(zU_profile_window{i}(j,k))
                
                %extract time interval
                [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData_profile{i}{k},StartTime,EndTime,'u','int','int');

                %use only longest interval
                ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
                IntervalN = IntervalN(ind_longest); %number of interval
                IntervalInd = IntervalInd{ind_longest}; %indices for interval
                %t_wind = WindData_profile{i}{k}(IntervalN).t.int(IntervalInd); %times for longest interval

                %further reduce IntervalInd based on eliminating error times
                [~, ErrInd, ~] = intersect(WindData_profile{i}{k}(IntervalN).t.int,WindData_profile{i}{k}(IntervalN).t.err);
                ind_wind_err = intersect(ErrInd,IntervalInd)-IntervalInd(1); %get adjusted indices of error times
                IntervalInd_noerr = setdiff(IntervalInd,ErrInd); %indices for interval excluding error times
                %t_wind_noerr = WindData_profile{i}{k}(IntervalN).t.int(IntervalInd_noerr); %times for interval excluding error times

                %get velocity values using no error points
                u = WindData_profile{i}{k}(IntervalN).u.int(IntervalInd_noerr);
                v = WindData_profile{i}{k}(IntervalN).v.int(IntervalInd_noerr);
                w = WindData_profile{i}{k}(IntervalN).w.int(IntervalInd_noerr);
                Temp = WindData_profile{i}{k}(IntervalN).T.int(IntervalInd_noerr);

                %remove additional error points based on large deviations in wind velocity
                u_total = sqrt(u.^2+v.^2+w.^2); %total wind
                u_total_max = mean(u_total)+u_sigma_max*std(u_total); %maximum total wind based on multiple of std dev
                ind_good_wind = find(u_total<=u_total_max); %get indices of points with total wind below upper limit
                %t_wind_noerr = t_wind_noerr(ind_good_wind); %keep only good t
                u = u(ind_good_wind); %keep only good u
                v = v(ind_good_wind); %keep only good v
                w = w(ind_good_wind); %keep only good w
                Temp = Temp(ind_good_wind); %keep only good temperatures

                %compute wind angle
                theta = atan(mean(v./u))*180/pi; %wind angle

                %rotate instrument, call these 'rot' values
                [u_rot, v_rot, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

                %add to window lists - mean values
                theta_profile_window{i}(j,k) = theta; %wind angle
                u_profile_window{i}(j,k) = mean(u_rot); %rotated streamwise wind
                v_profile_window{i}(j,k) = mean(v_rot); %rotated transverse wind
                w_profile_window{i}(j,k) = mean(w_rot); %rotated vertical wind
                Temp_profile_window{i}(j,k) = mean(Temp); %temperature
            end
        end
    end
end

% SAVE DATA
save(SaveDataWindows_Path,'Sites','SiteNames','N_Sites','AnemometerName_base','AnemometerName_profile','*window'); %save primary data reduced file in ProcessedData folder