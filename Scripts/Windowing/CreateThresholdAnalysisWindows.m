%% SCRIPT TO OBTAIN EFFECTIVE THRESHOLDS BY TFEM METHOD
% Deltat = measurement interval duration (longer)
% deltat = sampling interval duration (shorter)

%% initialize
clearvars;

%% measurement interval times
Deltat_all = ...
    [duration(0,0,30),...
    duration(0,1,0),...
    duration(0,2,0),...
    duration(0,5,0),...
    duration(0,10,0)];
N_Deltat = length(Deltat_all); %number of measurement intervals
subwindow_per_window = duration(0,30,0)./Deltat_all; %get number of measurement interval subwindows per window

%% sampling interval times
deltat_all = ... %time intervals for window averaging
    [duration(0,0,0.5),...
    duration(0,0,0.72),...
    duration(0,0,1),...
    duration(0,0,1.4),...
    duration(0,0,2)];
N_deltat = length(deltat_all); %number of sampling intervals

%% calculate samples per subwindow
N_samples_per_subwindow = zeros(N_Deltat,N_deltat);
for m = 1:N_Deltat
    for s = 1:N_deltat
        N_samples_per_subwindow(m,s) = floor(Deltat_all(m)/deltat_all(s));
    end
end

%% parameter values
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
z0 = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m)
g = 9.8; %gravity m/s^2
kappa = 0.4; %von Karman parameter
fD_min = 0.005; %minimum detection rate for zero flux

%% load data 
folder_LoadData_1 = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
LoadData_1_Path = strcat(folder_LoadData_1,'TimeWindows_30min.mat'); %path for loading time window data
load(LoadData_1_Path);
folder_LoadData_2 = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
LoadData_2_Path = strcat(folder_LoadData_2,'DataWindows_30min.mat'); %path for loading full data windows
load(LoadData_2_Path);
folder_LoadData_3 = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
LoadData_3_Path = strcat(folder_LoadData_3,'FluxLawWindows_30min.mat'); %path for loading flux law windows
load(LoadData_3_Path);

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% specify where to save and plot data
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
SaveData_Path = strcat(folder_SaveData,'ThresholdAnalysisWindows'); %path for saving output data
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% initialize variable arrays
fD_all = cell(N_Sites,1); %Wenglor detection frequency array
fQ_all = cell(N_Sites,1); %Wenglor transport frequency array
lambda_all = cell(N_Sites,1); %estimated particle arrival rate array
zU_all = cell(N_Sites,1); %anemometer height
zs_all = cell(N_Sites,1); %roughness height array
uth_all = cell(N_Sites,1); %effective threshold wind
ustth_all = cell(N_Sites,1); %effective threshold shear velocity
tauth_all = cell(N_Sites,1); %effective threshold stress
timeofday_all = cell(N_Sites,1); %time of day array
starttime_all = cell(N_Sites,1); %start time array
endtime_all = cell(N_Sites,1); %end time array
t_all = cell(N_Sites,1); %timeseries times array
u_all = cell(N_Sites,1); %wind timeseries array
n_all = cell(N_Sites,1); %counts timeseries array

% %% GO THROUGH SITES
for i = 1:N_Sites
    %% get information about windows and subwindows
    N_windows = length(StartTime_window{i}); %get number of windows
    N_subwindows = N_windows*subwindow_per_window; %get number of subwindows for each measurement interval

    %% initialize variable arrays - by measurement interval and sampling interval
    fD_all{i} = cell(N_Deltat,N_deltat); %Wenglor detection frequency array
    fQ_all{i} = cell(N_Deltat,N_deltat); %Wenglor transport frequency array
    lambda_all{i} = cell(N_Deltat,N_deltat); %estimated particle arrival rate array
    zU_all{i} = cell(N_Deltat,N_deltat); %anemometer height
    zs_all{i} = cell(N_Deltat,N_deltat); %roughness height array
    uth_all{i} = cell(N_Deltat,N_deltat); %effective threshold wind
    ustth_all{i} = cell(N_Deltat,N_deltat); %effective threshold shear velocity
    tauth_all{i} = cell(N_Deltat,N_deltat); %effective threshold stress
    timeofday_all{i} = cell(N_Deltat,N_deltat); %time of day array
    starttime_all{i} = cell(N_Deltat,N_deltat); %start time array
    endtime_all{i} = cell(N_Deltat,N_deltat); %end time array
    t_all{i} = cell(N_Deltat,N_deltat); %timeseries times array
    u_all{i} = cell(N_Deltat,N_deltat); %wind timeseries array
    n_all{i} = cell(N_Deltat,N_deltat); %counts timeseries array
    
    %% initialize subwindow values
    for m = 1:N_Deltat
        for s = 1:N_deltat
            fD_all{i}{m,s} = zeros(N_subwindows(m),1); %Wenglor detection frequency array
            fQ_all{i}{m,s} = zeros(N_subwindows(m),1); %Wenglor transport frequency array
            lambda_all{i}{m,s} = zeros(N_subwindows(m),1); %estimated particle arrival rate array
            zU_all{i}{m,s} = zeros(N_subwindows(m),1); %anemometer height array
            zs_all{i}{m,s} = zeros(N_subwindows(m),1); %roughness height array
            uth_all{i}{m,s} = zeros(N_subwindows(m),1); %effective threshold wind
            ustth_all{i}{m,s} = zeros(N_subwindows(m),1); %effective threshold shear velocity
            tauth_all{i}{m,s} = zeros(N_subwindows(m),1); %effective threshold stress
            timeofday_all{i}{m,s} = zeros(N_subwindows(m),1); %time of day array
            starttime_all{i}{m,s} = datetime(zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1)); %start time array
            endtime_all{i}{m,s} = datetime(zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1),zeros(N_subwindows(m),1)); %end time array            
            t_all{i}{m,s} = cell(N_subwindows(m),1); %timeseries times array
            u_all{i}{m,s} = cell(N_subwindows(m),1); %wind timeseries array
            n_all{i}{m,s} = cell(N_subwindows(m),1); %counts timeseries array
        end
    end
  
    %% GO THROUGH 30-MINUTE WINDOWS
    for j = 1:N_windows
        
        %% get anemometer and roughness height
        zU = zU_window{i}(j); %anemometer height
        zs = zs_window{i}(j); %roughness height
        
        %% get needed flux values for analysis
        n = sum(n_int_window{i}{j},2); %interpolated total counts rate (counts/s)
        t_flux = t_flux_int_window{i}{j}; %times for interpolated flux
        ind_flux_err = ind_flux_err_window{i}{j}; %list of error time indices for flux
        flux_err_binary = zeros(length(n),1); %initialize timeseries with 0's and 1's for error points
        flux_err_binary(ind_flux_err) = 1; %set error indices to 1
        
        %% get needed wind values for analysis
        u = u_int_window{i}{j}; %rotated interpolated streamwise velocities
        t_wind = t_wind_int_window{i}{j}; %times for interpolated wind
        ind_wind_err = ind_wind_err_window{i}{j}; %list of error time indices for wind
        wind_err_binary = zeros(length(u),1); %initialize timeseries with 0's and 1's for error points
        wind_err_binary(ind_wind_err) = 1; %set error indices to 1
        
        %% GO THROUGH SAMPLING INTERVALS
        for s = 1:N_deltat
        
            %% compute window averages
            [n_deltat, t_deltat_n] = window_average(n, t_flux, deltat_all(s)); %counts rate timeseries (counts/s)
            [u_deltat, t_deltat_u] = window_average(u, t_wind, deltat_all(s)); %wind timeseries
            flux_err_binary_deltat = window_average(flux_err_binary, t_flux, deltat_all(s)); %binary flux errors
            wind_err_binary_deltat = window_average(wind_err_binary, t_wind, deltat_all(s)); %binary wind errors

            %% get common t, reduce window average timeseries based on these values
            [t_deltat,ind_n,ind_u] = intersect(t_deltat_n,t_deltat_u);
            n_deltat = n_deltat(ind_n);
            u_deltat = u_deltat(ind_u);
            flux_err_binary_deltat = flux_err_binary_deltat(ind_n);
            wind_err_binary_deltat = wind_err_binary_deltat(ind_u);
            err_binary_deltat = flux_err_binary_deltat+wind_err_binary_deltat; %combine window averages

            %% indices of non-error points
            ind_noerr_deltat = find(err_binary_deltat==0); %non-error indices are window-averaged points == 0
            
            %% go through measurement interval durations
            for m = 1:N_Deltat
                
                %% get subwindow times for measurement interval
                StartTimes_Deltat = StartTime_window{i}(j)+((1:subwindow_per_window(m))-1)*Deltat_all(m);
                EndTimes_Deltat = StartTime_window{i}(j)+(1:subwindow_per_window(m))*Deltat_all(m);
                
                %% go through subwindows
                for k = 1:subwindow_per_window(m)
                    
                    %get all data in subwindow - including error points
                    window_ind_subwindow = find(t_deltat>=StartTimes_Deltat(k)&t_deltat<EndTimes_Deltat(k)); %get indices of all subwindow points within window
                    T_subwindow = length(window_ind_subwindow); %number of timesteps
                    t_subwindow = t_deltat(window_ind_subwindow); %times for subwindow
                    n_subwindow = n_deltat(window_ind_subwindow); %counts for subwindow
                    u_subwindow = u_deltat(window_ind_subwindow); %winds for subwindow

                    %get data in subwindow - no error points
                    ind_subwindow_noerr = intersect(window_ind_subwindow,ind_noerr_deltat); %get indices only of subwindow points with no error (replace existing)
                    T_subwindow_noerr = length(ind_subwindow_noerr); %number of timesteps
                    t_subwindow_noerr = t_deltat(ind_subwindow_noerr); %times for subwindow
                    n_subwindow_noerr = n_deltat(ind_subwindow_noerr); %counts for subwindow
                    u_subwindow_noerr = u_deltat(ind_subwindow_noerr); %winds for subwindow
                    
                    %calculate detection rate and mean counts rate (only non-error points)
                    fD = sum(n_subwindow_noerr>0)/T_subwindow_noerr; %detection rate
                    n_bar = mean(n_subwindow_noerr); %mean counts rate
                    
                    %estimate particle arrival rate per averaging window
                    if fD<=fD_min %set to zero if below detection limit
                        lambda = 0;
                    else %otherwise estimate arrival rate based on fD
                        lambda = n_bar*seconds(deltat_all(s))/fD;
                    end

                    %estimate flux activity from flux detection rate and particle arrival rate
                    if lambda==0 %if no (or negligible) flux, set frequencies to zero
                        fQ = 0;
                    elseif fD==1 %if fD = 1, set fQ to 1
                        fQ = 1;
                    else %otherwise, estimate fQ based on other parameters
                        fQ = fD/(1-exp(-lambda)); %calculate fQ
                        if fQ>1
                            fQ = 1; %if correction gives fQ>1, just set as 1
                        end
                    end
           
                    %determine wind corresponding to fQ (only non-error points)
                    u_sort = sort(u_subwindow_noerr); %sort u's
                    ind_uth_fQ = round((1-fQ)*T_subwindow_noerr); %get index in list of u's corresponding to threshold
                    if (ind_uth_fQ==0||ind_uth_fQ==T_subwindow_noerr)||isnan(ind_uth_fQ) %points don't count if fQ = 0, 1, or is undefined
                        uth = NaN;
                        ustth = NaN;
                        tauth = NaN;
                    else
                        uth = u_sort(ind_uth_fQ); %threshold wind speed
                        ustth = (kappa*uth)/log(zU_window{i}(j)/z0(i)); %threshold shear velocity
                        tauth = rho_a(i)*ustth^2; %threshold shear stress
                    end
                    
                    %calculate time of day
                    timeofday = hour(mean(t_subwindow))+minute(mean(t_subwindow))/60+second(mean(t_subwindow))/3600;
                    
                    %get indices for adding values to array
                    array_ind_subwindow = (j-1)*subwindow_per_window(m)+k;
                    
                    %add to arrays of values
                    fD_all{i}{m,s}(array_ind_subwindow) = fD; %detection frequency
                    fQ_all{i}{m,s}(array_ind_subwindow) = fQ; %flux frequency
                    lambda_all{i}{m,s}(array_ind_subwindow) = lambda; %estimated particle detection rate
                    zU_all{i}{m,s}(array_ind_subwindow) = zU; %anemometer height - from 30-minute window
                    zs_all{i}{m,s}(array_ind_subwindow) = zs; %roughness height - from 30-minute window
                    uth_all{i}{m,s}(array_ind_subwindow) = uth; %threshold wind from flux frequency
                    tauth_all{i}{m,s}(array_ind_subwindow) = tauth; %threshold wind from flux frequency
                    ustth_all{i}{m,s}(array_ind_subwindow) = ustth; %threshold wind from flux frequency
                    timeofday_all{i}{m,s}(array_ind_subwindow) = timeofday; %time of day in hours
                    starttime_all{i}{m,s}(array_ind_subwindow) = StartTimes_Deltat(k); %start time
                    endtime_all{i}{m,s}(array_ind_subwindow) = EndTimes_Deltat(k); %start time
                    t_all{i}{m,s}{array_ind_subwindow} = t_subwindow; %timeseries times array
                    u_all{i}{m,s}{array_ind_subwindow} = u_subwindow; %wind timeseries array
                    n_all{i}{m,s}{array_ind_subwindow} = n_subwindow; %counts timeseries array
                end
            end
        end
    end
end

%% save analysis data
save(SaveData_Path,'SiteNames','Sites','N_Sites','*all'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create sample timeseries plot %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% info about time scales for sample plot
deltat_analysis = duration(0,0,1); %sampling interval for analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
ind_Deltat = find(Deltat_all==Deltat_analysis); %get measurement interval
ind_deltat = find(deltat_all==deltat_analysis); %get sampling interval

%% info about site and times for sample plot
Site_sampleplot = 'RanchoGuadalupe';
starttime_sampleplot = datetime(2015,3,24,13,31,0);
endtime_sampleplot = datetime(2015,3,24,13,36,0);
T_sampleplot = seconds(endtime_sampleplot-starttime_sampleplot); %duration of sample plot

%% get indices for windows, subwindows, and site
ind_Site = find(strcmp(Sites,Site_sampleplot)); %index for site
ind_window = min(find(StartTime_window{ind_Site}>=starttime_sampleplot)); %index for window
flux_ind_window = intersect(find(t_flux_int_window{ind_Site}{ind_window}>=starttime_sampleplot),...
    find(t_flux_int_window{ind_Site}{ind_window}<=endtime_sampleplot));
wind_ind_window = intersect(find(t_wind_int_window{ind_Site}{ind_window}>=starttime_sampleplot),...
    find(t_wind_int_window{ind_Site}{ind_window}<=endtime_sampleplot));
starttimes_subwindows = starttime_all{ind_Site}{ind_Deltat,ind_deltat}; %start times of contained subwindows
endtimes_subwindows = endtime_all{ind_Site}{ind_Deltat,ind_deltat}; %end times for contained subwindows
ind_subwindows = (min(find(starttimes_subwindows>=starttime_sampleplot))):...
    max((find(endtimes_subwindows<=endtime_sampleplot))); %indices for contained subwindows
N_subwindows = length(ind_subwindows); %get number of subwindows

%% get raw values for plotting
t_flux_raw = seconds(t_flux_int_window{ind_Site}{ind_window}(flux_ind_window)-starttime_sampleplot); %flux times
n_raw = sum(n_int_window{ind_Site}{ind_window}(flux_ind_window,:),2); %counts rate
t_wind_raw = seconds(t_wind_int_window{ind_Site}{ind_window}(wind_ind_window)-starttime_sampleplot); %wind times
u_raw = u_int_window{ind_Site}{ind_window}(wind_ind_window,:); %wind speed

%% get sample-averaged times series and measurement interval values
t_sample = []; %times
n_sample = []; %counts
u_sample = []; %wind speeds
tstart_measurement = []; %start times for measurement interval subwindows
tend_measurement = []; %start times for measurement interval subwindows
for k = ind_subwindows
    t_measurement = seconds(t_all{ind_Site}{ind_Deltat,ind_deltat}{k}-starttime_sampleplot);
    tstart_measurement = [tstart_measurement; min(t_measurement)];
    tend_measurement = [tend_measurement; max(t_measurement)];
    t_sample = [t_sample; t_measurement];
    n_sample = [n_sample; n_all{ind_Site}{ind_Deltat,ind_deltat}{k}];
    u_sample = [u_sample; u_all{ind_Site}{ind_Deltat,ind_deltat}{k}];
end
tmid_measurement = mean([tstart_measurement,tend_measurement],2); %midpoint of each measurement interval
fQ_measurement = fQ_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindows); %activities for each measurement interval
uth_measurement = uth_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindows); %thresholds for each measurement interval

%% start plotting
figure(1); clf; %initialize figure
PlotFont = 14; %set plotting font

%plot raw flux
subplot(2,3,1); cla
plot(t_flux_raw,n_raw);
xlim([0 T_sampleplot]); %set plot xlimits
ylim_n = ylim; %get plot y-limits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('counts, $$N$$ (s$$^{-1}$$)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_n(2),'(a)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot raw wind
subplot(2,3,4); cla
plot(t_wind_raw,u_raw);
xlim([0 T_sampleplot]); %set plot xlimits
ylim_u = [floor(min(u_raw)),ceil(max(u_raw))]; %get plot y-limits
ylim(ylim_u); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_u(2),'(b)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot sample-averaged flux
subplot(2,3,2); cla; hold on;
plot(t_sample,n_sample);
for k = 2:N_subwindows;
    plot(tstart_measurement(k)*[1,1],ylim_n,'k--');
end
xlim([0 T_sampleplot]); %set plot xlimits
ylim(ylim_n); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('counts, $$N$$ (s$$^{-1}$$)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_n(2),'(c)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot sample-averaged wind
subplot(2,3,5); cla; hold on;
plot(t_sample,u_sample);
for k = 2:N_subwindows;
    plot(tstart_measurement(k)*[1,1],ylim_u,'k--');
end
xlim([0 T_sampleplot]); %set plot xlimits
ylim(ylim_u); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_u(2),'(d)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot transport activity
subplot(2,3,3); cla; hold on;
bar(tmid_measurement,fQ_measurement,'FaceColor','c','EdgeColor','b');
xlim([0 T_sampleplot]);
ylim([0 1]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('transport activity, $$f_Q$$', 'Interpreter','Latex');
text(T_sampleplot/30,0.95,'(e)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot threshold wind
subplot(2,3,6); cla; hold on;
bar(tmid_measurement,uth_measurement,'FaceColor','c','EdgeColor','b');
xlim([0 T_sampleplot]);
ylim(ylim_u);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('threshold wind, $$u_{th}$$ (m/s)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_u(2),'(f)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
print([folder_Plots,'effectivethreshold_sample.png'],'-dpng');