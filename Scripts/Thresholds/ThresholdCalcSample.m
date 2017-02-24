%% SCRIPT TO CREATE SAMPLE PLOT FOR THRESHOLD CALCULATION

%% initialize
clearvars;

%% info about site and times for sample plot
Site_sampleplot = 'RanchoGuadalupe';
StartTime_sampleplot = datetime(2015,3,24,13,31,0);
EndTime_sampleplot = datetime(2015,3,24,13,36,0);

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,2); %sampling interval for analysis

%% plotting info
PlotFont = 12;

%% folders for loading and saving data
folder_LoadFluxLawData = '../../AnalysisData/FluxLaw/'; %folder for loading 30 minute data
folder_LoadThresholdData = '../../AnalysisData/Windowing/'; %folder for retrieving subwindow calcs
folder_LoadSubwindowData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving subwindow values
folder_LoadWindowData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving subwindow values
folder_LoadRoughnessData = '../../AnalysisData/Roughness/'; %folder for loading roughness data
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% paths for loading and saving data
LoadFluxLawData_Path = strcat(folder_LoadFluxLawData,'FluxLawCalcs_30min_Restricted'); %path for 30 minute data
LoadThresholdData_Path = strcat(folder_LoadThresholdData,'DataSubwindowCalcs_30min_Unrestricted'); %path for loading subwindow calcs
LoadWindowData_Path = strcat(folder_LoadWindowData,'DataWindows_30min_Unrestricted'); %path for loading subwindow data
LoadSubwindowData_Path = strcat(folder_LoadSubwindowData,'DataSubwindows_30min_Unrestricted'); %path for loading subwindow data
LoadRoughnessData_Path = strcat(folder_LoadRoughnessData,'RoughnessCalcs_30min_Restricted'); %path for 30 minute data

%% load data and rename variables as necessary
load(LoadFluxLawData_Path); %load 30-minute values
load(LoadThresholdData_Path); %load primary data
load(LoadWindowData_Path); %load 30-minute values
load(LoadSubwindowData_Path); %load 30-minute values
load(LoadRoughnessData_Path); %load Roughness values
addpath(folder_Functions); %point MATLAB to location of functions

%% get indices for plotting
ind_Deltat = find(Deltat_all==Deltat_analysis); %get measurement interval
ind_deltat = find(deltat_all==deltat_analysis); %get sampling interval
ind_Site = find(strcmp(Sites,Site_sampleplot)); %index for site

%% get z0 for Sites
z0_Site = z0Re_Q_fit;
sigma_z0_Site = 10.^(sigma_z0Re_ln_Q_fit);

%% get flux law intercept values for Sites
tauit_intercept = tauit_linearfit_all;
sigma_tauit_intercept = sigma_tauit_linearfit_all;

%% get duration of sample plot
T_sampleplot = seconds(EndTime_sampleplot-StartTime_sampleplot); %duration of sample plot

%% for plots, get window times associated with specific site
StartTime_window_analysis = StartTime_window{ind_Site};
EndTime_window_analysis = EndTime_window{ind_Site};
t_flux_window_analysis = t_flux_window{ind_Site};
ntotal_window_analysis = ntotal_window{ind_Site};
t_wind_window_analysis = t_wind_window{ind_Site};
u_window_analysis = u_window{ind_Site};

%% for plots, get subwindow times associated with specific site
StartTime_subwindow_analysis = StartTime_subwindow{ind_Site}{ind_Deltat,ind_deltat};
EndTime_subwindow_analysis = EndTime_subwindow{ind_Site}{ind_Deltat,ind_deltat};

%% get raw wind and flux data for sample plot
ind_window = max(find(StartTime_window_analysis<=StartTime_sampleplot)):min(find(EndTime_window_analysis>=EndTime_sampleplot));
t_flux_sampleplot = seconds(t_flux_window_analysis{ind_window}-StartTime_sampleplot);
ntotal_sampleplot = ntotal_window_analysis{ind_window};
t_wind_sampleplot = seconds(t_wind_window_analysis{ind_window}-StartTime_sampleplot);
u_sampleplot = u_window_analysis{ind_window};

%% pick out specific times for plotting
ind_flux_sampleplot = intersect(find(t_flux_sampleplot>=0),find(t_flux_sampleplot<=T_sampleplot)); % restrict to times in sample plot window
t_flux_sampleplot = t_flux_sampleplot(ind_flux_sampleplot); %restricted t
ntotal_sampleplot = ntotal_sampleplot(ind_flux_sampleplot); %restricted ntotal
ind_wind_sampleplot = intersect(find(t_wind_sampleplot>=0),find(t_wind_sampleplot<=T_sampleplot)); % restrict to times in sample plot window
t_wind_sampleplot = t_wind_sampleplot(ind_wind_sampleplot); %restricted t
u_sampleplot = u_sampleplot(ind_wind_sampleplot); %resticted u

%% get window-avg data for sample plot
ind_subwindows = min(find(StartTime_subwindow_analysis>=StartTime_sampleplot)):max(find(EndTime_subwindow_analysis<=EndTime_sampleplot));
N_subwindows = length(ind_subwindows); %get number of subwindows
t_avg_sampleplot = [];
ntotal_avg_sampleplot = [];
u_avg_sampleplot = [];
for i = 1:N_subwindows
    t_avg_sampleplot = [t_avg_sampleplot; seconds(t_subwindow{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindows(i)} - StartTime_sampleplot)];
    ntotal_avg_sampleplot = [ntotal_avg_sampleplot; ntotal_subwindow{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindows(i)}];
    u_avg_sampleplot = [u_avg_sampleplot; u_subwindow{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindows(i)}];
end
ind_avg_sampleplot = intersect(find(t_avg_sampleplot>=0),find(t_avg_sampleplot<=T_sampleplot)); % restrict to times in sample plot window
t_avg_sampleplot = t_avg_sampleplot(ind_avg_sampleplot); %restricted t
ntotal_avg_sampleplot = ntotal_avg_sampleplot(ind_avg_sampleplot); %restricted ntotal
u_avg_sampleplot = u_avg_sampleplot(ind_avg_sampleplot); %resticted u
ylim_ntotal_sampleplot = ceil(max(ntotal_sampleplot)/10)*10; %get maximum ntotal for sample plot

%% get measurement interval start, end, and midpoint times
tstart_measurement_sampleplot = seconds(StartTime_subwindow_analysis(ind_subwindows) - StartTime_sampleplot);
tend_measurement_sampleplot = seconds(EndTime_subwindow_analysis(ind_subwindows) - StartTime_sampleplot);
tmid_measurement_sampleplot = (tend_measurement_sampleplot+tstart_measurement_sampleplot)/2;

%% get threshold info for measurement intervals
fQ_measurement_sampleplot = fQ_subwindow_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindows);
uth_measurement_sampleplot = uth_subwindow_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindows);

%% start plotting
figure(1); clf; %initialize figure

%plot raw flux
%subplot(2,3,1); cla
subplot('Position',[0.07 0.6 0.245 0.38]); hold on;
plot(t_flux_sampleplot,ntotal_sampleplot,'b');
xlim([0 T_sampleplot]); %set plot xlimits
ylim([0 ylim_ntotal_sampleplot]); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('counts, $$N$$ (s$$^{-1}$$)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_ntotal_sampleplot,'(a)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot raw wind
%subplot(2,3,4); cla
subplot('Position',[0.07 0.10 0.245 0.38]); hold on;
plot(t_wind_sampleplot,u_sampleplot,'b');
xlim([0 T_sampleplot]); %set plot xlimits
ylim_u = [floor(min(u_sampleplot)),ceil(max(u_sampleplot))]; %get plot y-limits
ylim(ylim_u); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_u(2),'(d)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot sample-averaged flux
%subplot(2,3,2); cla; hold on;
subplot('Position',[0.4 0.6 0.245 0.38]); hold on;
plot(t_avg_sampleplot,ntotal_avg_sampleplot,'b');
for k = 2:N_subwindows
    plot(tstart_measurement_sampleplot(k)*[1,1],[0 ylim_ntotal_sampleplot],'k--');
end
xlim([0 T_sampleplot]); %set plot xlimits
ylim([0 ylim_ntotal_sampleplot]); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('counts, $$N$$ (s$$^{-1}$$)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_ntotal_sampleplot,'(b)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot sample-averaged wind
subplot(2,3,5); cla; hold on;
subplot('Position',[0.4 0.10 0.245 0.38]); hold on;
plot(t_avg_sampleplot,u_avg_sampleplot,'b');
for k = 2:N_subwindows
    plot(tstart_measurement_sampleplot(k)*[1,1],ylim_u,'k--');
end
xlim([0 T_sampleplot]); %set plot xlimits
ylim(ylim_u); %set plot ylimits
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_u(2),'(e)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot transport activity
subplot(2,3,3); cla; hold on;
subplot('Position',[0.735 0.6 0.245 0.38]); hold on;
bar(tmid_measurement_sampleplot,fQ_measurement_sampleplot,'FaceColor','c','EdgeColor','b');
xlim([0 T_sampleplot]);
ylim([0 1]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('transport activity, $$f_Q$$', 'Interpreter','Latex');
text(T_sampleplot/30,0.95,'(c)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%plot threshold wind
subplot(2,3,6); cla; hold on;
subplot('Position',[0.735 0.10 0.245 0.38]); hold on;
bar(tmid_measurement_sampleplot,uth_measurement_sampleplot,'FaceColor','c','EdgeColor','b');
xlim([0 T_sampleplot]);
ylim(ylim_u);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('time, $$t$$ (s)','Interpreter','Latex');
ylabel('threshold wind, $$u_{th}$$ (m/s)','Interpreter','Latex');
text(T_sampleplot/30,0.95*ylim_u(2),'(f)','FontSize',PlotFont)
set(gca, 'FontSize', PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[9 5],'PaperPosition',[0 0 9 5],'PaperPositionMode','Manual');
print([folder_Plots,'ThresholdCalcSample.png'],'-dpng');
print([folder_Plots,'ThresholdCalcSample.tif'],'-dtiff','-r600');