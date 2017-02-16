%% SCRIPT TO PLOT VALUES FROM DATA WINDOWS

%% initialize
clearvars;

%% parameters for sample plots
dt_plot = duration(0,0,5); %time for window-average for plot

%% folders for loading and saving data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

%% paths for loading and plotting data - restricted
LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Restricted'); %path for 30 minute data - for flux law analysis
folder_TimeseriesPlot = '../../PlotOutput/Timeseries/30min_Restricted/'; %folder for sample timeseries plots

% %% paths for loading and plotting data - unrestricted
% LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis
% folder_TimeseriesPlot = '../../PlotOutput/Timeseries/30min_Unrestricted/'; %folder for sample timeseries plots

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize sample plot
h1 = figure(1);
set(h1,'visible','off');

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
    
    WindowStartTimes = StartTime_window{i};
    WindowEndTimes = EndTime_window{i};
    N_Windows = length(WindowStartTimes);
        
    %% go through time windows
    for j = 1:N_Windows

        %% display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]

        %% get time information
        StartTime = WindowStartTimes(j); %get specific start time
        EndTime = WindowEndTimes(j); %get specific end time
        T_interval = seconds(EndTime-StartTime); % get duration of interval in seconds
      
        %% get timeseries info
        t_wind = t_wind_int_window{i}{j}; %wind times
        u_int = u_int_window{i}{j}; %interpolated streamwise wind
        w_int = w_int_window{i}{j}; %interpolated vertical wind
        t_flux = t_flux_int_window{i}{j}; %flux times (including interpolation)
        q_int = q_int_window{i}{j}; %partial flux
        zW = zW_window{i}{j}; %Wenglor heights
        N_zW = length(zW); %number of Wenglor heights
        
        %% get window averages
        t_plot = t_wind(1):dt_plot:t_wind(end-1); %times
        u_plot = window_average(u_int, t_wind, dt_plot); %compute window average
        w_plot = window_average(w_int, t_wind, dt_plot); %compute window average
        q_plot = zeros(length(t_plot),N_zW);
        for k=1:N_zW
            if ~isnan(q_int(:,k))
                q_plot(:,k) = window_average(q_int(:,k), t_flux, dt_plot); %compute window average
            else
                q_plot(:,k) = zeros(size(t_plot)); %if empty, set to zero
            end
        end

        %% create plot of streamwise wind
        subplot(6,1,1:2);
        plot(t_plot,u_plot);
        ylabel('streamwise wind, $$u$$ (ms$$^{-1}$$)','Interpreter','Latex');
        title([SiteNames{i},', ',datestr(StartTime, 'yyyy-mm-dd HH:MM'),' - ',datestr(EndTime, 'HH:MM')]);
        set(gca,'FontSize',14,'XMinorTick','On','YMinorTick','On','Box','On');
        ylims = [floor(min(u_plot)) ceil(max(u_plot))];
        ylim(ylims);
        %text(min(t_plot)+duration(0,0,30),ylims(2)*0.97,'(a)','FontSize',14);

        %% create plot of vertical wind
        subplot(6,1,3);
        plot(t_plot,w_plot);
        ylabel('vert. wind, $$w$$ (ms$$^{-1}$$)','Interpreter','Latex');
        set(gca,'FontSize',14,'XMinorTick','On','YMinorTick','On','Box','On');
        ylims = [floor(min(w_plot)*10) ceil(max(w_plot)*10)]/10;
        ylim(ylims);
        %text(min(t_plot)+duration(0,0,30),ylims(2)*0.65,'(b)','FontSize',14);

        %% create plot of partial sediment flux
        subplot(6,1,4:6); cla; hold on;
        legend_items = cell(N_zW,1);
        for k=1:N_zW
            plot(t_plot,q_plot(:,k));
            legend_items{k} = ['$$z_',int2str(k),'$$ = ',num2str(round(zW(k)*100,1)),' cm'];
        end
        legend(legend_items,'Location','NorthWest','Interpreter','Latex');
        ylabel('partial salt. flux, $$q_i$$ (gm$$^{-2}$$s$$^{-1}$$)','Interpreter','Latex');
        set(gca,'FontSize',14,'XMinorTick','On','YMinorTick','On','Box','On');
        ylims = [0 ceil(max(max(q_plot))/10)]*10;
        if range(ylims)==0; ylims = [0 10]; end; %adjust ylims if no saltation
        ylim(ylims);
        %text(min(t_plot)+duration(0,7,0),ylims(2)*0.95,'(c)','FontSize',14);
        
        %% print plots
        set(gcf, 'PaperPosition',[0 0 9 11]);
        plotname = ['SampleTimeseries_',Sites{i},'_',int2str(j),'_',datestr(StartTime,'yyyymmdd_HHMM'),'-',datestr(EndTime,'HHMM'),'.png'];
        print([folder_TimeseriesPlot,plotname],'-dpng');
    end
end