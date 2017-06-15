% %% INITIALIZATION
% %initialize
% clearvars;
% close all;
% 
% %% PARAMETERS AND INPUTS
% Site_SamplePlot = 'Oceano';
% StartTime_SamplePlot = datetime(2015,5,16,12,40,0);
% EndTime_SamplePlot = datetime(2015,5,16,12,45,0);
% T_subwindow_SamplePlot = duration(0,0,1);
% 
% %% INFO FOR PLOTTING
% PlotFont = 14;
% 
% %% LOAD DATA AND FUNCTIONS
% %folders for loading data, saving data, and functions
% folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
% folder_Plots = '../../PlotOutput/WenglorCalibration/'; %folder containing plot output
% folder_Functions = '../Functions/'; %folder with functions
% 
% %paths for loading and saving data
% LoadData_Path = strcat(folder_LoadData,'DataFullSubwindowCalcs_30min_Restricted'); %path for 30 minute data
% 
% %load data
% load(LoadData_Path); %load window data
% 
% %load functions
% addpath(folder_Functions); %point MATLAB to location of functions
% 
% %% GET INFO FOR SAMPLE PLOT
% i = find(strcmp(Sites,Site_SamplePlot));
% m = find(T_subwindow == T_subwindow_SamplePlot);
% j = intersect(...
%     find(StartTime_subwindow{i}{m}>=StartTime_SamplePlot),...
%     find(StartTime_subwindow{i}{m}<=EndTime_SamplePlot));
% 
% %% GET VALUES FOR SAMPLE PLOT
% u = ubar_subwindow{i}{m}(j);
% t = StartTime_subwindow{i}{m}(j);
% Qsum = Qsum_subwindow{i}{m}(j);
% Qfit = Qfit_subwindow{i}{m}(j);

%% PLOT VALUES
figure(1); clf;

%initialize subplot for wind speeds
subplot('Position',[0.11 0.52 0.84 0.45]); hold on;

%plot wind speed
plot(t,u,'k','LineWidth',1);

% format plot
%title([Site_SamplePlot,', ',datestr(StartTime_SamplePlot,'yyyy-mm-dd HH:MM'),'-',datestr(EndTime_SamplePlot,'HH:MM')])
ylabel('Horizontal wind speed, $$u$$ (m s$$^{-1}$$)', 'Interpreter','Latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XTickLabel',[]);
set(gca,'FontSize',PlotFont);
ylims = ylim;
text(t(10),ylims(1)+0.95*range(ylims),'(a)','FontSize',PlotFont);

%initialize subplot for saltation fluxes
subplot('Position',[0.11 0.04 0.84 0.45]); hold on;

%plot saltation flux
h1 = plot(t,Qsum,'b','LineWidth',2)
h2 = plot(t,Qfit,'r','LineWidth',1);
set(gca,'FontSize',PlotFont);

% format plot
ylabel('Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','Latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);
legend([h2, h1],{'fitting method, Q_{fit}','summation method, Q_{sum}'},'Location','North');
ylims = ylim;
ylim([0 ylims(2)]);
text(t(10),ylims(1)+0.95*range(ylims),'(b)','FontSize',PlotFont);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 8],'PaperPosition',[0 0 6 8],'PaperPositionMode','Manual');
print([folder_Plots,'SampleTimeseries_TotalFlux.png'],'-dpng');