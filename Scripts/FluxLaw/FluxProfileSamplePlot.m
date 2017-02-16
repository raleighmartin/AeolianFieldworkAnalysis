%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;
close all;

%% information for sample plot
Site_sampleplot = 'RanchoGuadalupe';
StartTime_sampleplot = datetime(2015,3,24,14,40,0);

%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for retrieving data for this analysis
folder_SaveData = '../../AnalysisData/FluxLaw/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/FluxLaw/'; %folder for plots

%% paths for loading and saving data - restricted
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data

%load data
load(LoadData_Path);
addpath(folder_Functions); %point MATLAB to location of functions

%% get information about plotting
ind_Site = find(strcmp(Site_sampleplot,Sites));
ind_Window = find(StartTimes_all{ind_Site}==StartTime_sampleplot);
StartTime = StartTimes_all{ind_Site}(ind_Window);
EndTime = EndTimes_all{ind_Site}(ind_Window);
zW = zW_all{ind_Site}{ind_Window};
qbar = qbar_all{ind_Site}{ind_Window};
sigma_qbar = sigma_qbar_all{ind_Site}{ind_Window};
zq = zq_all{ind_Site}(ind_Window);
q0 = Q_all{ind_Site}(ind_Window)/zq;
qpred = q0*exp(-zW/zq);

%% plot flux profile fit
figure(1); hold on;

%plotting
errorbar(zW,qbar,sigma_qbar,'b+','MarkerSize',10);
plot(zW,qpred,'k');

%plotting parameters
set(gca,'FontSize',16,'XMinorTick','On','YMinorTick','On','Box','On','YScale','log');
xlabel('Wenglor height, $$z_i$$ (m)','Interpreter','Latex');
ylabel('30-minute partial saltation flux, $$q_i$$ (g m$$^{-2}$$ s$$^{-1}$$)','Interpreter','Latex');
legend('data','fit','Location','NorthEast');
title([Site_sampleplot,', ',datestr(StartTime, 'yyyy-mm-dd HH:MM'),' - ',datestr(EndTime, 'HH:MM')]);

%print plot
set(gcf,'PaperPosition',[0 0 8 6]);
print([folder_Plots,'WenglorFluxProfile_',Site_sampleplot,'_',...
    datestr(StartTime,'yyyymmdd_HHMM'),'_',datestr(EndTime,'HHMM'),...
    '.png'],'-dpng');