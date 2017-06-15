%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;


%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

%% info about publication plot
StartTime_PublicationPlot = datetime(2015,3,24,14,40,0);

%% paths for loading and saving data - restricted
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data - for flux law analysis
folder_WenglorPlots = '../../PlotOutput/WenglorFluxProfiles_Restricted/'; %folder for plots

% %% paths for loading and saving data - unrestricted
% LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis
% folder_WenglorPlots = '../../PlotOutput/WenglorFluxProfiles_Unrestricted/'; %folder for plots

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%set info about plot sizing
PlotFont = 12;
MarkerSize_plot = 6;
LineWidth_plot = 0.5;

%initialize profile plots
close all;
h1 = figure(1);
set(h1,'visible','off');
set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On','YScale','log');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('Wenglor height, $$z$$ (m)','Interpreter','Latex');
ylabel('Height-specific saltation flux, $$q$$ (g m$$^{-2}$$ s$$^{-1}$$)','Interpreter','Latex');
set(gcf,'PaperUnits','inches','PaperSize',[5 4.5],'PaperPosition',[0 0 5 4.5],'PaperPositionMode','Manual');

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% get number of windows
    WindowStartTimes = StartTimes_all{i};
    WindowEndTimes = EndTimes_all{i};
    N_Windows = length(WindowStartTimes);
        
    %% go through time blocks
    for j = 1:N_Windows

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]

        %get specific start and end time
        StartTime = WindowStartTimes(j);
        EndTime = WindowEndTimes(j);

        %% FLUX CALCULATIONS FOR INTERVAL

        %get total flux
        Q = Q_all{i}(j);
        
        %only perform plot if non-zero flux
        if Q>0

            %get Wenglor heights
            zW = zW_all{i}{j};
            sigma_zW = sigma_zW_all{i}{j};
            
            %get Wenglor profile
            qbar = qbar_all{i}{j}; %get mean qz profile
            sigma_qbar = sigma_qbar_all{i}{j}; %get uncertainty in qz profile
            Q = Q_all{i}(j); %get total flux
            zq = zq_all{i}(j); %get flux height
            
            %get predicted q profile
            q_pred = (Q/zq)*exp(-zW/zq);

            %get min and max q
            q_min = max([min([min(qbar-sigma_qbar),min(q_pred)]),0]); %minimum q (or zero)
            q_min_plot = 10^(floor(log10(q_min)*5)/5); %set minimum of q range to nearest fifth of decade below
            q_max = max([max(qbar+sigma_qbar),max(q_pred)]); %maximum q
            q_max_plot = 10^(ceil(log10(q_max)*5)/5); %set maximum of q range to nearest fifth of decade above
            
            %plot flux profile fit
            cla; hold on;
            errorbar(zW,qbar,sigma_qbar,'b+','MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);
            plot(zW,q_pred,'k');
            ylim([q_min_plot q_max_plot]);
            legend('data','fit','Location','NorthEast');
            title([SiteNames{i},', ',datestr(StartTime, 'yyyy-mm-dd HH:MM'),' - ',datestr(EndTime, 'HH:MM')]);
            
            %print plot
            plotname = ['WenglorFluxProfile_',Sites{i},'_',int2str(j),'_',datestr(StartTime,'yyyymmdd_HHMM'),'-',datestr(EndTime,'HHMM')];
            print([folder_WenglorPlots,plotname,'.png'],'-dpng'); %draft plot
            if StartTime == StartTime_PublicationPlot
                print([folder_WenglorPlots,plotname,'.tif'],'-dtiff','-r600'); %publication plot (for selected time only)
            end
        end
    end
end