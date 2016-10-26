%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_WenglorPlots = '../../PlotOutput/WenglorFluxProfiles_Unrestricted/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%% paths for loading and saving data - unrestricted
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis

% %% paths for loading and saving data - restricted
% LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data - for thresholds analysis

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%initialize profile plots
close all;
h1 = figure(1);
set(h1,'visible','off');
set(gca,'FontSize',16,'XMinorTick','On','YMinorTick','On','Box','On','YScale','log');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('Wenglor height, $$z_i$$ (m)','Interpreter','Latex');
ylabel('30-minute partial saltation flux, $$\tilde{q}_i$$ (g m$$^{-2}$$ s$$^{-1}$$)','Interpreter','Latex');
set(gcf, 'PaperPosition',[0 0 8 6]);

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
        
        %only perform plot if there is flux data in window
        if ~isnan(Q)

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

            %plot flux profile fit
            figure(1);
            cla; hold on;
            errorbar(zW,qbar,sigma_qbar,'b+','MarkerSize',10);
            plot(zW,q_pred,'k');
            legend('data','fit','Location','NorthEast');
            title([SiteNames{i},', ',datestr(StartTime, 'yyyy-mm-dd HH:MM'),' - ',datestr(EndTime, 'HH:MM')]);
            print([folder_WenglorPlots,'WenglorFluxProfile_',Sites{i},'_',int2str(j),'.png'],'-dpng');
        end
    end
end