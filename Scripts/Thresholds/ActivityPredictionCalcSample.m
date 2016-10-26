%% SCRIPT TO COMPARE ACTIVITY PREDICTION TO OBSERVATION

%% initialize
clearvars;

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,1); %sampling interval for analysis

%% info for sample plot
%Site_sampleplot = 'Jericoacoara';
%Site_sampleplot = 'RanchoGuadalupe';
Site_sampleplot = 'Oceano';
fQ_sampleplot = [0.2, 0.5, 0.8]; %fQ for timeseries sample plot

%% plotting info
PlotFont = 14;

%% folders for loading and saving data
folder_LoadSubwindowData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_LoadSubwindowCalcs = '../../AnalysisData/Windowing/'; %folder for retrieving subwindow calcs
folder_LoadThresholdData = '../../AnalysisData/Thresholds/'; %folder for loading threshold values
folder_LoadActivityPredictionData = '../../AnalysisData/Thresholds/'; %folder for loading activity prediction
folder_Functions = '../Functions/'; %folder with functions
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% paths for loading and saving data
LoadSubwindowData_Path = strcat(folder_LoadSubwindowData,'DataSubwindows_30min_Unrestricted'); %path for loading time window data
LoadSubwindowCalcs_Path = strcat(folder_LoadSubwindowCalcs,'DataSubwindowCalcs_30min_Unrestricted'); %path for loading time window data
LoadThresholdData_Path = strcat(folder_LoadThresholdData,'ThresholdAnalysisData'); %path for threshold data
LoadActivityPredictionData_Path = strcat(folder_SaveData,'ActivityPredictionData'); %folder for loading activity prediction
SaveData_Path = strcat(folder_SaveData,'ActivityPredictionData'); %path for saving output data

%% load data
load(LoadSubwindowData_Path); %load primary data
load(LoadSubwindowCalcs_Path); %load subwindow values
load(LoadThresholdData_Path); %load threshold values
load(LoadActivityPredictionData_Path); %load fQ prediction values
addpath(folder_Functions); %point MATLAB to location of functions

%% information about measurement and sampling interval indices for analysis
ind_Deltat = find(Deltat_all == Deltat_analysis); %index for measurement interval
ind_deltat = find(deltat_all == deltat_analysis); %index for sampling interval
ind_Site = find(strcmp(Sites,Site_sampleplot)); %index for site

%%%%%%%%%
% PLOTS %
%%%%%%%%%

%initialize lists for sample plot
N_fQ_sampleplot = length(fQ_sampleplot); %get number of fQ's for plotting
ind_fQ_sampleplot = zeros(N_fQ_sampleplot,1); %initialize indices of subwindows for plotting
t_sampleplot = cell(N_fQ_sampleplot,1); %initialize t timeseries for sampleplot
u_sampleplot = cell(N_fQ_sampleplot,1); %initialize u timeseries for sampleplot
u_min_sampleplot = zeros(N_fQ_sampleplot,1); %initialize mininum u for each sampleplot
u_max_sampleplot = zeros(N_fQ_sampleplot,1); %initialize maximum u for each sampleplot
uft_sampleplot = zeros(N_fQ_sampleplot,1); %initialize uit's for sampleplot
uit_sampleplot = zeros(N_fQ_sampleplot,1); %initialize uit's for sampleplot
ind_ft_sampleplot = cell(N_fQ_sampleplot,1); %initialize indices of activity for fluid threshold hypothesis
ind_it_sampleplot = cell(N_fQ_sampleplot,1); %initialize indices of activity for impact threshold hypothesis
ind_hyst_sampleplot = cell(N_fQ_sampleplot,1); %initialize indices of activity for hysteresis hypothesis

%get values for sample plot
StartTime = StartTime_subwindow{ind_Site}{ind_Deltat,ind_deltat};
EndTime = EndTime_subwindow{ind_Site}{ind_Deltat,ind_deltat};
fQ = fQ_subwindow_all{ind_Site}{ind_Deltat,ind_deltat};
u = u_subwindow{ind_Site}{ind_Deltat,ind_deltat};
t = t_subwindow{ind_Site}{ind_Deltat,ind_deltat};
zU = zU_subwindow{ind_Site}{ind_Deltat,ind_deltat};
ustit = sqrt(tauit_all(ind_Site)/rho_Site(ind_Site));
uit = (ustit/kappa).*log(zU./z0_Site(ind_Site));
ustft = sqrt(tauft_all(ind_Site)/rho_Site(ind_Site));
uft = (ustft/kappa).*log(zU./z0_Site(ind_Site));
N_subwindows = length(StartTime);

for j = 1:N_fQ_sampleplot
    ind_fQ_sampleplot_all = find(abs(fQ_sampleplot(j)-fQ)==min(abs(fQ_sampleplot(j)-fQ))); %get all possible indices for fQ for sample plot
    fQ_diff = abs(fQpred_hyst_all{ind_Site}(ind_fQ_sampleplot_all)-fQ_sampleplot(j)); %get the differences between predicted fQs for these indices and sample value
    ind_fQ_sampleplot(j) = ind_fQ_sampleplot_all(find(fQ_diff==min(fQ_diff),1)); %get the index of the fQpred that is closest to fQ for sample plot
    t_datetime = t{ind_fQ_sampleplot(j)}; %get times for sample plot (in datetime format)
    t_sampleplot{j} = seconds(t_datetime-t_datetime(1)+deltat_analysis/2); %get times for sample plot (in seconds)
    u_sampleplot{j} = u{ind_fQ_sampleplot(j)}; %get wind speeds for sample plot
    u_min_sampleplot(j) = min(u_sampleplot{j}); %get mininum wind speed
    u_max_sampleplot(j) = max(u_sampleplot{j}); %get maximum wind speed
    uft_sampleplot(j) = uft(ind_fQ_sampleplot(j)); %get fluid threshold
    uit_sampleplot(j) = uit(ind_fQ_sampleplot(j)); %get impact threshold
    ind_ft_sampleplot{j} = ind_fplus_all{ind_Site}{ind_fQ_sampleplot(j)}; %get indices of activity for fluid threshold hypothesis
    ind_it_sampleplot{j} = union(ind_fplus_all{ind_Site}{ind_fQ_sampleplot(j)},ind_fint_all{ind_Site}{ind_fQ_sampleplot(j)}); %get indices of activity for impact threshold hypothesis
    ind_hyst_sampleplot{j} = union(ind_fplus_all{ind_Site}{ind_fQ_sampleplot(j)},ind_fint_down_all{ind_Site}{ind_fQ_sampleplot(j)}); %get indices of activity for hysteresis hypothesis
end

%generate plot
figure(1); clf;
paneltext = {'(a)','(b)','(c)';'(d)','(e)','(f)'}; %panel labels

xlims_sampleplot = [0 seconds(Deltat_analysis)];
ylims_sampleplot = [floor(min(u_min_sampleplot)), ceil(max(u_max_sampleplot))];

for j = 1:N_fQ_sampleplot
    %plot timeseries
    subplot(2,N_fQ_sampleplot,j); hold on;
    plot(t_sampleplot{j},u_sampleplot{j},'k','LineWidth',1);
    plot(xlims_sampleplot,uft_sampleplot(j)*[1 1],'b-.','LineWidth',2);
    plot(xlims_sampleplot,uit_sampleplot(j)*[1 1],'r--','LineWidth',2);
    plot(t_sampleplot{j},u_sampleplot{j},'k','LineWidth',1);
    text(2,ylims_sampleplot(2)-0.5,paneltext{1,j},'FontSize',PlotFont)
    xlim(xlims_sampleplot);
    ylim(ylims_sampleplot);
    xlabel('time (s)');
    if j==1
        ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
        h_legend = legend('wind speed','fluid thres., u_{ft}','impact thres., u_{it}');
        set(h_legend,'FontSize',PlotFont-4,'Location','NorthEast');
    end
    title(['$$f_{Q} = ',num2str(fQ_sampleplot(j),1),'$$'],'Interpreter','Latex');
    set(gca,'XMinorTick','On','YTick',ylims_sampleplot(1):ylims_sampleplot(2),'YMinorTick','On','Box','On','FontSize',PlotFont);

    %get info for area plot
    t_area = [];
    y_it = [];
    y_ft = [];
    y_hyst = [];
    for k = 1:length(t_sampleplot{j})
        t_area = [t_area; t_sampleplot{j}(k)-seconds(deltat_analysis)/2; t_sampleplot{j}(k)+seconds(deltat_analysis)/2];
        if ~isempty(intersect(k,ind_it_sampleplot{j}))
            y_it = [y_it; 1; 1];
        else
            y_it = [y_it; 0; 0];
        end

        if ~isempty(intersect(k,ind_ft_sampleplot{j}))
            y_ft = [y_ft; 1; 1];
        else
            y_ft = [y_ft; 0; 0];
        end

        if ~isempty(intersect(k,ind_hyst_sampleplot{j}))
            y_hyst = [y_hyst; 1; 1];
        else
            y_hyst = [y_hyst; 0; 0];
        end
    end

    %make area plot
    subplot(2,N_fQ_sampleplot,j+N_fQ_sampleplot); hold on;
    h=area(t_area,[y_it,y_hyst,y_ft]);
    if(min(ind_ft_sampleplot{j}>=10))
        text(2,2.75,paneltext{2,j},'FontSize',PlotFont)
    else
        text(15,2.75,paneltext{2,j},'FontSize',PlotFont)
    end
    h(1).FaceColor = 'r';
    h(2).FaceColor = 'm';
    h(3).FaceColor = 'b';
    xlabel('time (s)');
    set(gca,'yticklabel','','ytick',[],'YMinorTick','On','Box','On','FontSize',PlotFont);
    if j==1
        ylabel('transport expectation');
        h_legend = legend([h(3) h(2) h(1)],'fluid thres.','hysteresis','impact thres.','Location','NorthEast');
        set(h_legend,'FontSize',PlotFont-4)
    end

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 7]);
    print([folder_Plots,'fQpred_sampleplot.png'],'-dpng');
end