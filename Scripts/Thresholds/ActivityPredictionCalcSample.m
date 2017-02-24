%% SCRIPT TO COMPARE ACTIVITY PREDICTION TO OBSERVATION

%% initialize
clearvars;

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,2); %sampling interval for analysis
deltat_analysis_s = seconds(deltat_analysis); %sampling interval for analysis in seconds

%% info for sample plot
Site_sampleplot = 'Oceano';
StartTime_sampleplot = [datetime(2015,5,16,15,4,0);...
    datetime(2015,5,18,16,42,0);...
    datetime(2015,5,27,16,38,0)];
EndTime_sampleplot = [datetime(2015,5,16,15,5,0);...
    datetime(2015,5,18,16,43,0);...
    datetime(2015,5,27,16,39,0)];

%% plotting info
PlotFont = 12;
Inactive_LineWidth = 1;
Active_LineWidth = 2;

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

%%%%%%%%%%%%%%%%%%%%%
% GET INFO FOR PLOT %
%%%%%%%%%%%%%%%%%%%%%

% information about measurement and sampling interval indices for analysis
ind_Deltat = find(Deltat_all == Deltat_analysis); %index for measurement interval
ind_deltat = find(deltat_all == deltat_analysis); %index for sampling interval
ind_Site = find(strcmp(Sites,Site_sampleplot)); %index for site

%initialize lists for sample plot
N_sampleplot = length(StartTime_sampleplot); %get number of intervals for plotting
ind_sampleplot = zeros(N_sampleplot,1); %initialize indices of subwindows for plotting
fQ_sampleplot = zeros(N_sampleplot,1); %initialize fQ for sampleplot
t_sampleplot = cell(N_sampleplot,1); %initialize t timeseries for sampleplot
u_sampleplot = cell(N_sampleplot,1); %initialize u timeseries for sampleplot
u_min_sampleplot = zeros(N_sampleplot,1); %initialize mininum u for each sampleplot
u_max_sampleplot = zeros(N_sampleplot,1); %initialize maximum u for each sampleplot
n_sampleplot = cell(N_sampleplot,1); %initialize n timeseries for sampleplot
zU_sampleplot = zeros(N_sampleplot,1); %initialize zU for each sampleplot
uft_sampleplot = zeros(N_sampleplot,1); %initialize uit's for sampleplot
uit_sampleplot = zeros(N_sampleplot,1); %initialize uit's for sampleplot
ind_ft_sampleplot = cell(N_sampleplot,1); %initialize indices of activity for fluid threshold hypothesis
ind_it_sampleplot = cell(N_sampleplot,1); %initialize indices of activity for impact threshold hypothesis
ind_hyst_sampleplot = cell(N_sampleplot,1); %initialize indices of activity for hysteresis hypothesis
ind_n_sampleplot = cell(N_sampleplot,1); %initialize indices of actual activity

%get threshold values for site
ustit = sqrt(tauit_all(ind_Site)/rho_Site(ind_Site));
ustft = sqrt(tauft_all(ind_Site)/rho_Site(ind_Site));

%get values for area plot
t_area = cell(N_Sites,1);
y_area_it = cell(N_Sites,1);
y_area_ft = cell(N_Sites,1);
y_area_hyst = cell(N_Sites,1);

%get values for each sample plot interval
for j = 1:N_sampleplot
    ind_sampleplot(j) = find(StartTime_subwindow{ind_Site}{ind_Deltat,ind_deltat}==StartTime_sampleplot(j));
    fQ_sampleplot(j) = fQ_subwindow_all{ind_Site}{ind_Deltat,ind_deltat}(ind_sampleplot(j));
    t_datetime = t_subwindow{ind_Site}{ind_Deltat,ind_deltat}{ind_sampleplot(j)}; %get times for sample plot (in datetime format);
    t_sampleplot{j} = seconds(t_datetime-t_datetime(1)+deltat_analysis/2); %get times for sample plot (in seconds)
    u_sampleplot{j} = u_subwindow{ind_Site}{ind_Deltat,ind_deltat}{ind_sampleplot(j)}; %get wind speeds for sample plot
    u_min_sampleplot(j) = min(u_sampleplot{j}); %get mininum wind speed
    u_max_sampleplot(j) = max(u_sampleplot{j}); %get maximum wind speed
    n_sampleplot{j} = ntotal_subwindow{ind_Site}{ind_Deltat,ind_deltat}{ind_sampleplot(j)}; %initialize n timeseries for sampleplot
    zU_sampleplot = zU_subwindow{ind_Site}{ind_Deltat,ind_deltat}(ind_sampleplot(j));
    uit_sampleplot(j) = (ustit/kappa).*log(zU_sampleplot./z0_Site(ind_Site));
    uft_sampleplot(j) = (ustft/kappa).*log(zU_sampleplot./z0_Site(ind_Site));
    ind_ft_sampleplot{j} = ind_fplus_all{ind_Site}{ind_sampleplot(j)}; %get indices of activity for fluid threshold hypothesis
    ind_it_sampleplot{j} = union(ind_fplus_all{ind_Site}{ind_sampleplot(j)},ind_fint_all{ind_Site}{ind_sampleplot(j)}); %get indices of activity for impact threshold hypothesis
    ind_hyst_sampleplot{j} = union(ind_fplus_all{ind_Site}{ind_sampleplot(j)},ind_fint_down_all{ind_Site}{ind_sampleplot(j)}); %get indices of activity for hysteresis hypothesis
    ind_n_sampleplot{j} = find(n_sampleplot{j}>0); %get indices of actual activity
    
    %get number of times
    N_t = length(t_sampleplot{j});
    
    %get interpolated start times for impact threshold
    ind_start_it = ind_it_sampleplot{j}(find(diff(ind_it_sampleplot{j})>1)+1);
    if min(ind_it_sampleplot{j})>1
        ind_start_it = [min(ind_it_sampleplot{j}); ind_start_it];
    end
    t_start_it = t_sampleplot{j}(ind_start_it)-deltat_analysis_s*...
        (u_sampleplot{j}(ind_start_it)-uit_sampleplot(j))./...
        (u_sampleplot{j}(ind_start_it)-u_sampleplot{j}(ind_start_it-1));
    
    %get interpolated end times for impact threshold
    ind_end_it = ind_it_sampleplot{j}(find(diff(ind_it_sampleplot{j})>1));
    if max(ind_it_sampleplot{j})<N_t
        ind_end_it = [ind_end_it; max(ind_it_sampleplot{j})];
    end
    t_end_it = t_sampleplot{j}(ind_end_it)+deltat_analysis_s*...
        (uit_sampleplot(j)-u_sampleplot{j}(ind_end_it))./...
        (u_sampleplot{j}(ind_end_it+1)-u_sampleplot{j}(ind_end_it));

    %get interpolated start times for fluid threshold
    ind_start_ft = ind_ft_sampleplot{j}(find(diff(ind_ft_sampleplot{j})>1)+1);
    if min(ind_ft_sampleplot{j})>1
        ind_start_ft = [min(ind_ft_sampleplot{j}); ind_start_ft];
    end
    t_start_ft = t_sampleplot{j}(ind_start_ft)-deltat_analysis_s*...
        (u_sampleplot{j}(ind_start_ft)-uft_sampleplot(j))./...
        (u_sampleplot{j}(ind_start_ft)-u_sampleplot{j}(ind_start_ft-1));
    
    %get interpolated end times for fluid threshold
    ind_end_ft = ind_ft_sampleplot{j}(find(diff(ind_ft_sampleplot{j})>1));
    if max(ind_ft_sampleplot{j})<N_t
        ind_end_ft = [ind_end_ft; max(ind_ft_sampleplot{j})];
    end
    t_end_ft = t_sampleplot{j}(ind_end_ft)+deltat_analysis_s*...
        (uft_sampleplot(j)-u_sampleplot{j}(ind_end_ft))./...
        (u_sampleplot{j}(ind_end_ft+1)-u_sampleplot{j}(ind_end_ft));
    
    %get interpolated start times for hysteresis
    ind_start_hyst = ind_hyst_sampleplot{j}(find(diff(ind_hyst_sampleplot{j})>1)+1);
    if min(ind_hyst_sampleplot{j})>1
        ind_start_hyst = [min(ind_hyst_sampleplot{j}); ind_start_hyst];
    end
    t_start_hyst = t_sampleplot{j}(ind_start_hyst)-deltat_analysis_s*...
        (u_sampleplot{j}(ind_start_hyst)-uft_sampleplot(j))./...
        (u_sampleplot{j}(ind_start_hyst)-u_sampleplot{j}(ind_start_hyst-1));
    
    %get interpolated end times for hysteresis
    ind_end_hyst = ind_hyst_sampleplot{j}(find(diff(ind_hyst_sampleplot{j})>1));
    if max(ind_hyst_sampleplot{j})<N_t
        ind_end_hyst = [ind_end_hyst; max(ind_hyst_sampleplot{j})];
    end
    t_end_hyst = t_sampleplot{j}(ind_end_hyst)+deltat_analysis_s*...
        (uit_sampleplot(j)-u_sampleplot{j}(ind_end_hyst))./...
        (u_sampleplot{j}(ind_end_hyst+1)-u_sampleplot{j}(ind_end_hyst));
    
    %initialize t_area with first and last times
    t_area{j} = t_sampleplot{j}([1; end]);
    
    %get first and last values for area plot - it
    if min(ind_it_sampleplot{j})==1
        if max(ind_it_sampleplot{j})==N_t
            y_area_it{j} = [1; 1];
        else
            y_area_it{j} = [1; 0];
        end
    else
        if max(ind_it_sampleplot{j})==N_t
            y_area_it{j} = [0; 1];
        else
            y_area_it{j} = [0; 0];
        end
    end

    %get first and last values for area plot - ft
    if min(ind_ft_sampleplot{j})==1
        if max(ind_ft_sampleplot{j})==N_t
            y_area_ft{j} = [1; 1];
        else
            y_area_ft{j} = [1; 0];
        end
    else
        if max(ind_ft_sampleplot{j})==N_t
            y_area_ft{j} = [0; 1];
        else
            y_area_ft{j} = [0; 0];
        end
    end

    %get first and last values for area plot - hyst
    if min(ind_hyst_sampleplot{j})==1
        if max(ind_hyst_sampleplot{j})==N_t
            y_area_hyst{j} = [1; 1];
        else
            y_area_hyst{j} = [1; 0];
        end
    else
        if max(ind_hyst_sampleplot{j})==N_t
            y_area_hyst{j} = [0; 1];
        else
            y_area_hyst{j} = [0; 0];
        end
    end
    
    %get remaining values - it start and end times
    for k = 1:length(t_start_it)
        t_area{j} = [t_area{j}; t_start_it(k); t_start_it(k)+0.001];
        y_area_it{j} = [y_area_it{j}; 0; 1];
        y_area_ft{j} = [y_area_ft{j}; 0; 0];
        y_area_hyst{j} = [y_area_hyst{j}; 0; 0];
    end
    for k = 1:length(t_end_it)
        t_area{j} = [t_area{j}; t_end_it(k); t_end_it(k)+0.001];
        y_area_it{j} = [y_area_it{j}; 1; 0];
        y_area_ft{j} = [y_area_ft{j}; 0; 0];
        if ~isempty(find(t_end_hyst==t_end_it(k)))
            y_area_hyst{j} = [y_area_hyst{j}; 1; 0];
        else
            y_area_hyst{j} = [y_area_hyst{j}; 0; 0];
        end
    end 

    %get remaining values - ft start and end times
    for k = 1:length(t_start_ft)
        t_area{j} = [t_area{j}; t_start_ft(k); t_start_ft(k)+0.001];
        y_area_it{j} = [y_area_it{j}; 1; 1];
        y_area_ft{j} = [y_area_ft{j}; 0; 1];
        if ~isempty(find(t_start_hyst==t_start_ft(k)))
            y_area_hyst{j} = [y_area_hyst{j}; 0; 1];
        else
            y_area_hyst{j} = [y_area_hyst{j}; 1; 1];
        end
    end
    for k = 1:length(t_end_ft)
        t_area{j} = [t_area{j}; t_end_ft(k); t_end_ft(k)+0.001];
        y_area_it{j} = [y_area_it{j}; 1; 1];
        y_area_ft{j} = [y_area_ft{j}; 1; 0];
        y_area_hyst{j} = [y_area_hyst{j}; 1; 1];
    end
    
    %sort out values
    [t_area{j}, ind_sort] = sort(t_area{j});
    y_area_it{j} = y_area_it{j}(ind_sort);
    y_area_hyst{j} = y_area_hyst{j}(ind_sort);    
    y_area_ft{j} = y_area_ft{j}(ind_sort);
    
    %display values for paper
    StartTime = StartTime_sampleplot(j)
    fQpred_ft = round(length(ind_ft_sampleplot{j})/length(t_sampleplot{j}),2)
    fQpred_it = round(length(ind_it_sampleplot{j})/length(t_sampleplot{j}),2)
    fQpred_hyst = round(length(ind_hyst_sampleplot{j})/length(t_sampleplot{j}),2)
end

%%%%%%%%
% PLOT %
%%%%%%%%
fig = figure(1); clf;
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]); 
paneltext = {'(a)','(b)','(c)'}; %panel labels

xlims_sampleplot = [0 seconds(Deltat_analysis)];
ylims_sampleplot = [floor(min(u_min_sampleplot))-1, ceil(max(u_max_sampleplot))];

for j = 1:N_sampleplot
        
    %initialize subplot
    subplot('Position',[-0.25+0.31*j 0.11 0.28 0.84]); hold on;
    
    %% left axis -- timeseries
    yyaxis left
    
    %plot fluid and impact threshold
    h2=plot(xlims_sampleplot,uft_sampleplot(j)*[1 1],'b-.','LineWidth',2);
    h3=plot(xlims_sampleplot,uit_sampleplot(j)*[1 1],'g--','LineWidth',2);
    
    %plot timeseries
    h1=plot(t_sampleplot{j},u_sampleplot{j},'k','LineWidth',1);

    %format plot
    text(2,ylims_sampleplot(2)-0.25,paneltext{1,j},'FontSize',PlotFont)
    xlim(xlims_sampleplot);
    ylim(ylims_sampleplot);
    xlabel('time (s)');
    if j==1
        ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
    end
    title(['$$f_{Q} = ',num2str(fQ_sampleplot(j),2),'$$'],'Interpreter','Latex');
    set(gca,'XMinorTick','On','YMinorTick','On','YTick',ylims_sampleplot(1)+1:ylims_sampleplot(2),'Box','On','FontSize',PlotFont);
    
    %% right axis - transport occurrence
    yyaxis right
    
    %plot transport expectation
    h4=area(t_area{j},[y_area_it{j},y_area_hyst{j},y_area_ft{j}]);
    h4(1).FaceColor = 'g';
    h4(2).FaceColor = 'm';
    h4(3).FaceColor = 'b';
    
    %format plot
    xlim(xlims_sampleplot);
    ylim([0 12]);
    if j==1
       text(27.5,11.8,'16 May 2015, 15:04-15:05','FontSize',PlotFont-4)
    elseif j==2
       text(27.5,11.8,'18 May 2015, 16:42-16:43','FontSize',PlotFont-4)
       ylabel('predicted transport','Interpreter','Latex','Rotation',270,'Position',[135, 1.5, 0]);
    elseif j==3
       text(27.5,11.8,'27 May 2015, 16:38-16:39','FontSize',PlotFont-4)
       ylabel('occurrence','Interpreter','Latex','Rotation',270,'Position',[65, 1.5, 0]);
    end
    set(gca,'yticklabel','','ytick',[],'YMinorTick','On','Box','On','FontSize',PlotFont);
    
    if j==2
        h_legend = legend([h1 h2 h3 h4(3) h4(2) h4(1)],...
            'obs. wind','fluid thres., u_{ft}','impact thres., u_{it}',...
            'fluid thres. pred.','dual thres. pred.','impact thres. pred.',...
            'Location','NorthEast');
        set(h_legend,'FontSize',PlotFont)
    end
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[9 5],'PaperPosition',[0 0 9 5],'PaperPositionMode','Manual');
print([folder_Plots,'fQpred_sampleplot.png'],'-dpng');
print([folder_Plots,'fQpred_sampleplot.tif'],'-dtiff','-r600');

% %%%%%%%%
% % PLOT %
% %%%%%%%%
% figure(1); clf;
% paneltext = {'(a)','(b)','(c)';'(d)','(e)','(f)'}; %panel labels
% 
% xlims_sampleplot = [0 seconds(Deltat_analysis)];
% ylims_sampleplot = [floor(min(u_min_sampleplot)), ceil(max(u_max_sampleplot))];
% 
% for j = 1:N_sampleplot
%     %plot timeseries
%     subplot(2,N_sampleplot,j); hold on;
%     plot(t_sampleplot{j},u_sampleplot{j},'k','LineWidth',1);
%     plot(xlims_sampleplot,uft_sampleplot(j)*[1 1],'b-.','LineWidth',2);
%     plot(xlims_sampleplot,uit_sampleplot(j)*[1 1],'r--','LineWidth',2);
%     plot(t_sampleplot{j},u_sampleplot{j},'k','LineWidth',1);
%     text(2,ylims_sampleplot(2)-0.5,paneltext{1,j},'FontSize',PlotFont)
%     xlim(xlims_sampleplot);
%     ylim(ylims_sampleplot);
%     xlabel('time (s)');
%     if j==1
%         ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
%     elseif j==2
%         h_legend = legend('wind speed','fluid thres., u_{ft}','impact thres., u_{it}');
%         set(h_legend,'FontSize',PlotFont-4,'Location','NorthEast');
%     end
%     title(['$$f_{Q} = ',num2str(fQ_sampleplot(j),1),'$$'],'Interpreter','Latex');
%     set(gca,'XMinorTick','On','YTick',ylims_sampleplot(1):ylims_sampleplot(2),'YMinorTick','On','Box','On','FontSize',PlotFont);
% 
%     %get info for area plot
%     t_area = [];
%     y_it = [];
%     y_ft = [];
%     y_hyst = [];
%     for k = 1:length(t_sampleplot{j})
%         t_area = [t_area; t_sampleplot{j}(k)-seconds(deltat_analysis)/2; t_sampleplot{j}(k)+seconds(deltat_analysis)/2];
%         if ~isempty(intersect(k,ind_it_sampleplot{j}))
%             y_it = [y_it; 1; 1];
%         else
%             y_it = [y_it; 0; 0];
%         end
% 
%         if ~isempty(intersect(k,ind_ft_sampleplot{j}))
%             y_ft = [y_ft; 1; 1];
%         else
%             y_ft = [y_ft; 0; 0];
%         end
% 
%         if ~isempty(intersect(k,ind_hyst_sampleplot{j}))
%             y_hyst = [y_hyst; 1; 1];
%         else
%             y_hyst = [y_hyst; 0; 0];
%         end
%     end
% 
%     %make area plot
%     subplot(2,N_sampleplot,j+N_sampleplot); hold on;
%     h=area(t_area,[y_it,y_hyst,y_ft]);
%     if(min(ind_ft_sampleplot{j})>=10 || min(ind_ft_sampleplot{j})<3)
%         text(2,2.75,paneltext{2,j},'FontSize',PlotFont)
%     else
%         text(15,2.75,paneltext{2,j},'FontSize',PlotFont)
%     end
%     h(1).FaceColor = 'r';
%     h(2).FaceColor = 'm';
%     h(3).FaceColor = 'b';
%     xlabel('time (s)');
%     set(gca,'yticklabel','','ytick',[],'YMinorTick','On','Box','On','FontSize',PlotFont);
%     if j==1
%         ylabel('transport expectation');
%     elseif j==2
%         h_legend = legend([h(3) h(2) h(1)],'fluid thres.','dual thres.','impact thres.','Location','NorthEast');
%         set(h_legend,'FontSize',PlotFont-4)
%     end
%     ylim([0 3]);
%     
%     %print plot
%     set(gca, 'LooseInset', get(gca,'TightInset'));
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 7]);
%     print([folder_Plots,'fQpred_sampleplot.png'],'-dpng');
% end


% %%%%%%%%%%%%%%%%%%
% % ALTERNATE PLOT %
% %%%%%%%%%%%%%%%%%%
% 
% figure(1); clf;
% paneltext = {'(a)','(b)','(c)'}; %panel labels
% 
% xlims_sampleplot = [0 seconds(Deltat_analysis)];
% ylims_sampleplot = [floor(min(u_min_sampleplot))-1, ceil(max(u_max_sampleplot))];
% 
% for j = 1:N_sampleplot
%     
%     %get times and wind speeds for starting line from it
%     ind_start_it = ind_it_sampleplot{j}(find(diff(ind_it_sampleplot{j})>1)+1);
%     if min(ind_it_sampleplot{j})>1
%         ind_start_it = [min(ind_it_sampleplot{j}); ind_start_it];
%     end
%     t_start_it = t_sampleplot{j}(ind_start_it)-...
%         deltat_analysis_s*...
%         (u_sampleplot{j}(ind_start_it)-uit_sampleplot(j))./...
%         (u_sampleplot{j}(ind_start_it)-u_sampleplot{j}(ind_start_it-1));
%     u_start_it = uit_sampleplot(j)*ones(size(t_start_it));
%     
%     %get times and wind speeds for ending line from it
%     ind_end_it = ind_it_sampleplot{j}(find(diff(ind_it_sampleplot{j})>1));
%     if max(ind_it_sampleplot{j})<length(t_sampleplot{j})
%         ind_end_it = [ind_end_it; max(ind_it_sampleplot{j})];
%     end
%     t_end_it = t_sampleplot{j}(ind_end_it)+...
%         deltat_analysis_s*...
%         (uit_sampleplot(j)-u_sampleplot{j}(ind_end_it))./...
%         (u_sampleplot{j}(ind_end_it+1)-u_sampleplot{j}(ind_end_it));
%     u_end_it = uit_sampleplot(j)*ones(size(t_end_it));
% 
%     %get times and wind speeds for starting line from ft
%     ind_start_ft = ind_ft_sampleplot{j}(find(diff(ind_ft_sampleplot{j})>1)+1);
%     if min(ind_ft_sampleplot{j})>1
%         ind_start_ft = [min(ind_ft_sampleplot{j}); ind_start_ft];
%     end
%     t_start_ft = t_sampleplot{j}(ind_start_ft)-...
%         deltat_analysis_s*...
%         (u_sampleplot{j}(ind_start_ft)-uft_sampleplot(j))./...
%         (u_sampleplot{j}(ind_start_ft)-u_sampleplot{j}(ind_start_ft-1));
%     u_start_ft = uft_sampleplot(j)*ones(size(t_start_ft));
%     
%     %get times and wind speeds for ending line from ft
%     ind_end_ft = ind_ft_sampleplot{j}(find(diff(ind_ft_sampleplot{j})>1));
%     if max(ind_ft_sampleplot{j})<length(t_sampleplot{j})
%         ind_end_ft = [ind_end_ft; max(ind_ft_sampleplot{j})];
%     end
%     t_end_ft = t_sampleplot{j}(ind_end_ft)+...
%         deltat_analysis_s*...
%         (uft_sampleplot(j)-u_sampleplot{j}(ind_end_ft))./...
%         (u_sampleplot{j}(ind_end_ft+1)-u_sampleplot{j}(ind_end_ft));
%     u_end_ft = uft_sampleplot(j)*ones(size(t_end_ft));
%     
%     %get all times for it line
%     t_it = [t_sampleplot{j}; t_start_it; t_end_it]; %all times
%     [t_it,ind_sort_it] = sort(t_it); %all times, sorted
% 
%     %get all times for ft line
%     t_ft = [t_sampleplot{j}; t_start_ft; t_end_ft]; %all times
%     [t_ft,ind_sort_ft] = sort(t_ft); %all times, sorted
% 
%     %get all times for hyst line
%     t_hyst = [t_sampleplot{j}; t_start_ft; t_end_it]; %all times
%     [t_hyst,ind_sort_hyst] = sort(t_hyst); %all times, sorted
%     
%     %get wind speeds for it line
%     u_it = zeros(size(u_sampleplot{j}))*NaN; %initialize wind speeds as NaN
%     u_it(ind_it_sampleplot{j})=u_sampleplot{j}(ind_it_sampleplot{j}); %keep only those that are not NaN
%     u_it = [u_it; u_start_it; u_end_it]; %add in addition wind speeds for starting and ending
%     u_it = u_it(ind_sort_it); %sort these according to sorting of t
% 
%     %get wind speeds for ft line
%     u_ft = zeros(size(u_sampleplot{j}))*NaN; %initialize wind speeds as NaN
%     u_ft(ind_ft_sampleplot{j})=u_sampleplot{j}(ind_ft_sampleplot{j}); %keep only those that are not NaN
%     u_ft = [u_ft; u_start_ft; u_end_ft]; %add in addition wind speeds for starting and ending
%     u_ft = u_ft(ind_sort_ft); %sort these according to sorting of t
% 
%     %get wind speeds for hyst line
%     u_hyst = zeros(size(u_sampleplot{j}))*NaN; %initialize wind speeds as NaN
%     u_hyst(ind_hyst_sampleplot{j})=u_sampleplot{j}(ind_hyst_sampleplot{j}); %keep only those that are not NaN
%     u_hyst = [u_hyst; u_start_ft; u_end_it]; %add in addition wind speeds for starting and ending
%     u_hyst = u_hyst(ind_sort_hyst); %sort these according to sorting of t
%     
%     %get info for area plot of actual transport occurence
%     t_area = [];
%     y_area = [];
%     for k = 1:length(t_sampleplot{j})
%         t_area = [t_area; t_sampleplot{j}(k)-seconds(deltat_analysis)/2; t_sampleplot{j}(k)+seconds(deltat_analysis)/2];
%         if ~isempty(intersect(k,ind_n_sampleplot{j}))
%             y_area = [y_area; ylims_sampleplot(1)+1; ylims_sampleplot(1)+1];
%         else
%             y_area = [y_area; 0; 0];
%         end
%     end
%         
%     %initialize subplot
% %    subplot(1,N_sampleplot,j); hold on;
%     subplot('Position',[-0.26+0.32*j 0.14 0.28 0.82]); hold on;
%        
%     %plot fluid and impact thresholds
%     plot(xlims_sampleplot,uit_sampleplot(j)*[1 1],'r--','LineWidth',Active_LineWidth);
%     plot(xlims_sampleplot,uft_sampleplot(j)*[1 1],'b-.','LineWidth',Active_LineWidth);
% 
%     %plot wind timeseries
%     plot(t_sampleplot{j},u_sampleplot{j},'k','LineWidth',Inactive_LineWidth);
%     
%     %make timeseries bold for predictions
%     plot(t_it,u_it,'r','LineWidth', Active_LineWidth); %impact threshold timeseries
%     plot(t_hyst,u_hyst,'m','LineWidth', Active_LineWidth); %hysteresis timeseries
%     plot(t_ft,u_ft,'b','LineWidth', Active_LineWidth); %fluid threshold timeseries
%    
%     %plot times of observed transport
%     area(t_area,y_area,'FaceColor','k');
%     
%     %format plot
%     if j==1
%         ylabel('wind speed, $$u$$ (m/s)','Interpreter','Latex');
%     end
%     if j==2
%         h_legend = legend('impact thres','fluid thres','no trans pred','it only pred','it/hyst pred','it/hyst/ft pred','obs trans');
%         set(h_legend,'FontSize',PlotFont-4,'Location','North');
%     end
%     text(2,ylims_sampleplot(2)-0.5,paneltext{j},'FontSize',PlotFont)
%     xlim(xlims_sampleplot);
%     ylim(ylims_sampleplot);
%     xlabel('time (s)');
%     %title(['$$f_{Q} = ',num2str(fQ_sampleplot(j),1),'$$'],'Interpreter','Latex');
%     set(gca,'XMinorTick','On','YTick',ylims_sampleplot(1):ylims_sampleplot(2),'YMinorTick','On','Box','On','FontSize',PlotFont);
%   
%     %print plot
%     set(gca, 'LooseInset', get(gca,'TightInset'));
%     set(gcf,'PaperUnits','inches','PaperSize',[9 3.5],'PaperPosition',[0 0 9 3.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'fQpred_sampleplot_line.png'],'-dpng');
% end