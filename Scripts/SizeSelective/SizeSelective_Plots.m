%% initialize
clearvars;
close all;

%% select plot type
plot_type = 'paper';
%plot_type = 'presentation';

%%
%%%%%%%%%%%%%%%%
% DATA SOURCES %
%%%%%%%%%%%%%%%%

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/SizeSelective/'; %folder for saving analysis data
AnalysisData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis');
CorrectionData_Path = strcat(folder_AnalysisData,'PSDcorrectionFactor.dat');
folder_Plots = '../../PlotOutput/SizeSelective/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(AnalysisData_Path); %load analysis data
addpath(folder_Functions); %point MATLAB to location of functions

%% load correction data
fileID = fopen(CorrectionData_Path,'r');
formatSpec = '%f%f%f%[^\n\r]';
delimiter = '\t';
startRow = 2;
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
d_mid_correction = dataArray{1}/1000; %mid point diameter (mm)
fC_4cm = dataArray{2}; %correction factor for measurements above 4 cm
fC_7cm = dataArray{3}; %correction factor for measurements above 7 cm

%%
%%%%%%%%%%%%%%%%%
% PLOTTING INFO %
%%%%%%%%%%%%%%%%%

%% plotting information
PlotFont = 12; %font for labels
LineWidth_Plot = 1; %width of lines
Marker_Cluster = {'s','d','o','p','h','^','v','>'}; %markers for Clusters
Color_Cluster = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]}; %colors for sites / clusters
Label_Cluster = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'}; %markers for Clusters
Marker_bin = {'o','s','d','^','v','>','<','p','h','x','*','+'}; %markers for bins

%% colors for bins by wind strength
Color_bin = cell(N_bins,1);
for i = 1:N_bins
    Color_bin{i} = [0, (i-1)/(N_bins-1), 1-(i-1)/(N_bins-1)];
end
Color_ust_bin = cell(N_ust_bins,1);
for i = 1:N_ust_bins
    Color_ust_bin{i} = [(i-1)/(N_ust_bins-1), 0, 1-(i-1)/(N_ust_bins-1)];
end
Color_ustnorm_bin = cell(N_ustnorm_bins,1);
for i = 1:N_ustnorm_bins
    Color_ustnorm_bin{i} = [(i-1)/(N_ustnorm_bins-1), 0, 1-(i-1)/(N_ustnorm_bins-1)];
end
Color_taunorm_bin = cell(N_taunorm_bins,1);
for i = 1:N_taunorm_bins
    Color_taunorm_bin{i} = [(i-1)/(N_taunorm_bins-1), 0, 1-(i-1)/(N_taunorm_bins-1)];
end

%
%%%%%%%%%%%%%%%%%%
% MAIN TEXT PLOTS %
%%%%%%%%%%%%%%%%%%

%% Plot taunorm-conditioned airborne versus surface size distributions - surface only
figure(1); clf;

%initialize subplots
h_subplot = gobjects(N_Cluster,1);

for i = 1:N_Cluster
    
    %initialize subplot
    if strcmp(plot_type,'presentation')
        %initialize subplot - landscape
        if N_Cluster == 6 %defined subplot sizes for six clusters - landscape
            if i == 1
                h_subplot(1) = subplot('position',[0.1 0.54 0.25 0.40]); hold on;
            elseif i == 2
                h_subplot(2) = subplot('position',[0.41 0.54 0.25 0.40]); hold on;
            elseif i == 3
                h_subplot(3) = subplot('position',[0.72 0.54 0.25 0.40]); hold on;
            elseif i == 4
                h_subplot(4) = subplot('position',[0.1 0.06 0.25 0.40]); hold on;
            elseif i == 5
                h_subplot(5) = subplot('position',[0.41 0.06 0.25 0.40]); hold on;
            else
                h_subplot(6) = subplot('position',[0.72 0.06 0.25 0.40]); hold on;
            end
        else %otherwise, automated subplot sizes
            h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
        end
    elseif strcmp(plot_type,'paper')
        %initialize subplot - portrait
        if N_Cluster == 6 %defined subplot sizes for six clusters
            if i == 1
                h_subplot(1) = subplot('position',[0.12 0.7 0.38 0.24]); hold on;
            elseif i == 2
                h_subplot(2) = subplot('position',[0.58 0.7 0.38 0.24]); hold on;
            elseif i == 3
                h_subplot(3) = subplot('position',[0.12 0.38 0.38 0.24]); hold on;
            elseif i == 4
                h_subplot(4) = subplot('position',[0.58 0.38 0.38 0.24]); hold on;
            elseif i == 5
                h_subplot(5) = subplot('position',[0.12 0.06 0.38 0.24]); hold on;
            else
                h_subplot(6) = subplot('position',[0.58 0.06 0.38 0.24]); hold on;
            end
        else %otherwise, automated subplot sizes
            h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
        end
    end
        
    %plot surface distribution
    plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
        
    %format plot
    set(gca,'XScale','log','YScale','log','YMinorTick','On');
    xlim([0.06, 2]);
    set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
    set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'});
    ylim([1e-4 1e1]);

    %label plot
    htitle = title(ClusterNames{i});
    set(htitle,'Position',[0.35,3.5]); %set title below edge of box to accommodate second axis

    %labels 
    if strcmp(plot_type,'presentation')
        %xlabel - landscape
        if i>=2*round(N_Cluster/2)-2
            xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
        end
        %ylabel - landscape
        if mod(i,3) == 1
            ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
        end
    elseif strcmp(plot_type,'paper')
        %xlabel - portrait
        if i>=2*round(N_Cluster/2)-1
            xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
        end
        %ylabel - portrait
        if mod(i,2) == 1
            ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
        end
        %label for subplot id
        text(0.07, 6, Label_Cluster{i},'FontSize',14);
    end
    
    %add second axis
    ax1 = gca; %get handle for first axis
    set(ax1,'XColor','k','YColor','k');
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k',...
        'YColor','k');
    cla;
    xmin = 0.06/d50_bar_surface_Cluster(i);
    xmax = 2/d50_bar_surface_Cluster(i);
    line([xmin xmax],[0 0],'Color','k','Parent',ax2);
    set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2]);
    xmaxtick = floor(xmax*10)/10;
    set(ax2,'XTick',0.2:0.1:xmaxtick);
    XTickLabels = cell(8+xmaxtick*10,1);
    XTickLabels{1} = '0.2';
    XTickLabels{3} = '0.4';
    XTickLabels{6} = '0.7';
    for j = 1:floor(xmaxtick)
        XTickLabels{j*10-1} = int2str(j);
    end
    set(ax2,'XTickLabel',XTickLabels);
    set(ax2,'YTick',[]);
    
    %secondary xlabel
    if strcmp(plot_type,'presentation')
        %xlabel - landscape
        if i<=3
            xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
        end
    elseif strcmp(plot_type,'paper')
        %xlabel - portrait
        if i<=2
            xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
        end
    end
end

if strcmp(plot_type,'presentation')
    %print plot - landscape
    set(gcf,'PaperUnits','inches','PaperSize',[8 7],'PaperPosition',[0 0 8 7],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surfaceonly_landscape.png'],'-dpng');
elseif strcmp(plot_type,'paper')
    %print plot - portrait
    set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surfaceonly_portrait.png'],'-dpng');
end


%% Plot taunorm-conditioned airborne versus surface size distributions - surface and 1<tau<1.5 only
figure(2); clf;

%initialize subplots
h_subplot = gobjects(N_Cluster,1);

for i = 1:N_Cluster
    
    %initialize subplot
    if strcmp(plot_type,'presentation')
        %initialize subplot - landscape
        if N_Cluster == 6 %defined subplot sizes for six clusters - landscape
            if i == 1
                h_subplot(1) = subplot('position',[0.1 0.54 0.25 0.40]); hold on;
            elseif i == 2
                h_subplot(2) = subplot('position',[0.41 0.54 0.25 0.40]); hold on;
            elseif i == 3
                h_subplot(3) = subplot('position',[0.72 0.54 0.25 0.40]); hold on;
            elseif i == 4
                h_subplot(4) = subplot('position',[0.1 0.06 0.25 0.40]); hold on;
            elseif i == 5
                h_subplot(5) = subplot('position',[0.41 0.06 0.25 0.40]); hold on;
            else
                h_subplot(6) = subplot('position',[0.72 0.06 0.25 0.40]); hold on;
            end
        else %otherwise, automated subplot sizes
            h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
        end
    elseif strcmp(plot_type,'paper')
        %initialize subplot - portrait
        if N_Cluster == 6 %defined subplot sizes for six clusters
            if i == 1
                h_subplot(1) = subplot('position',[0.12 0.7 0.38 0.24]); hold on;
            elseif i == 2
                h_subplot(2) = subplot('position',[0.58 0.7 0.38 0.24]); hold on;
            elseif i == 3
                h_subplot(3) = subplot('position',[0.12 0.38 0.38 0.24]); hold on;
            elseif i == 4
                h_subplot(4) = subplot('position',[0.58 0.38 0.38 0.24]); hold on;
            elseif i == 5
                h_subplot(5) = subplot('position',[0.12 0.06 0.38 0.24]); hold on;
            else
                h_subplot(6) = subplot('position',[0.58 0.06 0.38 0.24]); hold on;
            end
        else %otherwise, automated subplot sizes
            h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
        end
    end
        
    %plot surface distribution
    plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
    
    %plot airborne distribution - tau<2.2
    for k = 1:3
        plot(d_taunormbar_airborne_mid_drange_Cluster{i},...
            dVdlogd_taunormbar_airborne_drange_Cluster{i}(k,:),...
            ['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_taunorm_bin{k})
    end
    
    %format plot
    set(gca,'XScale','log','YScale','log','YMinorTick','On');
    xlim([0.06, 2]);
    set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
    set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'});
    ylim([1e-4 1e1]);

    %label plot
    htitle = title(ClusterNames{i});
    set(htitle,'Position',[0.35,3.5]); %set title below edge of box to accommodate second axis

    %labels 
    if strcmp(plot_type,'presentation')
        %xlabel - landscape
        if i>=2*round(N_Cluster/2)-2
            xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
        end
        %ylabel - landscape
        if mod(i,3) == 1
            ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
        end
    elseif strcmp(plot_type,'paper')
        %xlabel - portrait
        if i>=2*round(N_Cluster/2)-1
            xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
        end
        %ylabel - portrait
        if mod(i,2) == 1
            ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
        end
        %label for subplot id
        text(0.07, 6, Label_Cluster{i},'FontSize',14);
    end
    
    %create legend
    legend_items = cell(4,1);
    legend_items{1} = 'Surface';
    for j = 1:3
        if j == 1
            legend_items{j+1} = ['\tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];       
        else
            legend_items{j+1} = [num2str(taunorm_min_bins(j),'%10.1f'),' < \tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];
        end
    end
    if i == 1
        h_legend = legend(legend_items,'Location','SouthWest');
        set(h_legend,'FontSize',8);
    end
    
    %add second axis
    ax1 = gca; %get handle for first axis
    set(ax1,'XColor','k','YColor','k');
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k',...
        'YColor','k');
    cla;
    xmin = 0.06/d50_bar_surface_Cluster(i);
    xmax = 2/d50_bar_surface_Cluster(i);
    line([xmin xmax],[0 0],'Color','k','Parent',ax2);
    set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2]);
    xmaxtick = floor(xmax*10)/10;
    set(ax2,'XTick',0.2:0.1:xmaxtick);
    XTickLabels = cell(8+xmaxtick*10,1);
    XTickLabels{1} = '0.2';
    XTickLabels{3} = '0.4';
    XTickLabels{6} = '0.7';
    for j = 1:floor(xmaxtick)
        XTickLabels{j*10-1} = int2str(j);
    end
    set(ax2,'XTickLabel',XTickLabels);
    set(ax2,'YTick',[]);
    
    %secondary xlabel
    if strcmp(plot_type,'presentation')
        %xlabel - landscape
        if i<=3
            xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
        end
    elseif strcmp(plot_type,'paper')
        %xlabel - portrait
        if i<=2
            xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
        end
    end
end

if strcmp(plot_type,'presentation')
    %print plot - landscape
    set(gcf,'PaperUnits','inches','PaperSize',[8 7],'PaperPosition',[0 0 8 7],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surface_lowtau_only_landscape.png'],'-dpng');
elseif strcmp(plot_type,'paper')
    %print plot - portrait
    set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surface_lowtau_only_portrait.png'],'-dpng');
end


%% Plot taunorm-conditioned airborne versus surface size distributions - landscape
figure(3); clf;

%initialize subplots
h_subplot = gobjects(N_Cluster,1);

for i = 1:N_Cluster
    
    %initialize subplot
    if strcmp(plot_type,'presentation')
        %initialize subplot - landscape
        if N_Cluster == 6 %defined subplot sizes for six clusters - landscape
            if i == 1
                h_subplot(1) = subplot('position',[0.1 0.54 0.25 0.40]); hold on;
            elseif i == 2
                h_subplot(2) = subplot('position',[0.41 0.54 0.25 0.40]); hold on;
            elseif i == 3
                h_subplot(3) = subplot('position',[0.72 0.54 0.25 0.40]); hold on;
            elseif i == 4
                h_subplot(4) = subplot('position',[0.1 0.06 0.25 0.40]); hold on;
            elseif i == 5
                h_subplot(5) = subplot('position',[0.41 0.06 0.25 0.40]); hold on;
            else
                h_subplot(6) = subplot('position',[0.72 0.06 0.25 0.40]); hold on;
            end
        else %otherwise, automated subplot sizes
            h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
        end
    elseif strcmp(plot_type,'paper')
        %initialize subplot - portrait
        if N_Cluster == 6 %defined subplot sizes for six clusters
            if i == 1
                h_subplot(1) = subplot('position',[0.12 0.7 0.38 0.24]); hold on;
            elseif i == 2
                h_subplot(2) = subplot('position',[0.58 0.7 0.38 0.24]); hold on;
            elseif i == 3
                h_subplot(3) = subplot('position',[0.12 0.38 0.38 0.24]); hold on;
            elseif i == 4
                h_subplot(4) = subplot('position',[0.58 0.38 0.38 0.24]); hold on;
            elseif i == 5
                h_subplot(5) = subplot('position',[0.12 0.06 0.38 0.24]); hold on;
            else
                h_subplot(6) = subplot('position',[0.58 0.06 0.38 0.24]); hold on;
            end
        else %otherwise, automated subplot sizes
            h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
        end
    end
        
    %plot surface distribution
    plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
    
    %plot airborne distribution
    for k = 1:N_taunorm_bins
        plot(d_taunormbar_airborne_mid_drange_Cluster{i},...
            dVdlogd_taunormbar_airborne_drange_Cluster{i}(k,:),...
            ['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_taunorm_bin{k})
    end
    
    %format plot
    set(gca,'XScale','log','YScale','log','YMinorTick','On','FontSize',12);
    xlim([0.06, 2]);
    set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
    set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'},'FontSize',12);
    ylim([1e-4 1e1]);

    %label plot
    htitle = title(ClusterNames{i});
    set(htitle,'Position',[0.35,3.5],'FontSize',13); %set title below edge of box to accommodate second axis

    %labels 
    if strcmp(plot_type,'presentation')
        %xlabel - landscape
        if i>=2*round(N_Cluster/2)-2
            xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex','FontSize',13)
        end
        %ylabel - landscape
        if mod(i,3) == 1
            ylabel('Non-dim. volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex','FontSize',13);
        end
    elseif strcmp(plot_type,'paper')
        %xlabel - portrait
        if i>=2*round(N_Cluster/2)-1
            xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex','FontSize',13)
        end
        %ylabel - portrait
        if mod(i,2) == 1
            ylabel('Non-dim. volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex','FontSize',13);
        end
        %label for subplot id
        text(0.07, 4, Label_Cluster{i},'FontSize',20);
    end
    
    %create legend
    legend_items = cell(N_taunorm_bins+1,1);
    legend_items{1} = 'Surface';
    for j = 1:N_taunorm_bins
        if j == 1
            legend_items{j+1} = ['\tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];       
        else
            legend_items{j+1} = [num2str(taunorm_min_bins(j),'%10.1f'),' < \tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];
        end
    end
    if i == 1
        h_legend = legend(legend_items,'Location','SouthWest');
        set(h_legend,'FontSize',9);
    end
    
    %add second axis
    ax1 = gca; %get handle for first axis
    set(ax1,'XColor','k','YColor','k','FontSize',12);
    ax2 = axes('Position',get(ax1,'Position'),...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'XColor','k',...
        'YColor','k');
    cla;
    xmin = 0.06/d50_bar_surface_Cluster(i);
    xmax = 2/d50_bar_surface_Cluster(i);
    line([xmin xmax],[0 0],'Color','k','Parent',ax2);
    set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2],'FontSize',12);
    xmaxtick = floor(xmax*10)/10;
    set(ax2,'XTick',0.2:0.1:xmaxtick);
    XTickLabels = cell(8+xmaxtick*10,1);
    XTickLabels{1} = '0.2';
    XTickLabels{3} = '0.4';
    XTickLabels{6} = '0.7';
    for j = 1:floor(xmaxtick)
        XTickLabels{j*10-1} = int2str(j);
    end
    set(ax2,'XTickLabel',XTickLabels,'FontSize',12);
    set(ax2,'YTick',[]);
    
    %secondary xlabel
    if strcmp(plot_type,'presentation')
        %xlabel - landscape
        if i<=3
            xlabel('Non-dim. grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex','FontSize',13)
        end
    elseif strcmp(plot_type,'paper')
        %xlabel - portrait
        if i<=2
            xlabel('Non-dim. grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex','FontSize',13)
        end
    end
end

if strcmp(plot_type,'presentation')
    %print plot - landscape
    set(gcf,'PaperUnits','inches','PaperSize',[8 7],'PaperPosition',[0 0 8 7],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_landscape.png'],'-dpng');
elseif strcmp(plot_type,'paper')
    %print plot - portrait
    set(gcf,'PaperUnits','inches','PaperSize',[7 10],'PaperPosition',[0 0 7 10],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_portrait.png'],'-dpng');
end


%% Plot taunorm-conditioned airborne versus surface size distributions - N1 only - surface only
figure(11); clf; hold on;

i = 3; %only N1
        
%plot surface distribution
plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
        
%format plot
set(gca,'XScale','log','YScale','log','YMinorTick','On');
xlim([0.06, 2]);
set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'});
ylim([1e-4 1e1]);

%label plot
htitle = title(ClusterNames{i});
set(htitle,'Position',[0.35,3.5],'FontSize',13); %set title below edge of box to accommodate second axis

%labels 
if strcmp(plot_type,'presentation')
    %xlabel - landscape
    xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
elseif strcmp(plot_type,'paper')
    %xlabel - portrait
    xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    %ylabel - portrait
    ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
end
    
%add second axis
ax1 = gca; %get handle for first axis
set(ax1,'XColor','k','YColor','k');
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k',...
    'YColor','k');
cla;
xmin = 0.06/d50_bar_surface_Cluster(i);
xmax = 2/d50_bar_surface_Cluster(i);
line([xmin xmax],[0 0],'Color','k','Parent',ax2);
set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2]);
xmaxtick = floor(xmax*10)/10;
set(ax2,'XTick',0.2:0.1:xmaxtick);
XTickLabels = cell(8+xmaxtick*10,1);
XTickLabels{1} = '0.2';
XTickLabels{3} = '0.4';
XTickLabels{6} = '0.7';
for j = 1:floor(xmaxtick)
    XTickLabels{j*10-1} = int2str(j);
end
set(ax2,'XTickLabel',XTickLabels);
set(ax2,'YTick',[]);

%secondary xlabel
if strcmp(plot_type,'presentation')
    %xlabel - landscape
    xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
elseif strcmp(plot_type,'paper')
    %xlabel - portrait
    xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
end

if strcmp(plot_type,'presentation')
    %print plot - landscape
    set(gcf,'PaperUnits','inches','PaperSize',[7 6],'PaperPosition',[0 0 7 6],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surfaceonly_N1_landscape.png'],'-dpng');
elseif strcmp(plot_type,'paper')
    %print plot - portrait
    set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surfaceonly_N1_portrait.png'],'-dpng');
end


%% Plot taunorm-conditioned airborne versus surface size distributions - surface and 1<tau<1.5 only - N1 only
figure(12); clf; hold on;

i = 3; %only N1
        
%plot surface distribution
plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);

%plot airborne distribution - 1<tau<1.5
k = 2;
plot(d_taunormbar_airborne_mid_drange_Cluster{i},...
    dVdlogd_taunormbar_airborne_drange_Cluster{i}(k,:),...
    ['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_taunorm_bin{k})

%format plot
set(gca,'XScale','log','YScale','log','YMinorTick','On');
xlim([0.06, 2]);
set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'});
ylim([1e-4 1e1]);

%label plot
htitle = title(ClusterNames{i});
set(htitle,'Position',[0.35,3.5]); %set title below edge of box to accommodate second axis

%labels 
if strcmp(plot_type,'presentation')
    %xlabel - landscape
    xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    %ylabel - landscape
    ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
elseif strcmp(plot_type,'paper')
    %xlabel - portrait
    xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    %ylabel - portrait
    ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
end

%create legend
legend_items = cell(2,1);
legend_items{1} = 'Surface';
legend_items{2} = [num2str(taunorm_min_bins(2),'%10.1f'),' < \tau/\tau_{it} \leq ',num2str(taunorm_max_bins(2),'%10.1f')];
h_legend = legend(legend_items,'Location','SouthWest');
set(h_legend,'FontSize',10);

%add second axis
ax1 = gca; %get handle for first axis
set(ax1,'XColor','k','YColor','k');
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k',...
    'YColor','k');
cla;
xmin = 0.06/d50_bar_surface_Cluster(i);
xmax = 2/d50_bar_surface_Cluster(i);
line([xmin xmax],[0 0],'Color','k','Parent',ax2);
set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2]);
xmaxtick = floor(xmax*10)/10;
set(ax2,'XTick',0.2:0.1:xmaxtick);
XTickLabels = cell(8+xmaxtick*10,1);
XTickLabels{1} = '0.2';
XTickLabels{3} = '0.4';
XTickLabels{6} = '0.7';
for j = 1:floor(xmaxtick)
    XTickLabels{j*10-1} = int2str(j);
end
set(ax2,'XTickLabel',XTickLabels);
set(ax2,'YTick',[]);

%secondary xlabel
if strcmp(plot_type,'presentation')
    %xlabel - landscape
    xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
elseif strcmp(plot_type,'paper')
    %xlabel - portrait
    xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
end

if strcmp(plot_type,'presentation')
    %print plot - landscape
    set(gcf,'PaperUnits','inches','PaperSize',[7 6],'PaperPosition',[0 0 7 6],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surface_lowtau_N1_only_landscape.png'],'-dpng');
elseif strcmp(plot_type,'paper')
    %print plot - portrait
    set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_surface_lowtau_N1_only_portrait.png'],'-dpng');
end


%% Plot taunorm-conditioned airborne versus surface size distributions - landscape - N1 only
figure(13); clf; hold on;

i = 3; %only N1
k_bins = 2:N_taunorm_bins; %skip first bin

%plot surface distribution
plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);

%plot airborne distribution
for k = k_bins
    plot(d_taunormbar_airborne_mid_drange_Cluster{i},...
        dVdlogd_taunormbar_airborne_drange_Cluster{i}(k,:),...
        ['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_taunorm_bin{k})
end

%format plot
set(gca,'XScale','log','YScale','log','YMinorTick','On');
xlim([0.06, 2]);
set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'});
ylim([1e-4 1e1]);

%label plot
htitle = title(ClusterNames{i});
set(htitle,'Position',[0.35,3.5]); %set title below edge of box to accommodate second axis

%labels 
if strcmp(plot_type,'presentation')
    %xlabel - landscape
    xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    %ylabel - landscape
    ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
elseif strcmp(plot_type,'paper')
    %xlabel - portrait
    xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    %ylabel - portrait
    ylabel('Non-dimensionalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
end

%create legend
legend_items = cell(length(k_bins)+1,1);
legend_items{1} = 'Surface';
for j = k_bins
    if j == 1
        legend_items{j} = ['\tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];       
    else
        legend_items{j} = [num2str(taunorm_min_bins(j),'%10.1f'),' < \tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];
    end
end
h_legend = legend(legend_items,'Location','SouthWest');
set(h_legend,'FontSize',10);

%add second axis
ax1 = gca; %get handle for first axis
set(ax1,'XColor','k','YColor','k');
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k',...
    'YColor','k');
cla;
xmin = 0.06/d50_bar_surface_Cluster(i);
xmax = 2/d50_bar_surface_Cluster(i);
line([xmin xmax],[0 0],'Color','k','Parent',ax2);
set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2]);
xmaxtick = floor(xmax*10)/10;
set(ax2,'XTick',0.2:0.1:xmaxtick);
XTickLabels = cell(8+xmaxtick*10,1);
XTickLabels{1} = '0.2';
XTickLabels{3} = '0.4';
XTickLabels{6} = '0.7';
for j = 1:floor(xmaxtick)
    XTickLabels{j*10-1} = int2str(j);
end
set(ax2,'XTickLabel',XTickLabels);
set(ax2,'YTick',[]);

%secondary xlabel
if strcmp(plot_type,'presentation')
    %xlabel - landscape
    xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
elseif strcmp(plot_type,'paper')
    %xlabel - portrait
    xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex')
end

if strcmp(plot_type,'presentation')
    %print plot - landscape
    set(gcf,'PaperUnits','inches','PaperSize',[7 6],'PaperPosition',[0 0 7 6],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_N1only_landscape.png'],'-dpng');
elseif strcmp(plot_type,'paper')
    %print plot - portrait
    set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
    print([folder_Plots,'GSD_taunorm_N1only_portrait.png'],'-dpng');
end


%% PLOT mean fratio / fsurface - unbinned - limited dhat
figure(4); clf; hold on;

%plot fratio values
for i = 1:N_Cluster  
    if strcmp(dref_type,'d50') 
        plot(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i),f_ratio_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    elseif strcmp(dref_type,'dmodal')
        plot(d_f_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),f_ratio_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    end
end

%plot modified fratio values
for i = 1:N_Cluster
    if strcmp(ClusterNames{i}(1:6),'Oceano')
        %get interpolated correction factors
        N_d = length(d_f_mid_Cluster{i});
        fC_4cm_interpolated = zeros(1,N_d)*NaN;
        fC_7cm_interpolated = zeros(1,N_d)*NaN;
        for j = 1:N_d
            prev_ind = find(d_mid_correction<=d_f_mid_Cluster{i}(j), 1, 'last' ); %index of previous correction factor
            next_ind = find(d_mid_correction>=d_f_mid_Cluster{i}(j), 1 ); %index of next correction factor
            next_wt = (log10(d_f_mid_Cluster{i}(j))-log10(d_mid_correction(prev_ind)))./...
                (log10(d_mid_correction(next_ind))-log10(d_mid_correction(prev_ind))); %weight of next correction fact
            prev_wt = 1-next_wt; %weight of previous correction factor
            if ~isempty(next_ind)
                fC_4cm_interpolated(j) = fC_4cm(prev_ind)*prev_wt+fC_4cm(next_ind)*next_wt;
                fC_7cm_interpolated(j) = fC_7cm(prev_ind)*prev_wt+fC_7cm(next_ind)*next_wt;
            end
        end  
        
        if strcmp(dref_type,'d50') 
            %plot(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i),f_ratio_bar_Cluster{i}.*fC_4cm_interpolated,'--.','Color',Color_Cluster{i}); %z > 4 cm correction
            plot(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i),f_ratio_bar_Cluster{i}.*fC_7cm_interpolated,'--','Color',Color_Cluster{i},'LineWidth',1); %z > 7 cm correction
        elseif strcmp(dref_type,'dmodal')
            %plot(d_f_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),f_ratio_bar_Cluster{i}.*fC_4cm_interpolated,'--.','Color',Color_Cluster{i}); %z > 4 cm correction
            plot(d_f_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),f_ratio_bar_Cluster{i}.*fC_7cm_interpolated,'--','Color',Color_Cluster{i},'LineWidth',1); %z > 7 cm correction
        end
    end
end    

%plot error bars
for i = 1:N_Cluster
    for k = 1:length(d_f_mid_Cluster{i})
        if strcmp(dref_type,'d50')
            plot(d_f_mid_Cluster{i}(k)./d50_bar_surface_Cluster(i)*[1 1],exp(log(f_ratio_bar_Cluster{i}(k))*[1 1]+f_ratio_sigmalog_Cluster{i}(k)*[-1 1]),'Color',Color_Cluster{i});
        elseif strcmp(dref_type,'dmodal')
            plot(d_f_mid_Cluster{i}(k)./dmodal_bar_surface_Cluster(i)*[1 1],exp(log(f_ratio_bar_Cluster{i}(k))*[1 1]+f_ratio_sigmalog_Cluster{i}(k)*[-1 1]),'Color',Color_Cluster{i});
        end
    end
end

%plot fine, medium, and coarse ranges
text(0.32,7,'Fine','HorizontalAlignment','Center','FontSize',13);
text(0.56,7,'Medium','HorizontalAlignment','Center','FontSize',13);
%text(0.95,7,'Coarse','HorizontalAlignment','Center','FontSize',13);
text(1.2,7,'Coarse','HorizontalAlignment','Center','FontSize',13);
plot([0.4 0.4],[1e-4 1e1],'k-.','LineWidth',0.5);
plot([0.8 0.8],[1e-4 1e1],'k-.','LineWidth',0.5);

%format plot
if strcmp(dref_type,'d50')==1
    xlabel('Non-dimensionalized grain size, $$d_i / d_{50,bed}$$','Interpreter','Latex','FontSize',14);
elseif strcmp(dref_type,'dmodal')==1
    xlabel('Non-dimensionalized grain size, $$d_i / d_{modal,bed}$$','Interpreter','Latex','FontSize',14);
end

%generate legend items
legend_items = ClusterNames;
for i = 1:N_Cluster
    if strcmp(ClusterNames{i}(1:6),'Oceano')
        legend_items{end+1} = [ClusterNames{i}, ' corrected'];
    end
end

% labels and legends
ylabel({'Ratio mean airborne to bed frac., $$\langle f_{i,air} \rangle/f_{i,bed}$$'},'Interpreter','Latex','HorizontalAlignment','Center'); %no inset plot
h_legend = legend(legend_items,'Location','SouthWest'); 
% if strcmp(plot_type,'presentation') %no inset plot
%     ylabel({'Ratio of airborne to bed fraction, $$\langle f_{i,air} \rangle/f_{i,bed}$$'},'Interpreter','Latex','HorizontalAlignment','Center'); %no inset plot
%     h_legend = legend(legend_items,'Location','SouthWest'); 
% elseif strcmp(plot_type,'paper') %with inset plot
%     ylabel({'Ratio of mean airborne to bed surface volume fraction, $$\langle f_{i,air} \rangle/f_{i,bed}$$'},'Interpreter','Latex','HorizontalAlignment','Center'); %original plot
%     h_legend = legend(legend_items,'Location','NorthEast'); 
%     text(0.26,7,'(a)','FontSize',PlotFont);
% end
set(h_legend,'FontSize',13);
set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',13);
ylim([1e-4 1e1]);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual'); %original plot
print([folder_Plots,'fratiobar_dhat_withcorrection_',dref_type,'.png'],'-dpng');

% %include inset plot with dimensional values, if for paper
% if strcmp(plot_type,'paper') 
%     axes('Position',[.25 .25 .35 .38]); hold on;
% 
%     %plot fratio values
%     for i = 1:N_Cluster   
%         plot(d_f_mid_Cluster{i},f_ratio_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
%     end
% 
%     %plot error bars
%     for i = 1:N_Cluster
%         for k = 1:length(d_f_mid_Cluster{i})
%             if strcmp(dref_type,'d50')
%                 plot(d_f_mid_Cluster{i}(k)*[1 1],exp(log(f_ratio_bar_Cluster{i}(k))*[1 1]+f_ratio_sigmalog_Cluster{i}(k)*[-1 1]),'Color',Color_Cluster{i});
%             elseif strcmp(dref_type,'dmodal')
%                 plot(d_f_mid_Cluster{i}(k)*[1 1],exp(log(f_ratio_bar_Cluster{i}(k))*[1 1]+f_ratio_sigmalog_Cluster{i}(k)*[-1 1]),'Color',Color_Cluster{i});
%             end
%         end
%     end
% 
%     xlabel('Grain size, $$d_{i}$$ (mm)','Interpreter','Latex');
%     ylabel('$$\langle f_{i,air} \rangle /f_{i,bed}$$','Interpreter','Latex');
%     set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
%     xlim([d_min d_max]);
%     ylim([1e-4 1e1]);
%     text(0.8,5,'(b)','FontSize',PlotFont);
% end
%     
% %print plot
% if strcmp(plot_type,'presentation')
%     set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual'); %no inset plot
%     print([folder_Plots,'fratiobar_dhat_noinset_',dref_type,'.png'],'-dpng');
% elseif strcmp(plot_type,'paper')
%     set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual'); %original plot
%     print([folder_Plots,'fratiobar_dhat_withinset_',dref_type,'.png'],'-dpng');
% end

%% PLOT correction factor from Namikas
figure(5); clf; hold on;

% create plot
plot(d_mid_correction, fC_4cm, 'b+--','LineWidth',1);
plot(d_mid_correction, fC_7cm, 'ro-.','LineWidth',1);
 
%format plot
set(gca,'XScale','log','YScale','log','YMinorTick','On','FontSize',12);
xlim([0.1, 1]);
set(gca,'xtick',[0.1:0.1:2]);
set(gca,'xticklabels',{'0.1','0.2','','0.4','','','0.7','','','1'},'FontSize',12);
ylim([0.5 20]);
xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex','FontSize',14)
ylabel('Correction factor, $$k_i$$','Interpreter','Latex','FontSize',14);
h_legend = legend('$$z \geq 4$$ cm','$$z \geq 7$$ cm','Location','NorthWest','Interpreter','Latex');
set(h_legend,'FontSize',14);

%add second axis
set(gca,'Position',[0.1300 0.1400 0.7750 0.7]) %reposition to make room for second label
ax1 = gca; %get handle for first axis
set(ax1,'XColor','k','YColor','k');
ax2 = axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k',...
    'YColor','k');
cla;
xmin = 0.06/d50_bar_surface_Cluster(i);
xmax = 2/d50_bar_surface_Cluster(i);
line([xmin xmax],[0 0],'Color','k','Parent',ax2);
set(ax2,'XScale','log','XLim',[xmin xmax],'YLim',[1 2]);
xmaxtick = floor(xmax*10)/10;
set(ax2,'XTick',0.2:0.1:xmaxtick);
XTickLabels = cell(8+xmaxtick*10,1);
XTickLabels{1} = '0.2';
XTickLabels{3} = '0.4';
XTickLabels{6} = '0.7';
for j = 1:floor(xmaxtick)
    XTickLabels{j*10-1} = int2str(j);
end
set(ax2,'XTickLabel',XTickLabels,'FontSize',12);
set(ax2,'YTick',[]);

%secondary xlabel
xlabel('Non-dimensionalized grain diameter, $$d/d_{50,bed}$$','Interpreter','Latex','FontSize',14)

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
print([folder_Plots,'GSD_correctionfactor.png'],'-dpng');

%% PLOT mean zqnorm VS dhat
figure(6); clf;

%% subplot for zq values
subplot('position',[0.08 0.56 0.9 0.42]); hold on;

%plot values
for i = 1:N_Cluster   
    if strcmp(dref_type,'d50') 
        plot(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i),zqi_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    elseif strcmp(dref_type,'dmodal')
        plot(d_f_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),zqi_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    end
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:length(d_f_mid_Cluster{i})
        if strcmp(dref_type,'d50')
            plot(d_f_mid_Cluster{i}(k)./d50_bar_surface_Cluster(i)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
        elseif strcmp(dref_type,'dmodal')
            plot(d_f_mid_Cluster{i}(k)./dmodal_bar_surface_Cluster(i)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
        end
    end
end

%plot fine, medium, and coarse ranges
text(0.32,0.01,'Fine','HorizontalAlignment','Center','FontSize',PlotFont);
text(0.56,0.01,'Medium','HorizontalAlignment','Center','FontSize',PlotFont);
text(1.0,0.01,'Coarse','HorizontalAlignment','Center','FontSize',PlotFont);
plot([0.4 0.4],[0 0.15],'k--','LineWidth',1);
plot([0.8 0.8],[0 0.15],'k--','LineWidth',1);

%format plot
xlim([0.28 2]);
ylim([0 0.15]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
if strcmp(dref_type,'dmodal')
    xlabel('Non-dimensionalized grain size, $$d_{i} / d_{modal}$$','Interpreter','Latex');
else
    xlabel('Non-dimensionalized grain size, $$d_{i} / d_{50,bed}$$','Interpreter','Latex');
end
ylabel('Mean size-selective saltation height, $$\langle z_{q,i} \rangle$$ (m)','Interpreter','Latex');
text(1.8,0.14,'(a)','FontSize',PlotFont);

%create legend
h_legend = legend(ClusterNames,'Location','NorthWest');
set(h_legend,'FontSize',12);

%% subplot for zq/dref values
subplot('position',[0.08 0.06 0.9 0.42]); hold on;

%plot zqnorm values
for i = 1:N_Cluster   
    if strcmp(dref_type,'d50') 
        plot(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i),zqinorm_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    elseif strcmp(dref_type,'dmodal')
        plot(d_f_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),zqinorm_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    end
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:length(d_f_mid_Cluster{i})
        if strcmp(dref_type,'d50')
            plot(d_f_mid_Cluster{i}(k)./d50_bar_surface_Cluster(i)*[1 1],zqinorm_bar_Cluster{i}(k)*[1 1]+zqinorm_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
        elseif strcmp(dref_type,'dmodal')
            plot(d_f_mid_Cluster{i}(k)./dmodal_bar_surface_Cluster(i)*[1 1],zqinorm_bar_Cluster{i}(k)*[1 1]+zqinorm_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
        end
    end
end

%plot mean zq/d50
for i = 1:N_Cluster
    plot([0.28 2],(1000*zq_bar_Cluster(i)/d50_bar_airborne_Cluster(i))*[1 1],'Color',Color_Cluster{i})
end

%plot fine, medium, and coarse ranges
plot([0.4 0.4],[0 450],'k--','LineWidth',1);
plot([0.8 0.8],[0 450],'k--','LineWidth',1);

%format plot
xlim([0.28 2]);
ylim([0 450]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
text(0.29,425,'(b)','FontSize',PlotFont);
if strcmp(dref_type,'dmodal')
    xlabel('Non-dimensionalized grain size, $$d_{i} / d_{modal}$$','Interpreter','Latex');
else
    xlabel('Non-dimensionalized grain size, $$d_{i} / d_{50,bed}$$','Interpreter','Latex');
end
ylabel('Mean non-dim. size-selective saltation height, $$\langle z_{q,i} \rangle/d_{i}$$','Interpreter','Latex');

%% inset plot with zq versus d
axes('Position',[.76 .328 .205 .18]); hold on;

%plot data
for i = 1:N_Cluster
    plot(d_f_mid_Cluster{i},zqi_bar_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:length(d_f_mid_Cluster{i})
        plot(d_f_mid_Cluster{i}(k)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end

%format plot
xlim([0.12 1.0]);
ylim([0 0.15]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont,'BoxStyle','Full');
set(gca,'xtick',[0.2:0.1:1]);
set(gca,'xticklabels',{'0.2','','0.4','','','0.7','','','1'});
xlabel('Grain size, $$d_{i}$$ (mm)','Interpreter','Latex','FontSize',10);
ylabel('$$\langle z_{q,i} \rangle$$ (m)','Interpreter','Latex','FontSize',10);
text(0.14,0.14,'(c)','FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 10],'PaperPosition',[0 0 8 10],'PaperPositionMode','Manual');
print([folder_Plots,'zq_d_',dref_type,'.png'],'-dpng');

% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % SUPPORTING INFORMATION PLOTS %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%% plot variation in znorm of BSNEs
figure(12); clf;
for i = 1:N_Cluster
    %initialize subplot
    if N_Cluster == 6 %defined subplot sizes for four clusters
        if i == 1
            h_subplot(1) = subplot('position',[0.08 0.72 0.40 0.26]); hold on;
        elseif i == 2
            h_subplot(2) = subplot('position',[0.57 0.72 0.40 0.26]); hold on;
        elseif i == 3
            h_subplot(3) = subplot('position',[0.08 0.39 0.40 0.26]); hold on;
        elseif i == 4
            h_subplot(4) = subplot('position',[0.57 0.39 0.40 0.26]); hold on;
        elseif i == 5
            h_subplot(5) = subplot('position',[0.08 0.06 0.40 0.26]); hold on;
        else
            h_subplot(6) = subplot('position',[0.57 0.06 0.40 0.26]); hold on;
        end
    else %otherwise, automated subplot sizes
        h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
    end   
    
    %plot znorm - usable profiles
    ind_usable = ind_usable_profile_Cluster{i};
    znorm_plot = znorm_profile_Cluster{i}(ind_usable);
    Time_plot = Time_profile_Cluster{i}(ind_usable);
    N_plot = length(znorm_plot);
    for j = 1:N_plot
        Time_profile = []; for k = 1:length(znorm_plot{j}); Time_profile = [Time_profile; Time_plot(j)]; end; %get times for profile
        h1 = plot(Time_profile,znorm_plot{j},'bo-');
    end
    
    %plot znorm - unusable profiles
    ind_unusable = setdiff(1:length(znorm_profile_Cluster{i}),ind_usable);
    znorm_plot = znorm_profile_Cluster{i}(ind_unusable);
    Time_plot = Time_profile_Cluster{i}(ind_unusable);
    N_plot = length(znorm_plot);
    for j = 1:N_plot
        Time_profile = []; for k = 1:length(znorm_plot{j}); Time_profile = [Time_profile; Time_plot(j)]; end; %get times for profile
        h2 = plot(Time_profile,znorm_plot{j},'rx-');
    end
    
    if mod(i,2) == 1
        ylabel('Non-dim. BSNE trap height, $$z/z_q$$','Interpreter','Latex','FontSize',PlotFont);
    end

    if i == 1
        h_legend = legend([h1 h2],{'usable','unusable'},'Location','North');
        set(h_legend,'FontSize',PlotFont);
    end
    
    title(ClusterNames{i},'FontSize',PlotFont);
    ylim([0 7]);
    if i==1
        set(gca,'XTick',[datetime(2014,11,13):duration(48,0,0):datetime(2014,11,19),datetime(2014,11,21)]);
    elseif i>=5
        set(gca,'XTick',[datetime(2015,5,23):duration(72,0,0):datetime(2015,6,2),datetime(2015,6,5)]);
    elseif i>=3
        set(gca,'XTick',[datetime(2015,5,15):duration(48,0,0):datetime(2015,5,23)]);
    elseif i>=2
        set(gca,'XTick',[datetime(2015,3,23):duration(24,0,0):datetime(2015,3,25)]);
    end
    datetick('x','mmm dd','keepticks')
    set(gca,'YMinorTick','On','Box','On','FontSize',PlotFont);
    xlims = xlim;
    text(xlims(1)+range(xlims)*0.03,6.5,Label_Cluster{i},'FontSize',PlotFont)
    
    %set(gca,'XTick',xlims(1):(floor(days(range(xlims))/3)*duration(24,0,0)):xlims(2))
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 10],'PaperPosition',[0 0 8 10],'PaperPositionMode','Manual');
print([folder_Plots,'znorm_profile_BSNE.png'],'-dpng');

% %% plot variation in reference zqnorm with saltation flux - dhat values
% figure(13); clf;
% 
% %initialize subplots
% h_subplot = gobjects(N_Cluster,1);
% 
% for i = 1:N_Cluster
%     %initialize subplot
%     if N_Cluster == 6 %defined subplot sizes for six clusters
%         if i == 1
%             h_subplot(1) = subplot('position',[0.10 0.71 0.395 0.25]); hold on;
%         elseif i == 2
%             h_subplot(2) = subplot('position',[0.58 0.71 0.395 0.25]); hold on;
%         elseif i == 3
%             h_subplot(3) = subplot('position',[0.10 0.385 0.395 0.25]); hold on;
%         elseif i == 4
%             h_subplot(4) = subplot('position',[0.58 0.385 0.395 0.25]); hold on;
%         elseif i == 5
%             h_subplot(5) = subplot('position',[0.10 0.06 0.395 0.25]); hold on;
%         else
%             h_subplot(6) = subplot('position',[0.58 0.06 0.395 0.25]); hold on;
%         end
%     else %otherwise, automated subplot sizes
%         h_subplot(i) = subplot(round((N_Cluster+1)/2),2,i); hold on;
%     end
%     
%     ind_usable = ind_usable_profile_Cluster{i};
%     h_fine_air = plot(taunorm_profile_Cluster{i}(ind_usable),zqinorm_dhat_fine_Cluster{i}(ind_usable),'v');
%     h_medium_air = plot(taunorm_profile_Cluster{i}(ind_usable),zqinorm_dhat_medium_Cluster{i}(ind_usable),'o');
%     h_coarse_air = plot(taunorm_profile_Cluster{i}(ind_usable),zqinorm_dhat_coarse_Cluster{i}(ind_usable),'^');
%     c_fine = get(h_fine_air,'Color');
%     c_medium = get(h_medium_air,'Color');
%     c_coarse = get(h_coarse_air,'Color');
%         
%     %plot error bars
%     for k = 1:length(ind_usable)
%         plot(taunorm_profile_Cluster{i}(ind_usable(k))*[1 1],zqinorm_dhat_fine_Cluster{i}(ind_usable(k))*[1 1]+sigma_zqinorm_dhat_fine_Cluster{i}(ind_usable(k))*[-1 1],'Color',c_fine);
%         plot(taunorm_profile_Cluster{i}(ind_usable(k))*[1 1],zqinorm_dhat_medium_Cluster{i}(ind_usable(k))*[1 1]+sigma_zqinorm_dhat_medium_Cluster{i}(ind_usable(k))*[-1 1],'Color',c_medium);
%         plot(taunorm_profile_Cluster{i}(ind_usable(k))*[1 1],zqinorm_dhat_coarse_Cluster{i}(ind_usable(k))*[1 1]+sigma_zqinorm_dhat_coarse_Cluster{i}(ind_usable(k))*[-1 1],'Color',c_coarse);
%     end
%     
%     %plot trends - fine
%     if abs(b_zqinorm_taunorm_fine(i)./sigma_b_zqinorm_taunorm_fine(i))>=2
%         plot(taunorm_fit_zqinorm_fine{i},zqinorm_pred_fine{i},'-','Color',c_fine,'LineWidth',1);
%     else
%         plot(taunorm_fit_zqinorm_fine{i},zqinorm_pred_fine{i},'-.','Color',c_fine,'LineWidth',1);
%     end
%     
%     %plot trends - medium
%     if abs(b_zqinorm_taunorm_medium(i)./sigma_b_zqinorm_taunorm_medium(i))>=2
%         plot(taunorm_fit_zqinorm_medium{i},zqinorm_pred_medium{i},'-','Color',c_medium,'LineWidth',1);
%     else
%         plot(taunorm_fit_zqinorm_medium{i},zqinorm_pred_medium{i},'-.','Color',c_medium,'LineWidth',1);
%     end
%     
%     %plot trends - coarse
%     if abs(b_zqinorm_taunorm_coarse(i)./sigma_b_zqinorm_taunorm_coarse(i))>=2
%         plot(taunorm_fit_zqinorm_coarse{i},zqinorm_pred_coarse{i},'-','Color',c_coarse,'LineWidth',1);
%     else
%         plot(taunorm_fit_zqinorm_coarse{i},zqinorm_pred_coarse{i},'-.','Color',c_coarse,'LineWidth',1);
%     end
%     
%     %format plot
%     xlim([1 4]);
%     ylim([0 400]);
%     ylims = ylim;
%     text(1.05, ylims(1)+range(ylims)*0.94, Label_Cluster{i},'FontSize',14);
%     if i==N_Cluster || i==N_Cluster-1
%         xlabel('Non-dimensionalized shear stress, $$\tau/\tau_{it}$$','Interpreter','Latex')
%     end
%     if mod(i,2) == 1
%         ylabel('Norm. size-sel. salt. ht., $$z_{q,i}/d_{i}$$','Interpreter','Latex');
%     end
%     if i == 2
%         h_legend = legend([h_fine_air, h_medium_air, h_coarse_air],...
%             ['d_{i} /d_{50,bed}=',num2str(dhat_fine)],...
%             ['d_{i} /d_{50,bed}=',num2str(dhat_medium)],...
%             ['d_{i} /d_{50,bed}=',num2str(dhat_coarse)],...
%             'Location','SouthWest');
%         set(h_legend,'FontSize',10);
%     end
%     title(ClusterNames{i});
%     set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[8 9],'PaperPosition',[0 0 8 9],'PaperPositionMode','Manual');
% print([folder_Plots,'SaltationHeight_IndexDHat_NormShearStress.png'],'-dpng');