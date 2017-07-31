%% initialize
clearvars;
close all;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA SOURCES AND LOADING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% information about where to load data and save plots
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
GrainSizeData_Path = strcat(folder_GrainSizeData,'GrainSizeData'); %path for loading mean grain size data
folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%% load data and functions
load(GrainSizeData_Path); %grain size data
addpath(folder_Functions); %point MATLAB to location of functions

%% get information about sites
N_Sites = length(Sites);

%% Cluster surface data by location
% Location_Cluster = {...
%     {{'BSNE_A','BSNE_B','Tower_A','Tower_B','U1','Upwind_X0m'},{'Upwind_X5m','Upwind_X10m','Upwind_X15m','Upwind_X20m'}},...
%     {{'A','B','C','Wenglor','A','B','C','Wenglor'}},...
%     {{'A','B','A1','A2','A3','A4'},{'Wenglor'},{'C','D','B1','B2_B3','B4'}},...
%     };
% Name_Cluster = {...
%     {{'0'},{'Upwind'}},...
%     {{'All'}},...
%     {{'A'},{'Wenglor'},{'B'}},...
%     };
Location_Cluster = {...
    {{'BSNE_A','BSNE_B','Tower_A','Tower_B','U1','Upwind_X0m','Upwind_X5m','Upwind_X10m','Upwind_X15m','Upwind_X20m'}},...
    {{'A','B','C','Wenglor','A','B','C','Wenglor'}},...
    {{'A','B','A1','A2','A3','A4'},{'C','D','B1','B2_B3','B4'}},...
    };
Name_Cluster = {...
    {{'Jericoacoara'}},...
    {{'Rancho Guadalupe'}},...
    {{'Oceano N'},{'Oceano S'}},...
    };
N_Cluster_Site = zeros(N_Sites,1);
for i = 1:N_Sites
    N_Cluster_Site(i) = length(Location_Cluster{i});
end
Label_Cluster = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'}; %markers for Clusters

%% go through sites, make figure
figure(1);
N_Cluster = sum(N_Cluster_Site);

for i = 1:N_Sites
       
    %% load relevant information about surface samples
    GrainSize_surface = GrainSizeData_all{i}.GrainSize_Surface; %surface sample array
    Location_surface = {GrainSize_surface.Location}; %locations of surface samples
    Date_surface = [GrainSize_surface.Date]; %dates of surface samples
    Time_surface = [GrainSize_surface.CollectionTime]; %times of surface samples
    d10_surface = [GrainSize_surface.d_10_mm]; %d10 of surface samples
    d50_surface = [GrainSize_surface.d_50_mm]; %d50 of surface samples
    d90_surface = [GrainSize_surface.d_90_mm]; %d90 of surface samples    
    
    for j = 1:N_Cluster_Site(i)
            
        %get values for cluster
        ind_Cluster = []; %inialize list of indices associated with cluster
        for k = 1:length(Location_Cluster{i}{j})
            ind_Cluster = [ind_Cluster, find(strcmp(Location_surface,Location_Cluster{i}{j}{k})==1)];
        end
        Time_Cluster = Time_surface(ind_Cluster); %get times for cluster
        [Time_Cluster, ind_sort] = sort(Time_Cluster); %sort times for cluster
        ind_Cluster = ind_Cluster(ind_sort); %sort associated cluster indices
        d10_Cluster = d10_surface(ind_Cluster); %get d10 values for cluster
        d50_Cluster = d50_surface(ind_Cluster); %get d50 values for cluster
        d90_Cluster = d90_surface(ind_Cluster); %get d90 values for cluster
        
        %% initialize subplot
        ind_subplot = sum(N_Cluster_Site(1:i))-N_Cluster_Site(i)+j;
        
        if N_Cluster == 4 %defined subplot sizes for four clusters
            if ind_subplot == 1
                subplot('position',[0.07 0.57 0.4 0.39]); hold on;
            elseif ind_subplot == 2
                subplot('position',[0.56 0.57 0.4 0.39]); hold on;
            elseif ind_subplot == 3
                subplot('position',[0.07 0.07 0.4 0.39]); hold on;
            else
                subplot('position',[0.56 0.07 0.4 0.39]); hold on;
            end
        else %otherwise, automated subplot sizes
            subplot(ceil(N_Cluster/2),2,ind_subplot); hold on;
        end
        
        
        %plot cluster values through time
        plot(Time_Cluster, d90_Cluster, '-^'); %plot d90
        plot(Time_Cluster, d50_Cluster, '-o'); %plot d50
        plot(Time_Cluster, d10_Cluster, '-v'); %plot d10
        if mod(ind_subplot,2)==1
            ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex','FontSize',12);
        end
        
        %plot dividing line for Oceano
        if i == 3
            plot([datetime(2015,5,23),datetime(2015,5,23)],[0 1],'k--');
        end
            
        %format plot
        title(Name_Cluster{i}{j},'FontSize',12);
        ylim([0 1]);
        set(gca,'XMinorTick','On','YMinorTick','On','Box','On','FontSize',12);
        xlim([min(Date_surface),max(Date_surface)+duration(24,0,0)])
        text(min(Date_surface)+range(xlim)*0.03,0.93,Label_Cluster{ind_subplot},'FontSize',12)

        datetick('x','mmm dd','keeplimits','keepticks')
        if ind_subplot == 1
            h_legend = legend('d_{90}','d_{50}','d_{10}','Location','North');
            set(h_legend,'FontSize',12);
        end
    end
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[9 6],'PaperPosition',[0 0 9 6],'PaperPositionMode','Manual');
print([folder_Plots,'SurfaceGrainSizeThroughTime.png'],'-dpng');