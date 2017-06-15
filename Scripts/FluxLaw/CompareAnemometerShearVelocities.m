%% COMPARE WINDS AT DIFFERENT ANEMOMETER HEIGHTS

%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%

%initialize
clearvars;
close all;

%parameters
g = 9.8; %gravity (m/s^2)
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
fQ_Qtau_fit_min = 0.1; %minimum mean fQ for Q versus tau comparison (and determination of threshold)
N_sigma_tauit = 2; %std deviations from impact threshold for minimum tauex

%info for plotting
Markers_Plot = {'s','d','o','<','>','^'};
Colors_Plot = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};

%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for retrieving data for this analysis
folder_SaveData = '../../AnalysisData/WindAnalysis/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/WindAnalysis/'; %folder for plots

%% paths for loading and saving data - unrestricted
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'FluxLawCalcs_30min_Unrestricted'); %path for 30 minute data

%% load data
load(LoadData_Path);
addpath(folder_Functions); %point MATLAB to location of functions

%% go through sites

%initialize plots
figure(1); clf; %all anemometers
figure(2); clf; %lowest anemometers

for i = 1:N_Sites
    
    %% get anemometer information
    N_anemometers = size(zU_profile_all{i},2); %number of anemometers
    ind_lowest = find(~isnan(ustRe_profile_all{i}(:,1))); %indices for which lowest anemometer is valid
    
    %get list of indices for plotting
    ind_plot = cell(N_anemometers,1);
    for j = 2:N_anemometers
        ind_plot{j} = intersect(ind_lowest, ...
            find(~isnan(ustRe_profile_all{i}(:,j)) &...
            zU_profile_all{i}(:,j)>zU_profile_all{i}(:,1)));
    end
    
    %get mean heights for anemometers and labels for legend
    zbar_anemometer = zeros(N_anemometers,1);
    zbar_anemometer(1) = round(mean(zU_profile_all{i}(ind_lowest,1)),2);
    legend_labels = cell(N_anemometers-1,1);
    legend_labels{1} = ['1-1 line'];
    for j = 2:N_anemometers
        zbar_anemometer(j) = round(mean(zU_profile_all{i}(ind_plot{j},j)),2);
        legend_labels{j} = ['z = ',num2str(zbar_anemometer(j)),' m'];
    end
    
    %% generate plot - all
    figure(1); subplot(1,N_Sites,i); hold on;
    plot([0 max(max(ustRe_profile_all{i}))],[0 max(max(ustRe_profile_all{i}))],'k--'); %1-1 line
    for j = 2:N_anemometers %main points
        plot(ustRe_profile_all{i}(ind_plot{j},1),ustRe_profile_all{i}(ind_plot{j},j),Markers_Plot{j-1},'Color',Colors_Plot{j-1},'MarkerSize',3);
    end
    
    %error bars
    for j = 2:N_anemometers %error bars
        for k = 1:ind_plot{j}
            plot(ustRe_profile_all{i}(k,1)+[1 -1]*sigma_ustRe_profile_all{i}(k,1),ustRe_profile_all{i}(k,j)*[1 1],'Color',Colors_Plot{j-1});
            plot(ustRe_profile_all{i}(k,1)*[1 1],ustRe_profile_all{i}(k,j)+[1 -1]*sigma_ustRe_profile_all{i}(k,j),'Color',Colors_Plot{j-1});
        end
    end 

    %format plot
    xlabel(['u_{*} (m/s) for low anemometer, z = ',num2str(zbar_anemometer(1)),' m']);
    ylabel('u_{*} (m/s) for high anemometer');
    legend(legend_labels,'Location','SouthEast');
    title(SiteNames{i});
    
    %print plot
    set(gca,'LooseInset',get(gca,'TightInset'),'Box','On','FontSize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 5]);
    print([folder_Plots,'ust_comparison_all.png'],'-dpng');
    
    %% generate plot - lowest anemometers
    figure(2); subplot(1,N_Sites,i); hold on;
    plot([0 max(max(ustRe_profile_all{i}))],[0 max(max(ustRe_profile_all{i}))],'k--'); %1-1 line
    plot(ustRe_profile_all{i}(ind_plot{2},1),ustRe_profile_all{i}(ind_plot{2},2),Markers_Plot{2-1},'Color',Colors_Plot{2-1}); %main points
    
    %error bars
    for k = 1:ind_plot{2}
        plot(ustRe_profile_all{i}(k,1)+[1 -1]*sigma_ustRe_profile_all{i}(k,1),ustRe_profile_all{i}(k,2)*[1 1],'Color',Colors_Plot{2-1});
        plot(ustRe_profile_all{i}(k,1)*[1 1],ustRe_profile_all{i}(k,2)+[1 -1]*sigma_ustRe_profile_all{i}(k,2),'Color',Colors_Plot{2-1});
    end

    %format plot
    xlabel(['u_{*} (m/s) for low anemometer, z = ',num2str(zbar_anemometer(1)),' m']);
    ylabel(['u_{*} (m/s) for second anemometer, z = ',num2str(zbar_anemometer(2)),' m']);
    legend(legend_labels{1},'Location','SouthEast');
    title(SiteNames{i});
    
    %print plot
    set(gca,'LooseInset',get(gca,'TightInset'),'Box','On','FontSize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 5]);
    print([folder_Plots,'ust_comparison_lowest.png'],'-dpng');
end