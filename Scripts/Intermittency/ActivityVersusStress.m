%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;
close all;

%% parameters
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
theta_max = 20; %maximum absolute wind angle for calcs
zL_max = 0.2; %maximum absolute stability value for calcs

%% binning info - tau
tau_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
tau_bin_N_min = 3; %minimum number of entries for bin

%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for retrieving data for this analysis
folder_LoadSecondaryData = '../../AnalysisData/Thresholds/'; %folder for retrieving data for this analysis
folder_SaveData = '../../AnalysisData/Intermittency/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/Intermittency/'; %folder for plots

% %% paths for loading and saving data - restricted
% LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'RoughnessCalcs_30min_Restricted'); %path for 30 minute data

%% paths for loading and saving data - unrestricted
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data
LoadSecondaryData_Path = strcat(folder_LoadSecondaryData,'ThresholdAnalysisData'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'IntermittencyCalcs_30min_Unrestricted'); %path for 30 minute data

%% load data
load(LoadData_Path);
load(LoadSecondaryData_Path);
addpath(folder_Functions); %point MATLAB to location of functions

%% plotting info
PlotFont = 12;
PlotMarkers_Site = {'s','d','o'};
PlotColors_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

%% get 30 minute window wind directions adjusted for predominant wind during continuous transport
theta_adjusted_all = cell(N_Sites,1); %wind direction
for i = 1:N_Sites
    theta_Site = mean(theta_all{i}(fQ_all{i}==1)); %wind direction for site, based on mean wind during continuous transport
    theta_adjusted_all{i} = theta_all{i} - theta_Site;
end

%% for Jeri, convert fQ NaN values to zeros, assuming that times with no flux data are fQ = 0 (not valid for RG and not necessary for Oceano)
ind_Site = find(strcmp(SiteNames,'Jericoacoara'));
fQ_all{ind_Site}(isnan(fQ_all{ind_Site}))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMARY ANALYSIS BY SITE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get 30-minute binned values

%initialize binned variables - tau
tau_bin_avg_all = cell(N_Sites,1);
tau_bin_SE_all = cell(N_Sites,1);
fQ_tau_bin_avg_all = cell(N_Sites,1);
fQ_tau_bin_SE_all = cell(N_Sites,1);

for i = 1:length(Sites)
    
    %get indices for binning by tau
    ind_theta = find(abs(theta_adjusted_all{i})<=theta_max); %indices for theta range
    ind_zL = find(abs(zL_all{i})<=zL_max); %indices for stability range
    ind_wind = intersect(ind_theta,ind_zL); %indices for ok wind
    ind_fQ = find(~isnan(fQ_all{i})); %indices of ok fQ
    ind_binning = intersect(ind_wind,ind_fQ); %indices for binning
    
    %get values for binning
    tau_binning = tauRe_all{i}(ind_binning);
    fQ_binning = fQ_all{i}(ind_binning);

    %get binned values for tau
    [~, ~, tau_bin_min, tau_bin_max, tau_bin_avg, tau_bin_SE] = Binning(tau_binning, tau_bin_minrange, tau_bin_maxrange, tau_bin_N_min);
    N_tau_bins = length(tau_bin_min);
    tau_bin_avg_all{i} = tau_bin_avg;
    tau_bin_SE_all{i} = tau_bin_SE;

    %get binned values for fQ
    fQ_tau_bin_avg_all{i} = zeros(N_tau_bins,1);
    fQ_tau_bin_SE_all{i} = zeros(N_tau_bins,1);
    for k=1:N_tau_bins
        tau_bin_ind = find(tau_binning>=tau_bin_min(k)&tau_binning<=tau_bin_max(k));
        fQ_tau_bin_avg_all{i}(k) = mean(fQ_binning(tau_bin_ind));
        fQ_tau_bin_SE_all{i}(k) = std(fQ_binning(tau_bin_ind))/sqrt(length(tau_bin_ind));
    end
end

%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%
save(SaveData_Path,'*all');

%%%%%%%%%
% PLOTS %
%%%%%%%%%

%% plot fQ versus tau
figure(2); clf; hold on; %initialize plot

%plot fQ average values, impact threshold, and fluid threshold
for i = 1:N_Sites
    plot(tau_bin_avg_all{i},fQ_tau_bin_avg_all{i},PlotMarkers_Site{i},'Color',PlotColors_Site{i});
    plot(tauit_all(i)*[1 1],[0 0.1],'--','Color',PlotColors_Site{i},'LineWidth',2); %vertical lines for impact threshold
    plot(tauft_all(i)*[1 1],[0.9 1],'-.','Color',PlotColors_Site{i},'LineWidth',2); %vertical lines for fluid threshold
end

%plot fQ SE values
for i = 1:N_Sites
    N_tau_bins = length(tau_bin_avg_all{i});
    for k = 1:N_tau_bins
        plot(tau_bin_avg_all{i}(k)+tau_bin_SE_all{i}(k)*[-1 1],fQ_tau_bin_avg_all{i}(k)*[1 1],'Color',PlotColors_Site{i});
        plot(tau_bin_avg_all{i}(k)*[1 1],fQ_tau_bin_avg_all{i}(k)+fQ_tau_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors_Site{i});
    end
end

%generate legend items
legend_items = cell(N_Sites*3,1);
for i = 1:N_Sites
    legend_items{3*i-2}=SiteNames{i};
    legend_items{3*i-1}=['\tau_{it}'];
    legend_items{3*i}=['\tau_{ft}'];
end

%annotate plot
ylim([0 1]);
xlim([0 0.3]);
legend(legend_items,'Location','SouthEast');
xlabel('shear stress, $$\tau$$ (Pa)','interpreter','latex');
ylabel('saltation activity, $$f_{Q}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
print([folder_Plots,'fQ_tau.png'],'-dpng');
print([folder_Plots,'fQ_tau.tif'],'-dtiff','-r600');