%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;
close all;

%% parameters
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
fQ_max_z0 = 0.05; %maximum activity for z0 calc
Q_max_fit = 30; %maximum Q for zs fit
theta_max = 20; %maximum absolute wind angle for calcs
zL_max = 0.2; %maximum absolute stability value for calcs

%% binning info - fQ
fQ_min = 0.05;
fQ_max = 0.95;
fQ_bin_minrange = 0.1;
fQ_bin_maxrange = 0.2;
fQ_bin_N_min = 3;

%% binning info - tau
tau_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
tau_bin_N_min = 3; %minimum number of entries for bin

%% binning info - Q
Q_min = 2; %mininum Q for zs fitting
Q_max = 60;
Q_bin_minrange = 1; %minimum difference between upper and lower value in bin
Q_bin_maxrange = 3; %maximum difference between upper and lower value in bin
Q_bin_N_min = 3; %minimum number of entries for bin

%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for retrieving data for this analysis
folder_SaveData = '../../AnalysisData/Misc/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/Roughness/'; %folder for plots

%% paths for loading and saving data - restricted
LoadProfileData_Path = strcat(folder_LoadData,'LogProfiles_30min_Restricted'); %path for 30 minute data on log profiles
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'RoughnessCalcs_30min_Restricted'); %path for 30 minute data

% %% paths for loading and saving data - unrestricted
% LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'RoughnessCalcs_30min_Unrestricted'); %path for 30 minute data

%% load data
load(LoadProfileData_Path);
load(LoadData_Path);
addpath(folder_Functions); %point MATLAB to location of functions

%% plotting info
PlotFont = 12;
PlotMarkers_Site = {'s','d','o'};
PlotColors_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

%% get 30 minute and subwindow wind directions adjusted for predominant wind during continuous transport
theta_adjusted_all = cell(N_Sites,1); %wind direction
for i = 1:N_Sites
    theta_Site = mean(theta_all{i}(fQ_all{i}==1)); %wind direction for site, based on mean wind during continuous transport
    theta_adjusted_all{i} = theta_all{i} - theta_Site;
end

%% for Jeri, convert fQ and Q NaN values to zeros, assuming that times with no flux data are Q = fQ = 0 (not valid for RG and not necessary for Oceano)
ind_Site = find(strcmp(SiteNames,'Jericoacoara'));
fQ_all{ind_Site}(isnan(fQ_all{ind_Site}))=0;
Q_all{ind_Site}(isnan(Q_all{ind_Site}))=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMARY ANALYSIS BY SITE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get 30-minute binned values for zs and perform fits versus fQ

%initialize binned variables - fQ
fQ_bin_avg_all = cell(N_Sites,1);
fQ_bin_SE_all = cell(N_Sites,1);
zsRe_ln_fQ_bin_avg_all = cell(N_Sites,1);
zsRe_ln_fQ_bin_SE_all = cell(N_Sites,1);
zsLog_ln_fQ_bin_avg_all = cell(N_Sites,1);
zsLog_ln_fQ_bin_SE_all = cell(N_Sites,1);
z0Re_ln_fQ_bin_avg_all = zeros(N_Sites,1);
z0Re_ln_fQ_bin_SE_all = zeros(N_Sites,1);
z0Log_ln_fQ_bin_avg_all = zeros(N_Sites,1);
z0Log_ln_fQ_bin_SE_all = zeros(N_Sites,1);

%initialize binned variables - tau
tau_bin_avg_all = cell(N_Sites,1);
tau_bin_SE_all = cell(N_Sites,1);
zsRe_ln_tau_bin_avg_all = cell(N_Sites,1);
zsRe_ln_tau_bin_SE_all = cell(N_Sites,1);
zsLog_ln_tau_bin_avg_all = cell(N_Sites,1);
zsLog_ln_tau_bin_SE_all = cell(N_Sites,1);

%initialize binned variables - Q
Q_bin_avg_all = cell(N_Sites,1);
Q_bin_SE_all = cell(N_Sites,1);
zsRe_ln_Q_bin_avg_all = cell(N_Sites,1);
zsRe_ln_Q_bin_SE_all = cell(N_Sites,1);
zsLog_ln_Q_bin_avg_all = cell(N_Sites,1);
zsLog_ln_Q_bin_SE_all = cell(N_Sites,1);
z0Re_ln_Q_bin_avg_all = zeros(N_Sites,1);
z0Re_ln_Q_bin_SE_all = zeros(N_Sites,1);
z0Log_ln_Q_bin_avg_all = zeros(N_Sites,1);
z0Log_ln_Q_bin_SE_all = zeros(N_Sites,1);

%initialize fit values - fQ
z0Re_fQ_fit = zeros(N_Sites,1);
z0Re_ln_fQ_fit = zeros(N_Sites,1);
sigma_z0Re_ln_fQ_fit = zeros(N_Sites,1);
dlnzsRe_dfQ_fit = zeros(N_Sites,1); %slope of fit line
sigma_dlnzsRe_dfQ_fit = zeros(N_Sites,1); %slope of fit line
z0Log_fQ_fit = zeros(N_Sites,1);
z0Log_ln_fQ_fit = zeros(N_Sites,1);
sigma_z0Log_ln_fQ_fit = zeros(N_Sites,1);
dlnzsLog_dfQ_fit = zeros(N_Sites,1); %slope of fit line
sigma_dlnzsLog_dfQ_fit = zeros(N_Sites,1); %slope of fit line

%initialize fit values - Q
z0Re_Q_fit = zeros(N_Sites,1);
z0Re_ln_Q_fit = zeros(N_Sites,1);
sigma_z0Re_ln_Q_fit = zeros(N_Sites,1);
dlnzsRe_dQ_fit = zeros(N_Sites,1); %slope of fit line
sigma_dlnzsRe_dQ_fit = zeros(N_Sites,1); %slope of fit line
z0Log_Q_fit = zeros(N_Sites,1);
z0Log_ln_Q_fit = zeros(N_Sites,1);
sigma_z0Log_ln_Q_fit = zeros(N_Sites,1);
dlnzsLog_dQ_fit = zeros(N_Sites,1); %slope of fit line
sigma_dlnzsLog_dQ_fit = zeros(N_Sites,1); %slope of fit line

for i = 1:length(Sites)
    
    %% binning by fQ
    %get indices for binning
    ind_theta = find(abs(theta_adjusted_all{i})<=theta_max); %indices for theta range
    ind_zL = find(abs(zL_all{i})<=zL_max); %indices for stability range
    ind_wind = intersect(ind_theta,ind_zL); %indices for ok wind
    ind_fQ = intersect(find(fQ_all{i}>=fQ_min),find(fQ_all{i}<=fQ_max)); %indices for fQ range
    ind_binning = intersect(ind_wind,ind_fQ); %indices for binning
    
    %perform binning of zs by fQ if there are values to bin
    if(~isempty(ind_binning))
        
        %get values for binning
        fQ_binning = fQ_all{i}(ind_binning);
        zsRe_binning = zsRe_all{i}(ind_binning);
        zsLog_binning = zsLog_all{i}(ind_binning);
        
        %get binned values for fQ
        [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
        N_fQ_bins = length(fQ_bin_min);
        fQ_bin_avg_all{i} = fQ_bin_avg;
        fQ_bin_SE_all{i} = fQ_bin_SE;
        
        %get binned values for zsRe
        zsRe_ln_fQ_bin_avg_all{i} = zeros(N_fQ_bins,1);
        zsRe_ln_fQ_bin_SE_all{i} = zeros(N_fQ_bins,1);
        for k=1:N_fQ_bins
            fQ_bin_ind = find(fQ_binning>=fQ_bin_min(k)&fQ_binning<=fQ_bin_max(k));
            zsRe_ln_fQ_bin_avg_all{i}(k) = mean(log(zsRe_binning(fQ_bin_ind)));
            zsRe_ln_fQ_bin_SE_all{i}(k) = std(log(zsRe_binning(fQ_bin_ind)))/sqrt(length(fQ_bin_ind));
        end
        
        %get binned values for zsLog
        zsLog_ln_fQ_bin_avg_all{i} = zeros(N_fQ_bins,1);
        zsLog_ln_fQ_bin_SE_all{i} = zeros(N_fQ_bins,1);
        for k=1:N_fQ_bins
            fQ_bin_ind = find(fQ_binning>=fQ_bin_min(k)&fQ_binning<=fQ_bin_max(k));
            zsLog_ln_fQ_bin_avg_all{i}(k) = mean(log(zsLog_binning(fQ_bin_ind)));
            zsLog_ln_fQ_bin_SE_all{i}(k) = std(log(zsLog_binning(fQ_bin_ind)))/sqrt(length(fQ_bin_ind));
        end
    end
    
    %get binned values for z0 calculation
    ind_fQ = find(fQ_all{i}<=fQ_max_z0);
    ind_binning = intersect(ind_wind,ind_fQ);
    
    %compute mean and SE of z0Re values, add to list
    z0Re_ln_fQ_bin_avg_all(i) = mean(log(zsRe_all{i}(ind_binning)));
    z0Re_ln_fQ_bin_SE_all(i) = std(log(zsRe_all{i}(ind_binning)))/sqrt(length(ind_binning));
    
    %compute mean and SE of z0Log values, add to list
    z0Log_ln_fQ_bin_avg_all(i) = mean(log(zsLog_all{i}(ind_binning)));
    z0Log_ln_fQ_bin_SE_all(i) = std(log(zsLog_all{i}(ind_binning)))/sqrt(length(ind_binning));

    %% binning by tau
    %get indices for binning
    ind_binning = ind_wind; %indices for binning
    
    %get values for binning
    tau_binning = tauRe_all{i}(ind_binning);
    zsRe_binning = zsRe_all{i}(ind_binning);
    zsLog_binning = zsLog_all{i}(ind_binning);

    %get binned values for fQ
    [~, ~, tau_bin_min, tau_bin_max, tau_bin_avg, tau_bin_SE] = Binning(tau_binning, tau_bin_minrange, tau_bin_maxrange, tau_bin_N_min);
    N_tau_bins = length(tau_bin_min);
    tau_bin_avg_all{i} = tau_bin_avg;
    tau_bin_SE_all{i} = tau_bin_SE;

    %get binned values for zsRe
    zsRe_ln_tau_bin_avg_all{i} = zeros(N_tau_bins,1);
    zsRe_ln_tau_bin_SE_all{i} = zeros(N_tau_bins,1);
    for k=1:N_tau_bins
        tau_bin_ind = find(tau_binning>=tau_bin_min(k)&tau_binning<=tau_bin_max(k));
        zsRe_ln_tau_bin_avg_all{i}(k) = mean(log(zsRe_binning(tau_bin_ind)));
        zsRe_ln_tau_bin_SE_all{i}(k) = std(log(zsRe_binning(tau_bin_ind)))/sqrt(length(tau_bin_ind));
    end

    %get binned values for zsLog
    zsLog_ln_tau_bin_avg_all{i} = zeros(N_tau_bins,1);
    zsLog_ln_tau_bin_SE_all{i} = zeros(N_tau_bins,1);
    for k=1:N_tau_bins
        tau_bin_ind = find(tau_binning>=tau_bin_min(k)&tau_binning<=tau_bin_max(k));
        zsLog_ln_tau_bin_avg_all{i}(k) = mean(log(zsLog_binning(tau_bin_ind)));
        zsLog_ln_tau_bin_SE_all{i}(k) = std(log(zsLog_binning(tau_bin_ind)))/sqrt(length(tau_bin_ind));
    end
    
    %% binning by Q
    %get indices for binning
    ind_Q = intersect(find(Q_all{i}>=Q_min),find(Q_all{i}<=Q_max)); %indices for Q > Q_min
    ind_binning = intersect(ind_wind,ind_Q); %indices for binning
    
    %get values for binning
    Q_binning = Q_all{i}(ind_binning);
    zsRe_binning = zsRe_all{i}(ind_binning);
    zsLog_binning = zsLog_all{i}(ind_binning);

    %get binned values for fQ
    [~, Q_bin_N, Q_bin_min, Q_bin_max, Q_bin_avg, Q_bin_SE] = Binning(Q_binning, Q_bin_minrange, Q_bin_maxrange, Q_bin_N_min);
    N_Q_bins = length(Q_bin_min);
    Q_bin_avg_all{i} = Q_bin_avg;
    Q_bin_SE_all{i} = Q_bin_SE;

    %get binned values for zsRe
    zsRe_ln_Q_bin_avg_all{i} = zeros(N_Q_bins,1);
    zsRe_ln_Q_bin_SE_all{i} = zeros(N_Q_bins,1);
    for k=1:N_Q_bins
        Q_bin_ind = find(Q_binning>=Q_bin_min(k)&Q_binning<=Q_bin_max(k));
        zsRe_ln_Q_bin_avg_all{i}(k) = mean(log(zsRe_binning(Q_bin_ind)));
        zsRe_ln_Q_bin_SE_all{i}(k) = std(log(zsRe_binning(Q_bin_ind)))/sqrt(length(Q_bin_ind));
    end
    
    %get binned values for zsLog
    zsLog_ln_Q_bin_avg_all{i} = zeros(N_Q_bins,1);
    zsLog_ln_Q_bin_SE_all{i} = zeros(N_Q_bins,1);
    for k=1:N_Q_bins
        Q_bin_ind = find(Q_binning>=Q_bin_min(k)&Q_binning<=Q_bin_max(k));
        zsLog_ln_Q_bin_avg_all{i}(k) = mean(log(zsLog_binning(Q_bin_ind)));
        zsLog_ln_Q_bin_SE_all{i}(k) = std(log(zsLog_binning(Q_bin_ind)))/sqrt(length(Q_bin_ind));
    end

    %adjust uncertainties for bins with only a single value
    ind_bin_full = find(Q_bin_N >= Q_bin_N_min);
    ind_bin_one = find(Q_bin_N == 1);
    zsRe_relerr_full = mean((zsRe_ln_Q_bin_SE_all{i}(ind_bin_full).*sqrt(Q_bin_N(ind_bin_full)))./zsRe_ln_Q_bin_avg_all{i}(ind_bin_full));
    zsRe_ln_Q_bin_SE_all{i}(ind_bin_one) = zsRe_relerr_full*zsRe_ln_Q_bin_avg_all{i}(ind_bin_one);
    zsLog_relerr_full = mean((zsLog_ln_Q_bin_SE_all{i}(ind_bin_full).*sqrt(Q_bin_N(ind_bin_full)))./zsLog_ln_Q_bin_avg_all{i}(ind_bin_full));
    zsLog_ln_Q_bin_SE_all{i}(ind_bin_one) = zsLog_relerr_full*zsLog_ln_Q_bin_avg_all{i}(ind_bin_one);
    
    %get binned values for z0 calculation
    ind_Q = find(Q_all{i}<=Q_min);
    ind_binning = intersect(ind_wind,ind_Q);
    
    %compute mean and SE of z0Re values, add to list
    z0Re_ln_Q_bin_avg_all(i) = mean(log(zsRe_all{i}(ind_binning)));
    z0Re_ln_Q_bin_SE_all(i) = std(log(zsRe_all{i}(ind_binning)))/sqrt(length(ind_binning));
    
    %compute mean and SE of z0Log values, add to list
    z0Log_ln_Q_bin_avg_all(i) = mean(log(zsLog_all{i}(ind_binning)));
    z0Log_ln_Q_bin_SE_all(i) = std(log(zsLog_all{i}(ind_binning)))/sqrt(length(ind_binning));
    
    %% fitting zs versus fQ
    %fQ values for both fits
    fQ_fit = [0; fQ_bin_avg_all{i}];
    
    %get values for fitting - Re
    zsRe_ln_fQ_fit = [z0Re_ln_fQ_bin_avg_all(i); zsRe_ln_fQ_bin_avg_all{i}];
    sigma_zsRe_ln_fQ_fit = [z0Re_ln_fQ_bin_SE_all(i); zsRe_ln_fQ_bin_SE_all{i}];
    ind_fQ_fit_Re = find(~isnan(zsRe_ln_fQ_fit)); %remove NaN values
    zsRe_ln_fQ_fit = zsRe_ln_fQ_fit(ind_fQ_fit_Re);
    sigma_zsRe_ln_fQ_fit = sigma_zsRe_ln_fQ_fit(ind_fQ_fit_Re);
    fQ_fit_Re = fQ_fit(ind_fQ_fit_Re);
    
    %get values for fitting - Log
    zsLog_ln_fQ_fit = [z0Log_ln_fQ_bin_avg_all(i); zsLog_ln_fQ_bin_avg_all{i}];
    sigma_zsLog_ln_fQ_fit = [z0Log_ln_fQ_bin_SE_all(i); zsLog_ln_fQ_bin_SE_all{i}];
    ind_fQ_fit_Log = find(~isnan(zsLog_ln_fQ_fit)); %remove NaN values
    zsLog_ln_fQ_fit = zsLog_ln_fQ_fit(ind_fQ_fit_Log);
    sigma_zsLog_ln_fQ_fit = sigma_zsLog_ln_fQ_fit(ind_fQ_fit_Log);
    fQ_fit_Log = fQ_fit(ind_fQ_fit_Log);
    
    %perform linear fit, compute and save z0, z0 uncertainty, and fit slope - Re
    [a, b, sigma_a, sigma_b] = linearfit(fQ_fit_Re, zsRe_ln_fQ_fit, sigma_zsRe_ln_fQ_fit);
    z0Re_fQ_fit(i) = exp(a);
    z0Re_ln_fQ_fit(i) = a;
    sigma_z0Re_ln_fQ_fit(i) = sigma_a;
    dlnzsRe_dfQ_fit(i) = b %display values for paper
    sigma_dlnzsRe_dfQ_fit(i) = sigma_b %display values for paper
    
    %perform linear fit, compute and save z0, z0 uncertainty, and fit slope - Log
    [a, b, sigma_a, sigma_b] = linearfit(fQ_fit_Log, zsLog_ln_fQ_fit, sigma_zsLog_ln_fQ_fit);
    z0Log_fQ_fit(i) = exp(a);
    z0Log_ln_fQ_fit(i) = a;
    sigma_z0Log_ln_fQ_fit(i) = sigma_a;
    dlnzsLog_dfQ_fit(i) = b;
    sigma_dlnzsLog_dfQ_fit(i) = sigma_b;
    
    %% fitting zs versus Q
    ind_Q_fit = find(Q_bin_avg_all{i}<Q_max_fit); %indices of data points to fit
    Q_fit = Q_bin_avg_all{i}(ind_Q_fit); %zs values only for nonzero Q
    
    %get values for fitting - Re
    zsRe_ln_Q_fit = zsRe_ln_Q_bin_avg_all{i}(ind_Q_fit); %zs values only for nonzero Q
    sigma_zsRe_ln_Q_fit = zsRe_ln_Q_bin_SE_all{i}(ind_Q_fit); %zs values only for nonzero Q
    ind_Q_fit_Re = find(~isnan(zsRe_ln_Q_fit)); %remove NaN values
    zsRe_ln_Q_fit = zsRe_ln_Q_fit(ind_Q_fit_Re);
    sigma_zsRe_ln_Q_fit = sigma_zsRe_ln_Q_fit(ind_Q_fit_Re);
    Q_fit_Re = Q_fit(ind_Q_fit_Re);
    
    %get values for fitting - Log
    zsLog_ln_Q_fit = zsLog_ln_Q_bin_avg_all{i}(ind_Q_fit); %zs values only for nonzero Q
    sigma_zsLog_ln_Q_fit = zsLog_ln_Q_bin_SE_all{i}(ind_Q_fit); %zs values only for nonzero Q
    ind_Q_fit_Log = find(~isnan(zsLog_ln_Q_fit)); %remove NaN values
    zsLog_ln_Q_fit = zsLog_ln_Q_fit(ind_Q_fit_Log);
    sigma_zsLog_ln_Q_fit = sigma_zsLog_ln_Q_fit(ind_Q_fit_Log);
    Q_fit_Log = Q_fit(ind_Q_fit_Log);
    
    %perform linear fit, compute and save z0, z0 uncertainty, and fit slope - Re
    [a, b, sigma_a, sigma_b] = linearfit(Q_fit_Re, zsRe_ln_Q_fit, sigma_zsRe_ln_Q_fit);
    z0Re_Q_fit(i) = exp(a) %display values for paper
    z0Re_ln_Q_fit(i) = a;
    sigma_z0Re_ln_Q_fit(i) = sigma_a %display values for paper
    dlnzsRe_dQ_fit(i) = b;
    sigma_dlnzsRe_dQ_fit(i) = sigma_b;
    
    %perform linear fit, compute and save z0, z0 uncertainty, and fit slope - Log
    [a, b, sigma_a, sigma_b] = linearfit(Q_fit_Log, zsLog_ln_Q_fit, sigma_zsLog_ln_Q_fit);
    z0Log_Q_fit(i) = exp(a);
    z0Log_ln_Q_fit(i) = a;
    sigma_z0Log_ln_Q_fit(i) = sigma_a;
    dlnzsLog_dQ_fit(i) = b;
    sigma_dlnzsLog_dQ_fit(i) = sigma_b;
end

%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%
save(SaveData_Path,'*fit');

%%%%%%%%%
% PLOTS %
%%%%%%%%%

%% plot zs versus fQ
figure(1); clf; %initialize plot

%% Reynolds values - first as separate plot, then joint plot
for n = 1:2
    if n==1
        hold on;
    elseif n==2
        %print plot
        set(gca, 'LooseInset', get(gca,'TightInset'));
        set(gcf,'PaperUnits','inches','PaperSize',[5 4],'PaperPosition',[0 0 5 4],'PaperPositionMode','Manual');
        print([folder_Plots,'zsRe_fQ.png'],'-dpng');
        print([folder_Plots,'zsRe_fQ.tif'],'-dtiff','-r600');

        %set up subplot
        clf;
        subplot(1,2,1);
        hold on;
       
        %plot z0 average values
        for i = 1:N_Sites
            plot(0.01,exp(z0Re_ln_fQ_bin_avg_all(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
        end

        %plot z0 SE values
        for i = 1:N_Sites
            plot([0.01 0.01],exp(z0Re_ln_fQ_bin_avg_all(i)+z0Re_ln_fQ_bin_SE_all(i)*[-1 1]),'Color',PlotColors_Site{i});
        end
    end
        
    %plot zs average values
    for i = 1:N_Sites
        plot(fQ_bin_avg_all{i},exp(zsRe_ln_fQ_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
    end

    %plot zs SE values
    for i = 1:N_Sites
        N_fQ_bins = length(fQ_bin_avg_all{i});
        for k = 1:N_fQ_bins
            plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],exp(zsRe_ln_fQ_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
            plot(fQ_bin_avg_all{i}(k)*[1 1],exp(zsRe_ln_fQ_bin_avg_all{i}(k)+zsRe_ln_fQ_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
        end
    end

    %plot fit line
    for i = 1:N_Sites
        plot([0 1],exp(z0Re_ln_fQ_fit(i)+dlnzsRe_dfQ_fit(i)*[0 1]),'Color',PlotColors_Site{i});
        plot([0 0],exp(z0Re_ln_fQ_fit(i)+sigma_z0Re_ln_fQ_fit(i)*[-1 1]),'Color',PlotColors_Site{i},'LineWidth',3);
    end

    %annotate plot
    legend(SiteNames,'Location','SouthEast');
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    set(gca,'YScale','log');
end
ylim([1e-7 1e-2]);
title('Reynolds stress derivation');

%% Log values
subplot(1,2,2); 
hold on;

%plot z0 average values
for i = 1:N_Sites
    plot(0.01,exp(z0Log_ln_fQ_bin_avg_all(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot z0 SE values
for i = 1:N_Sites
    plot([0.01 0.01],exp(z0Log_ln_fQ_bin_avg_all(i)+z0Log_ln_fQ_bin_SE_all(i)*[-1 1]),'Color',PlotColors_Site{i});
end

%plot zs average values
for i = 1:N_Sites
    plot(fQ_bin_avg_all{i},exp(zsLog_ln_fQ_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot zs SE values
for i = 1:N_Sites
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],exp(zsLog_ln_fQ_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
        plot(fQ_bin_avg_all{i}(k)*[1 1],exp(zsLog_ln_fQ_bin_avg_all{i}(k)+zsLog_ln_fQ_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
    end
end

%plot fit line
for i = 1:N_Sites
    plot([0 1],exp(z0Log_ln_fQ_fit(i)+dlnzsLog_dfQ_fit(i)*[0 1]),'Color',PlotColors_Site{i});
    plot([0 0],exp(z0Log_ln_fQ_fit(i)+sigma_z0Log_ln_fQ_fit(i)*[-1 1]),'Color',PlotColors_Site{i},'LineWidth',3);
end

%annotate plot
ylim([1e-7 1e-2]);
legend(SiteNames,'Location','SouthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
title('Log law derivation');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'YScale','log');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
print([folder_Plots,'zs_fQ.png'],'-dpng');
%print([folder_Plots,'zs_fQ.tif'],'-dpng','-dtiff');

%% plot zs versus tau
figure(2); clf; %initialize plot

%% Reynolds values
subplot(1,2,1);
hold on;

%plot zs average values
for i = 1:N_Sites
    plot(tau_bin_avg_all{i},exp(zsRe_ln_tau_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot zs SE values
for i = 1:N_Sites
    N_tau_bins = length(tau_bin_avg_all{i});
    for k = 1:N_tau_bins
        plot(tau_bin_avg_all{i}(k)+tau_bin_SE_all{i}(k)*[-1 1],exp(zsRe_ln_tau_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
        plot(tau_bin_avg_all{i}(k)*[1 1],exp(zsRe_ln_tau_bin_avg_all{i}(k)+zsRe_ln_tau_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
    end
end

%annotate plot
ylim([1e-7 1e-2]);
legend(SiteNames,'Location','SouthEast');
xlabel('shear stress, $$\tau$$ (Pa)','interpreter','latex');
ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
title('Reynolds stress derivation');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'YScale','log');

%% Log values
subplot(1,2,2); 
hold on;

%plot zs average values
for i = 1:N_Sites
    plot(tau_bin_avg_all{i},exp(zsLog_ln_tau_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot zs SE values
for i = 1:N_Sites
    N_tau_bins = length(tau_bin_avg_all{i});
    for k = 1:N_tau_bins
        plot(tau_bin_avg_all{i}(k)+tau_bin_SE_all{i}(k)*[-1 1],exp(zsLog_ln_tau_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
        plot(tau_bin_avg_all{i}(k)*[1 1],exp(zsLog_ln_tau_bin_avg_all{i}(k)+zsLog_ln_tau_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
    end
end

%annotate plot
ylim([1e-7 1e-2]);
legend(SiteNames,'Location','SouthEast');
xlabel('shear stress, $$\tau$$ (Pa)','interpreter','latex');
ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
title('Log law derivation');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'YScale','log');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
print([folder_Plots,'zs_tau.png'],'-dpng');
%print([folder_Plots,'zs_tau.tif'],'-dtiff','-r600');

%% plot zs versus Q
figure(3); clf; %initialize plot

%% Reynolds values - first as separate plot, then joint plot
for n = 1:2
    if n==1
        hold on;
    elseif n==2
        %print plot
        set(gca, 'LooseInset', get(gca,'TightInset'));
        set(gcf,'PaperUnits','inches','PaperSize',[5 4],'PaperPosition',[0 0 5 4],'PaperPositionMode','Manual');
        print([folder_Plots,'zsRe_Q.png'],'-dpng');
        print([folder_Plots,'zsRe_Q.tif'],'-dtiff','-r600');

        %set up subplot
        clf;
        subplot(1,2,1);
        hold on;

        %plot z0 average values
        for i = 1:N_Sites
            plot(0.01,exp(z0Re_ln_Q_bin_avg_all(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
        end

        %plot z0 SE values
        for i = 1:N_Sites
            plot([0.01 0.01],exp(z0Re_ln_Q_bin_avg_all(i)+z0Re_ln_Q_bin_SE_all(i)*[-1 1]),'Color',PlotColors_Site{i});
        end
    end
        
    %plot zs average values
    for i = 1:N_Sites
        plot(Q_bin_avg_all{i},exp(zsRe_ln_Q_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
    end

    %plot zs SE values
    for i = 1:N_Sites
        N_Q_bins = length(Q_bin_avg_all{i});
        for k = 1:N_Q_bins
            plot(Q_bin_avg_all{i}(k)+Q_bin_SE_all{i}(k)*[-1 1],exp(zsRe_ln_Q_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
            plot(Q_bin_avg_all{i}(k)*[1 1],exp(zsRe_ln_Q_bin_avg_all{i}(k)+zsRe_ln_Q_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
        end
    end

    %plot fit line
    for i = 1:N_Sites
        plot([0 Q_max_fit],exp(z0Re_ln_Q_fit(i)+dlnzsRe_dQ_fit(i)*[0 Q_max_fit]),'Color',PlotColors_Site{i});
        plot([0 0],exp(z0Re_ln_Q_fit(i)+sigma_z0Re_ln_Q_fit(i)*[-1 1]),'Color',PlotColors_Site{i},'LineWidth',3);
    end

    %annotate plot
    legend(SiteNames,'Location','SouthEast');
    xlabel('saltation flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)','interpreter','latex');
    ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    set(gca,'YScale','log');
end
title('Reynolds stress derivation');
ylim([1e-7 1e-2]);

%% Log values
subplot(1,2,2); 
hold on;

%plot z0 average values
for i = 1:N_Sites
    plot(0.01,exp(z0Log_ln_Q_bin_avg_all(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot z0 SE values
for i = 1:N_Sites
    plot([0.01 0.01],exp(z0Log_ln_Q_bin_avg_all(i)+z0Log_ln_Q_bin_SE_all(i)*[-1 1]),'Color',PlotColors_Site{i});
end

%plot zs average values
for i = 1:N_Sites
    plot(Q_bin_avg_all{i},exp(zsLog_ln_Q_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot zs SE values
for i = 1:N_Sites
    N_Q_bins = length(Q_bin_avg_all{i});
    for k = 1:N_Q_bins
        plot(Q_bin_avg_all{i}(k)+Q_bin_SE_all{i}(k)*[-1 1],exp(zsLog_ln_Q_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
        plot(Q_bin_avg_all{i}(k)*[1 1],exp(zsLog_ln_Q_bin_avg_all{i}(k)+zsLog_ln_Q_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
    end
end

%plot fit line
for i = 1:N_Sites
    plot([0 Q_max_fit],exp(z0Log_ln_Q_fit(i)+dlnzsLog_dQ_fit(i)*[0 Q_max_fit]),'Color',PlotColors_Site{i});
    plot([0 0],exp(z0Log_ln_Q_fit(i)+sigma_z0Log_ln_Q_fit(i)*[-1 1]),'Color',PlotColors_Site{i},'LineWidth',3);
end

%annotate plot
ylim([1e-7 1e-2]);
legend(SiteNames,'Location','SouthEast');
xlabel('saltation flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)','interpreter','latex');
ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
title('Log law derivation');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'YScale','log');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
print([folder_Plots,'zs_Q.png'],'-dpng');
%print([folder_Plots,'zs_Q.tif'],'-dtiff','-r600');


%% combined roughness plot for paper - zs versus Q and zs versus fQ
figure(4); clf; %initialize plot

%set up subplot - zs versus Q
subplot('Position',[0.08 0.1 0.43 0.88]); hold on;

%plot z0 average values
for i = 1:N_Sites
    plot(0.01,exp(z0Re_ln_Q_bin_avg_all(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot z0 SE values
for i = 1:N_Sites
    plot([0.01 0.01],exp(z0Re_ln_Q_bin_avg_all(i)+z0Re_ln_Q_bin_SE_all(i)*[-1 1]),'Color',PlotColors_Site{i});
end

%plot zs average values
for i = 1:N_Sites
    plot(Q_bin_avg_all{i},exp(zsRe_ln_Q_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot zs SE values
for i = 1:N_Sites
    N_Q_bins = length(Q_bin_avg_all{i});
    for k = 1:N_Q_bins
        plot(Q_bin_avg_all{i}(k)+Q_bin_SE_all{i}(k)*[-1 1],exp(zsRe_ln_Q_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
        plot(Q_bin_avg_all{i}(k)*[1 1],exp(zsRe_ln_Q_bin_avg_all{i}(k)+zsRe_ln_Q_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
    end
end

%plot fit line
for i = 1:N_Sites
    plot([0 Q_max_fit],exp(z0Re_ln_Q_fit(i)+dlnzsRe_dQ_fit(i)*[0 Q_max_fit]),'Color',PlotColors_Site{i});
    plot([0 0],exp(z0Re_ln_Q_fit(i)+sigma_z0Re_ln_Q_fit(i)*[-1 1]),'Color',PlotColors_Site{i},'LineWidth',3);
end

%annotate and format
legend(SiteNames,'Location','SouthEast');
xlabel('saltation flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)','interpreter','latex');
ylabel('effective roughness height, $$z_{s}$$ (m)','interpreter','latex');
text(2.5,8e-3,'(a)');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'YScale','log');
ylim([1e-5 1e-2]);

%set up subplot - zs versus fQ
subplot('Position',[0.56 0.1 0.43 0.88]); hold on;

%plot z0 average values
for i = 1:N_Sites
    plot(0.01,exp(z0Re_ln_fQ_bin_avg_all(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot z0 SE values
for i = 1:N_Sites
    plot([0.01 0.01],exp(z0Re_ln_fQ_bin_avg_all(i)+z0Re_ln_fQ_bin_SE_all(i)*[-1 1]),'Color',PlotColors_Site{i});
end

        
%plot zs average values
for i = 1:N_Sites
    plot(fQ_bin_avg_all{i},exp(zsRe_ln_fQ_bin_avg_all{i}),PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot zs SE values
for i = 1:N_Sites
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],exp(zsRe_ln_fQ_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors_Site{i});
        plot(fQ_bin_avg_all{i}(k)*[1 1],exp(zsRe_ln_fQ_bin_avg_all{i}(k)+zsRe_ln_fQ_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors_Site{i});
    end
end

%plot fit line
for i = 1:N_Sites
    plot([0 1],exp(z0Re_ln_fQ_fit(i)+dlnzsRe_dfQ_fit(i)*[0 1]),'Color',PlotColors_Site{i});
    plot([0 0],exp(z0Re_ln_fQ_fit(i)+sigma_z0Re_ln_fQ_fit(i)*[-1 1]),'Color',PlotColors_Site{i},'LineWidth',3);
end

%annotate and format plot
legend(SiteNames,'Location','SouthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
text(0.05,8e-3,'(b)');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','YTickLabel',[]);
set(gca,'YScale','log');
ylim([1e-5 1e-2]);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
print([folder_Plots,'zs_Q_z0_fQ.png'],'-dpng');
print([folder_Plots,'zs_Q_z0_fQ.tif'],'-dtiff','-r600');