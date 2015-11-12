%% ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP
% note: intervals 779:782 for Oceano look weird (June 1 14:28-15:13)

%initialize
clearvars;
close all;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_AnalysisData = '../AnalysisData/'; %folder for storing data output
folder_Plots = '../PlotOutput/FluxFrequency/'; %folder for plots

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows_all'));

%get info about number of sites
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%set outlier windows
ind_outlier = {[];[];[]};
%ind_outlier = {[];[];[779:782]};

%% PARAMETERS
kappa = 0.39; %von Karman
rho = 1.23; %air density, kg/m^3
ust_th = [0.35, 0.28, 0.28]; %assumed u* threshold (m/s) for each site
tau_th = rho*ust_th.^2; %assumed tau threshold (Pa) for each site
z0f = [9e-5, 1e-4, 2e-4]; %aerodynamic roughness length (m) at threshold
%z0f = [1e-5, 1e-5, 1e-5]; %aerodynamic roughness length (m) at threshold
%z0f = [3.7e-5, 3.7e-5, 8.5e-5]; %aerodynamic roughness length (m) at threshold

%% CALCULATIONS
%alternate calculation for z0 based only on lowest anemometer (standard is based on lowermost three anemometers)
%calculation of ust_th1 and tau_th1 (thresholds based on TFEM method for fQ1)
z0alt = cell(N_Sites,1);
ust_th1 = cell(N_Sites,1);
tau_th1 = cell(N_Sites,1);
for i=1:N_Sites
    z_lowest_matrix = cell2mat(z_lowest_all{i}')';
    z_lowest_list = z_lowest_matrix(:,1);
    z0alt{i} = z_lowest_list.*exp(-kappa*u_bar_cal_all{i}./ustRe_cal_all{i});
    ust_th1{i} = kappa*uth1_cal_all{i}./log(z_lowest_list./z0f(i)); %compute u* threshold, using constant z0
    %ust_th1{i} = kappa*uth1_cal_all{i}./log(z_lowest_list./z0alt{i}); %compute u* threshold, using observed z0
    ust_th1{i}(fQ1_all{i}==0) = NaN; %if fQ = 0 (no transport) then TFEM method is bunk and the values should be NaN
    tau_th1{i} = rho*ust_th1{i}.^2; %compute tau threshold
    %z0_bar = exp(mean(log(z0_cal_all{i}(Q_all{i}==0&~isnan(z0_cal_all{i})))));
    %z0alt_bar = exp(mean(log(z0alt{i}(Q_all{i}==0&~isnan(z0_cal_all{i})))));
end

%compute eta based on z0
etaz0_all = cell(N_Sites,1);
for i=1:N_Sites
    %etaz0_all{i} = real(1-log(ze_all{i}./z0_cal_all{i}).^2./log(ze_all{i}./z0f(i)).^2); %using original z0
    etaz0_all{i} = real(1-log(ze_all{i}./z0alt{i}).^2./log(ze_all{i}./z0f(i)).^2); %using alt z0
end

% %% COMPUTE USTEX AND TAUEX BASED ON TFEM THRESHOLD
% ustex_all = cell(N_Sites,1);
% tauex_all = cell(N_Sites,1);
% tauratio_all = cell(N_Sites,1);
% for i=1:N_Sites
%     ustex_all{i} = ustRe_cal_all{i}-ust_th1{i};
%     tauex_all{i} = tauRe_cal_all{i}-tau_th1{i};
%     tauratio_all{i} = tauex_all{i}./tauRe_cal_all{i};
% end

%COMPUTE USTEX AND TAUEX BASED ON THRESHOLD BY SITE
ustex_all = cell(N_Sites,1);
tauex_all = cell(N_Sites,1);
tauratio_all = cell(N_Sites,1);
for i=1:N_Sites
    ustex_all{i} = ustRe_cal_all{i}-ust_th(i);
    tauex_all{i} = tauRe_cal_all{i}-tau_th(i);
    tauratio_all{i} = tauex_all{i}./tauRe_cal_all{i};
end

% %COMPUTE USTEX AND TAUEX BASED ON THRESHOLDS BY DATE
% ustex_all = cell(N_Sites,1);
% tauex_all = cell(N_Sites,1);
% tauratio_all = cell(N_Sites,1);
% for i=1:N_Sites
%     ustex_all{i} = zeros(size(date_all{i}));
%     tauex_all{i} = zeros(size(date_all{i}));
%     tauratio_all{i} = zeros(size(date_all{i}));
%     
%     %apply empirical threshold - Jericoacoara - all dates
%     if i == 1;
%         ustex_all{i} = ustRe_cal_all{i}-0.3474;
%         tauex_all{i} = tauRe_cal_all{i}-0.14841;
%         tauratio_all{i} = tauex_all{i}/0.14841;
%     %apply empirical threshold - Rancho Guadalupe - all dates
%     elseif i == 2;
%         ustex_all{i} = ustRe_cal_all{i}-0.2839;
%         tauex_all{i} = tauRe_cal_all{i}-0.099111;
%         tauratio_all{i} = tauex_all{i}/0.099111;
%     elseif i == 3;
%         %apply empirical threshold - Oceano - May 15-19
%         date_ind = intersect(find(date_all{i}>=datetime(2015,5,15)),find(date_all{i}<=datetime(2015,5,19)));
%         ustex_all{i}(date_ind) = ustRe_cal_all{i}(date_ind)-0.2794;
%         tauex_all{i}(date_ind) = tauRe_cal_all{i}(date_ind)-0.096036;
%         tauratio_all{i}(date_ind) = tauex_all{i}(date_ind)/0.096036;
%         %apply empirical threshold - Oceano - May 23-31
%         date_ind = intersect(find(date_all{i}>=datetime(2015,5,23)),find(date_all{i}<=datetime(2015,5,31)));
%         ustex_all{i}(date_ind) = ustRe_cal_all{i}(date_ind)-0.2903;
%         tauex_all{i}(date_ind) = tauRe_cal_all{i}(date_ind)-0.10364;
%         tauratio_all{i}(date_ind) = tauex_all{i}(date_ind)/0.10364;
%         %apply empirical threshold - Oceano - June 1-4
%         date_ind = intersect(find(date_all{i}>=datetime(2015,6,1)),find(date_all{i}<=datetime(2015,6,4)));
%         ustex_all{i}(date_ind) = ustRe_cal_all{i}(date_ind)-0.2792;
%         tauex_all{i}(date_ind) = tauRe_cal_all{i}(date_ind)-0.095865;
%         tauratio_all{i}(date_ind) = tauex_all{i}(date_ind)/0.095865;
%     end
% end

%% CREATE BINS
%set frequency bins
fQ_bins_min = [0 0.3 0.6 0.9];
fQ_bins_max = [0.3 0.6 0.9 1];
fQ_bins_mid = mean([fQ_bins_min; fQ_bins_max]);
fQ_bins_legend_all = {'f_{Q} = 0-0.3','f_{Q} = 0.3-0.6', 'f_{Q} = 0.6-0.9', 'f_{Q} = 0.9-1'};
N_fQ_bins = length(fQ_bins_mid);
Markers_fQ_bins = {'rx','bv','g^','ko'};

%create u* bins
ust_bins_min = 0.025:.025:0.55;
ust_bins_max = 0.05:.025:0.575;
ust_bins_mid = mean([ust_bins_min;ust_bins_max]);
N_ust_bins = length(ust_bins_mid);

%tau_bins
% tau_bins_min = 0:.025:0.375;
% tau_bins_max = 0.025:.025:0.4;
tau_bins_min = 0:.0125:0.3875;
tau_bins_max = 0.0125:.0125:0.4;
tau_bins_mid = mean([tau_bins_min;tau_bins_max]);
N_tau_bins = length(tau_bins_mid);

%tauex_bins
% tauex_bins_min = -0.15:0.025:0.275;
% tauex_bins_max = -0.125:0.025:0.3;
tauex_bins_min = -0.1:0.0125:0.2375;
tauex_bins_max = -0.0875:0.0125:0.25;
tauex_bins_mid = mean([tauex_bins_min;tauex_bins_max]);
N_tauex_bins = length(tauex_bins_mid);

%etatau_bins
etatau_bins_min = 0:0.01:0.22;
etatau_bins_max = 0.01:0.01:0.23;
etatau_bins_mid = mean([etatau_bins_min;etatau_bins_max]);
N_etatau_bins = length(etatau_bins_mid);

%ustth_bins
ustth_bins_min = 0.21:0.01:0.38;
ustth_bins_max = 0.22:0.01:0.39;
ustth_bins_mid = mean([ustth_bins_min;ustth_bins_max]);
N_ustth_bins = length(ustth_bins_mid);

%% FOR EACH SITE, PERFORM BINNING BY fQ
fQ_bins_legend = cell(N_Sites,1); %get legend values for each site
ust_fQ_bin_values = cell(N_Sites,1); %u* into frequency bins
tau_fQ_bin_values = cell(N_Sites,1); %tau into frequency bins
tauex_fQ_bin_values = cell(N_Sites,1); %tauex into frequency bins
Q_fQ_bin_values = cell(N_Sites,1); %fluxes into frequency bins
eta_fQ_bin_values = cell(N_Sites,1); %etas into frequency bins
z0_fQ_bin_values = cell(N_Sites,1); %z0s into frequency bins
z0alt_fQ_bin_values = cell(N_Sites,1); %alternate calcs for z0 into frequency bins
etaz0_fQ_bin_values = cell(N_Sites,1); %etas based on z0 into frequency bins
ubarustd_fQ_bin_values = cell(N_Sites,1); %ubar/ustd into frequency bins

%apply values to bins
for i=1:N_Sites
    
    %initialize bins for each site
    fQ_bins_legend{i} = cell(N_Sites,1);
    ust_fQ_bin_values{i} = cell(N_fQ_bins,1);
    tau_fQ_bin_values{i} = cell(N_fQ_bins,1);
    tauex_fQ_bin_values{i} = cell(N_fQ_bins,1);
    Q_fQ_bin_values{i} = cell(N_fQ_bins,1);
    eta_fQ_bin_values{i} = cell(N_fQ_bins,1);
    z0_fQ_bin_values{i} = cell(N_fQ_bins,1);
    z0alt_fQ_bin_values{i} = cell(N_fQ_bins,1);
    etaz0_fQ_bin_values{i} = cell(N_fQ_bins,1);
    ubarustd_fQ_bin_values{i} = cell(N_fQ_bins,1);
    
    %go through frequency bins, get values
    for j=1:N_fQ_bins
        bin_ind = find(fQ1_all{i}>=fQ_bins_min(j)&fQ1_all{i}<=fQ_bins_max(j)); %use 1 second frequencies for binning
        %bin_ind = find(fQ_all{i}(:,1)>=fQ_bins_min(j)&fQ_all{i}(:,1)<=fQ_bins_max(j)); %use shortest time interval for binning
        bin_ind = setdiff(bin_ind,ind_outlier{i}); %remove outliers specified above
        fQ_bins_legend{i}{j} = strcat(fQ_bins_legend_all{j},'; N = ',int2str(length(bin_ind)));
        if ~isempty(bin_ind)
            ust_fQ_bin_values{i}{j} = ustRe_cal_all{i}(bin_ind);
            tau_fQ_bin_values{i}{j} = tauRe_cal_all{i}(bin_ind);
            tauex_fQ_bin_values{i}{j} = tauex_all{i}(bin_ind);
            Q_fQ_bin_values{i}{j} = Q_all{i}(bin_ind);
            eta_fQ_bin_values{i}{j} = eta_cal_all{i}(bin_ind);
            z0_fQ_bin_values{i}{j} = z0_cal_all{i}(bin_ind);
            z0alt_fQ_bin_values{i}{j} = z0alt{i}(bin_ind);
            etaz0_fQ_bin_values{i}{j} = etaz0_all{i}(bin_ind);
            ubarustd_fQ_bin_values{i}{j} = u_bar_cal_all{i}(bin_ind)./u_std_cal_all{i}(bin_ind);
        end
    end
end

%% FOR EACH SITE, PERFORM BINNING BY ETA
eta_bins_min = 0:.025:0.475;
eta_bins_max = 0.025:.025:0.5;
eta_bins_mid = mean([eta_bins_min;eta_bins_max]);
N_eta_bins = length(eta_bins_mid);

%separate fQ's into eta bins
fQ_eta_bin_values = cell(N_Sites,1); %fQ's into eta bins
fQ_eta_bin_avg = cell(N_Sites,1); %get average fQ for eta bin
fQ_eta_bin_std = cell(N_Sites,1); %get std fQ for eta bin
fQ_eta_bin_SE = cell(N_Sites,1); %get SE fQ for eta bin

for i = 1:N_Sites
    %initialize fQ's into eta bins
    fQ_eta_bin_values{i} = cell(N_eta_bins,1); 
    fQ_eta_bin_avg{i} = zeros(N_eta_bins,1)*NaN;
    fQ_eta_bin_std{i} = zeros(N_eta_bins,1)*NaN;
    fQ_eta_bin_SE{i} = zeros(N_eta_bins,1)*NaN;
    
    %get fQs for eta bin
    for j = 1:N_eta_bins
        bin_ind = find(eta_cal_all{i}>=eta_bins_min(j)&eta_cal_all{i}<=eta_bins_max(j));
        fQ_eta_bin_values{i}{j} = fQ1_all{i}(bin_ind);
        fQ_eta_bin_values{i}{j} = fQ_eta_bin_values{i}{j}(~isnan(fQ_eta_bin_values{i}{j}));
        fQ_eta_bin_avg{i}(j) = mean(fQ_eta_bin_values{i}{j});
        fQ_eta_bin_std{i}(j) = std(fQ_eta_bin_values{i}{j});
        fQ_eta_bin_SE{i}(j) = std(fQ_eta_bin_values{i}{j})./sqrt(length(~isnan(fQ_eta_bin_values{i}{j})));
    end
end

%% FOR EACH SITE, PERFORM BINNING BY TAUEX
%separate tauratios into tau_ex bins
tauratio_tauex_bin_values = cell(N_Sites,1); %tauratios into tau_ex bins
tauratio_tauex_bin_avg = cell(N_Sites,1); %get average tauratio for tau_ex bins
tauratio_tauex_bin_std = cell(N_Sites,1); %get std tauratio for tau_ex bins
tauratio_tauex_bin_SE = cell(N_Sites,1); %get SE tauratio for tau_ex bins

for i = 1:N_Sites
    %initialize tauratio's into tauex bins
    tauratio_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    tauratio_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    tauratio_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    tauratio_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %get tauratios for tauex bins
    for j = 1:N_tauex_bins
        bin_ind = find(tauex_all{i}>=tauex_bins_min(j)&tauex_all{i}<=tauex_bins_max(j));
        tauratio_tauex_bin_values{i}{j} = tauratio_all{i}(bin_ind);
        tauratio_tauex_bin_values{i}{j} = tauratio_tauex_bin_values{i}{j}(~isnan(tauratio_tauex_bin_values{i}{j}));
        tauratio_tauex_bin_avg{i}(j) = mean(tauratio_tauex_bin_values{i}{j});
        tauratio_tauex_bin_std{i}(j) = std(tauratio_tauex_bin_values{i}{j});
        tauratio_tauex_bin_SE{i}(j) = std(tauratio_tauex_bin_values{i}{j})./sqrt(length(~isnan(tauratio_tauex_bin_values{i}{j})));
    end
end

%% FOR EACH SITE, PERFORM BINNING BY ALTERNATE fQ BINS
%(more bins, for direct comparisons to frequency)
fQ_bins_alt_min = 0:0.05:0.95;
fQ_bins_alt_max = 0.05:0.05:1;
fQ_bins_alt_mid = mean([fQ_bins_alt_min; fQ_bins_alt_max]);
N_fQ_bins_alt = length(fQ_bins_alt_mid);

%separate u*th based on TFEM into frequency bins
ustth_fQ_bin_values = cell(N_Sites,1);
ustth_fQ_bin_avg = cell(N_Sites,1);
ustth_fQ_bin_std = cell(N_Sites,1);
ustth_fQ_bin_SE = cell(N_Sites,1);
tauth_fQ_bin_values = cell(N_Sites,1);
tauth_fQ_bin_avg = cell(N_Sites,1);
tauth_fQ_bin_std = cell(N_Sites,1);
tauth_fQ_bin_SE = cell(N_Sites,1);

for i = 1:N_Sites
    ustth_fQ_bin_values{i} = cell(N_fQ_bins_alt,1);
    ustth_fQ_bin_avg{i} = zeros(N_fQ_bins_alt,1)*NaN;
    ustth_fQ_bin_std{i} = zeros(N_fQ_bins_alt,1)*NaN;
    ustth_fQ_bin_SE{i} = zeros(N_fQ_bins_alt,1)*NaN;
    tauth_fQ_bin_values{i} = cell(N_fQ_bins_alt,1);
    tauth_fQ_bin_avg{i} = zeros(N_fQ_bins_alt,1)*NaN;
    tauth_fQ_bin_std{i} = zeros(N_fQ_bins_alt,1)*NaN;
    tauth_fQ_bin_SE{i} = zeros(N_fQ_bins_alt,1)*NaN;
    
    %get u*th's / tauth's for fQ bins
    for j = 1:N_fQ_bins_alt
        bin_ind = find(fQ1_all{i}>=fQ_bins_alt_min(j)&fQ1_all{i}<=fQ_bins_alt_max(j)); %use 1 second frequencies for binning
        ustth_fQ_bin_values{i}{j} = ust_th1{i}(bin_ind);
        ustth_fQ_bin_values{i}{j} = ustth_fQ_bin_values{i}{j}(~isnan(ustth_fQ_bin_values{i}{j}));
        ustth_fQ_bin_avg{i}(j) = mean(ustth_fQ_bin_values{i}{j});
        ustth_fQ_bin_std{i}(j) = std(ustth_fQ_bin_values{i}{j});
        ustth_fQ_bin_SE{i}(j) = std(ustth_fQ_bin_values{i}{j})/sqrt(length(~isnan(ustth_fQ_bin_values{i}{j})));
        tauth_fQ_bin_values{i}{j} = tau_th1{i}(bin_ind);
        tauth_fQ_bin_values{i}{j} = tauth_fQ_bin_values{i}{j}(~isnan(tauth_fQ_bin_values{i}{j}));
        tauth_fQ_bin_avg{i}(j) = mean(tauth_fQ_bin_values{i}{j});
        tauth_fQ_bin_std{i}(j) = std(tauth_fQ_bin_values{i}{j});
        tauth_fQ_bin_SE{i}(j) = std(tauth_fQ_bin_values{i}{j})/sqrt(length(~isnan(tauth_fQ_bin_values{i}{j})));
    end 
end

%% PERFORM SUB-BINNING WITHIN FREQUENCY BINS, CREATE PLOTS FOR EACH SITE
for i = 1:N_Sites
    
    %get indices of fQ bins that are not empty
    ind_fQ_bins_full = find(1-cellfun(@isempty,Q_fQ_bin_values{i}))';
    
    %fluxes and other values into u* bins within fQ bins
    Q_fQ_ust_bin_avg = cell(N_ust_bins,1);
    Q_fQ_ust_bin_std = cell(N_ust_bins,1);
    Q_fQ_ust_bin_SE = cell(N_ust_bins,1);
    eta_fQ_ust_bin_avg = cell(N_ust_bins,1);
    eta_fQ_ust_bin_std = cell(N_ust_bins,1);
    eta_fQ_ust_bin_SE = cell(N_ust_bins,1);
    z0_fQ_ust_bin_avg = cell(N_ust_bins,1);
    z0_fQ_ust_bin_std = cell(N_ust_bins,1);
    z0_fQ_ust_bin_SE = cell(N_ust_bins,1);
    z0alt_fQ_ust_bin_avg = cell(N_ust_bins,1);
    z0alt_fQ_ust_bin_std = cell(N_ust_bins,1);
    z0alt_fQ_ust_bin_SE = cell(N_ust_bins,1);
    etaz0_fQ_ust_bin_avg = cell(N_ust_bins,1);
    etaz0_fQ_ust_bin_std = cell(N_ust_bins,1);
    etaz0_fQ_ust_bin_SE = cell(N_ust_bins,1);
    ubarustd_fQ_ust_bin_avg = cell(N_ust_bins,1);
    ubarustd_fQ_ust_bin_std = cell(N_ust_bins,1);
    ubarustd_fQ_ust_bin_SE = cell(N_ust_bins,1);
    
    %fluxes into tau bins within fQ bins
    Q_fQ_tau_bin_avg = cell(N_tau_bins,1);
    Q_fQ_tau_bin_std = cell(N_tau_bins,1);
    Q_fQ_tau_bin_SE = cell(N_tau_bins,1);
    
    %fluxes and other values into tauex bins within fQ bins
    Q_fQ_tauex_bin_avg = cell(N_tauex_bins,1);
    Q_fQ_tauex_bin_std = cell(N_tauex_bins,1);
    Q_fQ_tauex_bin_SE = cell(N_tauex_bins,1);
    eta_fQ_tauex_bin_avg = cell(N_tauex_bins,1);
    eta_fQ_tauex_bin_std = cell(N_tauex_bins,1);
    eta_fQ_tauex_bin_SE = cell(N_tauex_bins,1);
    z0_fQ_tauex_bin_avg = cell(N_tauex_bins,1);
    z0_fQ_tauex_bin_std = cell(N_tauex_bins,1);
    z0_fQ_tauex_bin_SE = cell(N_tauex_bins,1);
    z0alt_fQ_tauex_bin_avg = cell(N_tauex_bins,1);
    z0alt_fQ_tauex_bin_std = cell(N_tauex_bins,1);
    z0alt_fQ_tauex_bin_SE = cell(N_tauex_bins,1);
    etaz0_fQ_tauex_bin_avg = cell(N_tauex_bins,1);
    etaz0_fQ_tauex_bin_std = cell(N_tauex_bins,1);
    etaz0_fQ_tauex_bin_SE = cell(N_tauex_bins,1);
    
    %fluxes into eta*tau bins within fQ bins
    Q_fQ_etatau_bin_avg = cell(N_etatau_bins,1);
    Q_fQ_etatau_bin_std = cell(N_etatau_bins,1);
    Q_fQ_etatau_bin_SE = cell(N_etatau_bins,1);
    Q_fQ_etaz0tau_bin_avg = cell(N_etatau_bins,1);
    Q_fQ_etaz0tau_bin_std = cell(N_etatau_bins,1);
    Q_fQ_etaz0tau_bin_SE = cell(N_etatau_bins,1);

    %go through each fQ bin
    for j = 1:N_fQ_bins
        
        %initialize u* bins for each fQ bin
        Q_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        Q_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        Q_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        eta_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        eta_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        eta_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        z0_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        z0_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        z0_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        z0alt_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        z0alt_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        z0alt_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        etaz0_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        etaz0_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        etaz0_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        ubarustd_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        ubarustd_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        ubarustd_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        
        %initialize tau bins for each fQ bin
        Q_fQ_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
        Q_fQ_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
        Q_fQ_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
  
        %initialize tau_ex bins for each fQ bin
        Q_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
        Q_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
        Q_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
        eta_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
        eta_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
        eta_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
        z0_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
        z0_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
        z0_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
        z0alt_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
        z0alt_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
        z0alt_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
        etaz0_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
        etaz0_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
        etaz0_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
        
        %initialize etatau bins for each fQ bin
        Q_fQ_etatau_bin_avg{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etatau_bin_std{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etatau_bin_SE{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etaz0tau_bin_avg{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etaz0tau_bin_std{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etaz0tau_bin_SE{j} = zeros(N_etatau_bins,1)*NaN;
        
        %go through each u* bin
        for k=1:N_ust_bins
            bin_ind = find(ust_fQ_bin_values{i}{j}>=ust_bins_min(k)&ust_fQ_bin_values{i}{j}<=ust_bins_max(k));
            if ~isempty(bin_ind)
                Q_fQ_ust_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
                Q_fQ_ust_bin_values = Q_fQ_ust_bin_values(~isnan(Q_fQ_ust_bin_values));
                Q_fQ_ust_bin_avg{j}(k) = mean(Q_fQ_ust_bin_values);
                Q_fQ_ust_bin_SE{j}(k) = std(Q_fQ_ust_bin_values)/sqrt(length(~isnan(Q_fQ_ust_bin_values)));
                Q_fQ_ust_bin_std{j}(k) = std(Q_fQ_ust_bin_values);
                
                eta_fQ_ust_bin_values = eta_fQ_bin_values{i}{j}(bin_ind);
                eta_fQ_ust_bin_values = eta_fQ_ust_bin_values(~isnan(eta_fQ_ust_bin_values));
                eta_fQ_ust_bin_avg{j}(k) = mean(eta_fQ_ust_bin_values);
                eta_fQ_ust_bin_SE{j}(k) = std(eta_fQ_ust_bin_values)/sqrt(length(~isnan(eta_fQ_ust_bin_values)));
                eta_fQ_ust_bin_std{j}(k) = std(eta_fQ_ust_bin_values);
                
                z0_fQ_ust_bin_values = z0_fQ_bin_values{i}{j}(bin_ind);
                z0_fQ_ust_bin_values = z0_fQ_ust_bin_values(~isnan(z0_fQ_ust_bin_values));
                z0_fQ_ust_bin_avg{j}(k) = exp(mean(log(z0_fQ_ust_bin_values))); %get mean in log space
                z0_fQ_ust_bin_SE{j}(k) = std(z0_fQ_ust_bin_values)/sqrt(length(~isnan(z0_fQ_ust_bin_values)));
                z0_fQ_ust_bin_std{j}(k) = std(z0_fQ_ust_bin_values);
                
                z0alt_fQ_ust_bin_values = z0alt_fQ_bin_values{i}{j}(bin_ind);
                z0alt_fQ_ust_bin_values = z0alt_fQ_ust_bin_values(~isnan(z0alt_fQ_ust_bin_values));
                z0alt_fQ_ust_bin_avg{j}(k) = exp(mean(log(z0alt_fQ_ust_bin_values))); %get mean in log space
                z0alt_fQ_ust_bin_SE{j}(k) = std(z0alt_fQ_ust_bin_values)/sqrt(length(~isnan(z0alt_fQ_ust_bin_values)));
                z0alt_fQ_ust_bin_std{j}(k) = std(z0alt_fQ_ust_bin_values);
                
                etaz0_fQ_ust_bin_values = etaz0_fQ_bin_values{i}{j}(bin_ind);
                etaz0_fQ_ust_bin_values = etaz0_fQ_ust_bin_values(~isnan(etaz0_fQ_ust_bin_values));
                etaz0_fQ_ust_bin_avg{j}(k) = mean(etaz0_fQ_ust_bin_values);
                etaz0_fQ_ust_bin_SE{j}(k) = std(etaz0_fQ_ust_bin_values)/sqrt(length(~isnan(etaz0_fQ_ust_bin_values)));
                etaz0_fQ_ust_bin_std{j}(k) = std(etaz0_fQ_ust_bin_values);
                
                ubarustd_fQ_ust_bin_values = ubarustd_fQ_bin_values{i}{j}(bin_ind);
                ubarustd_fQ_ust_bin_values = ubarustd_fQ_ust_bin_values(~isnan(ubarustd_fQ_ust_bin_values));
                ubarustd_fQ_ust_bin_avg{j}(k) = mean(ubarustd_fQ_ust_bin_values);
                ubarustd_fQ_ust_bin_SE{j}(k) = std(ubarustd_fQ_ust_bin_values)/sqrt(length(~isnan(ubarustd_fQ_ust_bin_values)));
                ubarustd_fQ_ust_bin_std{j}(k) = std(ubarustd_fQ_ust_bin_values);
            end
        end
        
        %go through each etatau bin
        for k=1:N_tau_bins
            bin_ind = find(tau_fQ_bin_values{i}{j}>=tau_bins_min(k)&tau_fQ_bin_values{i}{j}<=tau_bins_max(k));
            if ~isempty(bin_ind)
                Q_fQ_tau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
                Q_fQ_tau_bin_values = Q_fQ_tau_bin_values(~isnan(Q_fQ_tau_bin_values));
                Q_fQ_tau_bin_avg{j}(k) = mean(Q_fQ_tau_bin_values);
                Q_fQ_tau_bin_SE{j}(k) = std(Q_fQ_tau_bin_values)/sqrt(length(~isnan(Q_fQ_tau_bin_values)));
                Q_fQ_tau_bin_std{j}(k) = std(Q_fQ_tau_bin_values);
            end
        end
        
        %go through each tau_ex bin
        for k=1:N_tauex_bins
            bin_ind = find(tauex_fQ_bin_values{i}{j}>=tauex_bins_min(k)&tauex_fQ_bin_values{i}{j}<=tauex_bins_max(k));
            if ~isempty(bin_ind)
                Q_fQ_tauex_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
                Q_fQ_tauex_bin_values = Q_fQ_tauex_bin_values(~isnan(Q_fQ_tauex_bin_values));
                Q_fQ_tauex_bin_avg{j}(k) = mean(Q_fQ_tauex_bin_values);
                Q_fQ_tauex_bin_SE{j}(k) = std(Q_fQ_tauex_bin_values)/sqrt(length(~isnan(Q_fQ_tauex_bin_values)));
                Q_fQ_tauex_bin_std{j}(k) = std(Q_fQ_tauex_bin_values);
                
                eta_fQ_tauex_bin_values = eta_fQ_bin_values{i}{j}(bin_ind);
                eta_fQ_tauex_bin_values = eta_fQ_tauex_bin_values(~isnan(eta_fQ_tauex_bin_values));
                eta_fQ_tauex_bin_avg{j}(k) = mean(eta_fQ_tauex_bin_values);
                eta_fQ_tauex_bin_SE{j}(k) = std(eta_fQ_tauex_bin_values)/sqrt(length(~isnan(eta_fQ_tauex_bin_values)));
                eta_fQ_tauex_bin_std{j}(k) = std(eta_fQ_tauex_bin_values);
                
                z0_fQ_tauex_bin_values = z0_fQ_bin_values{i}{j}(bin_ind);
                eta_fQ_tauex_bin_values = eta_fQ_tauex_bin_values(~isnan(eta_fQ_tauex_bin_values));
                z0_fQ_tauex_bin_avg{j}(k) = exp(mean(log(z0_fQ_tauex_bin_values))); %get mean in log space
                z0_fQ_tauex_bin_SE{j}(k) = std(z0_fQ_tauex_bin_values)/sqrt(length(~isnan(z0_fQ_tauex_bin_values)));
                z0_fQ_tauex_bin_std{j}(k) = std(z0_fQ_tauex_bin_values);
                
                z0alt_fQ_tauex_bin_values = z0alt_fQ_bin_values{i}{j}(bin_ind);
                z0alt_fQ_tauex_bin_values = z0alt_fQ_tauex_bin_values(~isnan(z0alt_fQ_tauex_bin_values));
                z0alt_fQ_tauex_bin_avg{j}(k) = exp(mean(log(z0alt_fQ_tauex_bin_values))); %get mean in log space
                z0alt_fQ_tauex_bin_SE{j}(k) = std(z0alt_fQ_tauex_bin_values)/sqrt(length(~isnan(z0alt_fQ_tauex_bin_values)));
                z0alt_fQ_tauex_bin_std{j}(k) = std(z0alt_fQ_tauex_bin_values);
                
                etaz0_fQ_tauex_bin_values = etaz0_fQ_bin_values{i}{j}(bin_ind);
                etaz0_fQ_tauex_bin_values = etaz0_fQ_tauex_bin_values(~isnan(etaz0_fQ_tauex_bin_values));
                etaz0_fQ_tauex_bin_avg{j}(k) = mean(etaz0_fQ_tauex_bin_values);
                etaz0_fQ_tauex_bin_SE{j}(k) = std(etaz0_fQ_tauex_bin_values)/sqrt(length(~isnan(etaz0_fQ_tauex_bin_SE{j}(k))));
                etaz0_fQ_tauex_bin_std{j}(k) = std(etaz0_fQ_tauex_bin_values);

            end
        end
        
        %go through each etatau bin
        for k=1:N_etatau_bins
            %first, definition of eta based on fluctuations in u
            bin_ind = find(eta_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}>=etatau_bins_min(k)&...
                eta_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}<=etatau_bins_max(k));
            if ~isempty(bin_ind)
                Q_fQ_etatau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
                Q_fQ_etatau_bin_values = Q_fQ_etatau_bin_values(~isnan(Q_fQ_etatau_bin_values));
                Q_fQ_etatau_bin_avg{j}(k) = mean(Q_fQ_etatau_bin_values);
                Q_fQ_etatau_bin_SE{j}(k) = std(Q_fQ_etatau_bin_values)/sqrt(length(~isnan(Q_fQ_etatau_bin_values)));
                Q_fQ_etatau_bin_std{j}(k) = std(Q_fQ_etatau_bin_values);
            end
            %second, definition of eta based on fluctuations in z0
            bin_ind = find(etaz0_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}>=etatau_bins_min(k)&...
                etaz0_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}<=etatau_bins_max(k));
            if ~isempty(bin_ind)
                Q_fQ_etaz0tau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
                Q_fQ_etaz0tau_bin_values = Q_fQ_etaz0tau_bin_values(~isnan(Q_fQ_etaz0tau_bin_values));
                Q_fQ_etaz0tau_bin_avg{j}(k) = mean(Q_fQ_etaz0tau_bin_values);
                Q_fQ_etaz0tau_bin_SE{j}(k) = std(Q_fQ_etaz0tau_bin_values)/sqrt(length(~isnan(Q_fQ_etaz0tau_bin_values)));
                Q_fQ_etaz0tau_bin_std{j}(k) = std(Q_fQ_etaz0tau_bin_values);
            end
        end    
    end

    %% PLOTS %%
%     
%     %flux versus u*
%     figure(1); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,Q_fQ_ust_bin_avg{j},Q_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     ylim([0 max(Q_all{i})]);
%     xlabel('u_{*} (m/s)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'Flux_Ust_fQ_',Sites{i},'.png'],'-dpng');
% 
%     %flux versus tau
%     figure(2); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tau_bins_mid,Q_fQ_tau_bin_avg{j},Q_fQ_tau_bin_std{j},Markers_fQ_bins{j});
%        %plot(tau_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     ylim([0 max(Q_all{i})]);
%     xlabel('\tau (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'Flux_Tau_fQ_',Sites{i},'.png'],'-dpng'); 
% 
    %flux versus tau_ex
    figure(3); clf; hold on;
    for j = ind_fQ_bins_full
       errorbar(tauex_bins_mid,Q_fQ_tauex_bin_avg{j},Q_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
       %plot(tauex_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
    ylim([0 max(Q_all{i})]);
    xlabel('\tau_{ex} (Pa)');
    ylabel('Q (g m^{-1} s^{-1})');
    title(Sites{i});
    h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'Flux_Tauex_fQ_',Sites{i},'.png'],'-dpng');
  
%     %eta versus u*
%     figure(4); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,eta_fQ_ust_bin_avg{j},eta_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},eta_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     xlabel('u_{*} (m/s)');
%     ylabel('\eta = (\int u^2-u^2_{th} dt) / (\int u^2 dt)');
%     ylim([0 max(eta_cal_all{i})]);
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'eta_Ust_fQ_',Sites{i},'.png'],'-dpng');
% 
    %eta versus tau_ex
    figure(5); clf; hold on;
    for j = ind_fQ_bins_full
       errorbar(tauex_bins_mid,eta_fQ_tauex_bin_avg{j},eta_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
       %plot(tauex_fQ_bin_values{i}{j},eta_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
%     [tauex_sort, ind_sort] = sort(tauex_all{i});
%     tauratio_sort = tauratio_all{i}(ind_sort);
%     plot(tauex_sort,tauratio_sort,'k');
    %plot tauex/tau versus tauex
    plot(tauex_bins_mid,tauratio_tauex_bin_avg{i},'k');
    xlabel('\tau_{ex} (Pa)');
    ylabel('\eta = (\int u^2-u^2_{th} dt) / (\int u^2 dt)');
    tauratio_plot = [0 0.75];
    ylim([0 max(eta_cal_all{i})]);
    title(Sites{i});
    legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
    legend_items{length(legend_items)+1} = '\tau_{ex}/\tau';
    h_legend = legend(legend_items,'Location','SouthEast');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'eta_Tauex_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %z0 versus u*
%     figure(6); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,z0_fQ_ust_bin_avg{j},z0_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},z0_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     set(gca,'yscale','log');
%     if i==1
%         ylim([1e-6 1e-2]);
%     elseif i==3
%        ylim([1e-7 1e-2]);
%     end
%     xlabel('u_{*} (m/s)');
%     ylabel('z_{0} (m)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'z0_Ust_fQ_',Sites{i},'.png'],'-dpng');
% 
%     %z0 versus tau_ex
%     figure(7); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tauex_bins_mid,z0_fQ_tauex_bin_avg{j},z0_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
%        %plot(tauex_fQ_bin_values{i}{j},z0_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     set(gca,'yscale','log');
%     if i==1
%         ylim([1e-6 1e-2]);
%     elseif i==3
%        ylim([1e-7 1e-2]);
%     end
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('z_{0} (m)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'z0_Tauex_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %alternate z0 versus u*
%     figure(8); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,z0alt_fQ_ust_bin_avg{j},z0alt_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},z0alt_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     set(gca,'yscale','log');
%     if i==1
%         ylim([1e-6 1e-2]);
%     elseif i==3
%        ylim([1e-5 1e-3]);
%     end
%     xlabel('u_{*} (m/s)');
%     ylabel('z_{0} (m)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'z0alt_Ust_fQ_',Sites{i},'.png'],'-dpng');
%     
    %alternate z0 versus tau_ex
    figure(9); clf; hold on;
    for j = 1:N_fQ_bins
       errorbar(tauex_bins_mid,z0alt_fQ_tauex_bin_avg{j},z0alt_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
       %plot(tauex_fQ_bin_values{i}{j},z0alt_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
    plot(tauex_bins_mid,1e-4+40*1e-4*tauex_bins_mid);
%    plot(tauex_bins_mid,1e-4+20*2e-4*tauex_bins_mid);
    set(gca,'yscale','log');
    if i==1
        ylim([1e-6 1e-2]);
    elseif i==3
       ylim([1e-5 1e-3]);
    end
    xlabel('\tau_{ex} (Pa)');
    ylabel('z_{0} (m)');
    title(Sites{i});
    h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'z0alt_Tauex_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %etaz0 versus u*
%     figure(10); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,etaz0_fQ_ust_bin_avg{j},etaz0_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},etaz0_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     xlabel('u_{*} (m/s)');
%     ylabel('\eta = 1-ln^2(z_d/z_s)/ln^2(z_d/z_0)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'etaz0_Ust_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %etaz0 versus tau_ex
%     figure(11); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tauex_bins_mid,etaz0_fQ_tauex_bin_avg{j},etaz0_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
%        %plot(tauex_fQ_bin_values{i}{j},etaz0_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('\eta = 1-ln^2(z_d/z_s)/ln^2(z_d/z_0)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'etaz0_tauex_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %flux versus eta*tau
%     figure(12); clf; hold on; %initialize plot
%     for j = ind_fQ_bins_full
%        errorbar(etatau_bins_mid,Q_fQ_etatau_bin_avg{j},Q_fQ_etatau_bin_std{j},Markers_fQ_bins{j});
%        %plot(eta_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     %xlim([0, min(etatau_bins_mid(etatau_bins_mid>max(eta_cal_all{i}.*tauRe_cal_all{i})))]);
%     ylim([0, max(Q_all{i})]);
%     xlabel('\eta\tau (Pa)');
%     ylabel('Q (g/m/s)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(gca,'FontSize',16);
%     title([Sites{i}]);
%     print([folder_Plots,'flux_etatau_',Sites{i},'.png'],'-dpng');
% 
%     %flux versus etaz0*tau
%     figure(13); clf; hold on; %initialize plot
%     for j = ind_fQ_bins_full
%        errorbar(etatau_bins_mid,Q_fQ_etaz0tau_bin_avg{j},Q_fQ_etaz0tau_bin_std{j},Markers_fQ_bins{j});
%        %plot(etaz0_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     xlabel('\eta_{z0}\tau (Pa)');
%     ylabel('Q (g/m/s)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(gca,'FontSize',16);
%     xlim([0, max(etatau_bins_mid(etatau_bins_mid<max(etaz0_all{i}.*tauRe_cal_all{i})))]);
%     title([Sites{i}]);
%     print([folder_Plots,'flux_etaz0tau_',Sites{i},'.png'],'-dpng');
    
    %ubar/ustd versus u*
    figure(14); clf; hold on; %initialize plot
    for j = ind_fQ_bins_full
       errorbar(ust_bins_mid,ubarustd_fQ_ust_bin_avg{j},ubarustd_fQ_ust_bin_std{j},Markers_fQ_bins{j});
       %plot(ust_fQ_bin_values{i}{j},ubarustd_fQ_ust_bin_values{i}{j},Markers_fQ_bins{j});
    end
    xlabel('u_{*} (m/s)');
    ylabel('\mu_{u}/\sigma_{u}');
    title(Sites{i});
    h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'UbarUstd_Ust_fQ_',Sites{i},'.png'],'-dpng');
end

%flux frequency versus eta
figure(15); clf; hold on; %initialize plot
for i=1:N_Sites
    errorbar(eta_bins_mid,fQ_eta_bin_avg{i},fQ_eta_bin_std{i},Markers{i});
end
xlabel('\eta');
ylabel('f_{Q}');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'fQ_eta.png'],'-dpng');

%ustth versus fQ
figure(16); clf; hold on; %initialize plot
for i=1:N_Sites
    errorbar(fQ_bins_alt_mid,ustth_fQ_bin_avg{i},ustth_fQ_bin_std{i},Markers{i});
end
xlabel('f_{Q}');
ylabel('u_{*,th,TFEM}');
legend(Sites,'Location','SouthEast');
set(gca,'FontSize',16);
print([folder_Plots,'ustth_fQ.png'],'-dpng');

%tauth versus fQ
figure(17); clf; hold on; %initialize plot
%for i=1:N_Sites
    %errorbar(fQ_bins_alt_mid,tauth_fQ_bin_avg{i},tauth_fQ_bin_std{i},Markers{i});
%end
errorbar(fQ_bins_alt_mid,tauth_fQ_bin_avg{3},tauth_fQ_bin_std{3},'ob');
[a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_bins_alt_mid', tauth_fQ_bin_avg{3}, fQ_bins_alt_max'-fQ_bins_alt_min', tauth_fQ_bin_SE{3});
tauft = a;
tauit = a+b;
%plot([0 1],[0.14, 0.1]);
plot(fQ_bins_alt_mid,tauth_fit,'k','LineWidth',2);
plot(fQ_bins_alt_mid,tauth_fit+sigma_tauth_fit,'k--',fQ_bins_alt_mid,tauth_fit-sigma_tauth_fit,'k--');
xlabel('f_{Q}');
ylabel('\tau_{th,TFEM}');
title('Oceano');
%legend(Sites,'Location','SouthEast');
set(gca,'FontSize',16);
print([folder_Plots,'tauth_fQ.png'],'-dpng');