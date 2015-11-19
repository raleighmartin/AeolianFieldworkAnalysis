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
rho_a = 1.23; %air density, kg/m^3
tau_it = [0.13 0.13 0.087]; %impact threshold stress (Pa)
tau_ft = [0.18 0.18 0.13]; %fluid threshold stress (Pa)
ust_it = sqrt(tau_it./rho_a); %impact threshold shear velocity (m/s)
ust_ft = sqrt(tau_ft./rho_a); %fluid threshold shear velocity (m/s)
z0f = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m) at threshold
% tau_it = [0.1 0.1 0.1]; %impact threshold stress (Pa)
% tau_ft = [0.14 0.14 0.14]; %fluid threshold stress (Pa)
% ust_it = [0.35, 0.28, 0.28]; %assumed u* threshold (m/s) for each site
% tau_it = rho*ust_it.^2; %assumed tau threshold (Pa) for each site
% z0f = [9e-5, 1e-4, 2e-4]; %aerodynamic roughness length (m) at threshold
% z0f = [1e-5, 1e-5, 1e-5]; %aerodynamic roughness length (m) at threshold
% z0f = [3.7e-5, 3.7e-5, 8.5e-5]; %aerodynamic roughness length (m) at threshold

%% CALCULATIONS

%% COMPUTE U*TH_TFEM
ustth_TFEM_all = cell(N_Sites,1);
for i=1:N_Sites
    ustth_TFEM_all{i} = sqrt(tauth_TFEM_all{i}/rho_a);
end

%% COMPUTE USTEX, TAUEX, TAURATIO BASED ON THRESHOLD BY SITE
ustex_all = cell(N_Sites,1);
tauex_all = cell(N_Sites,1);
tauratio_all = cell(N_Sites,1);
for i=1:N_Sites
    ustex_all{i} = ustRe_all{i}-ust_it(i);
    tauex_all{i} = tauRe_all{i}-tau_it(i);
    tauratio_all{i} = tauex_all{i}./tauRe_all{i};
end

% %% COMPUTE USTEX AND TAUEX BASED ON TFEM THRESHOLD
% ustex_all = cell(N_Sites,1);
% tauex_all = cell(N_Sites,1);
% tauratio_all = cell(N_Sites,1);
% for i=1:N_Sites
%     ustex_all{i} = ustRe_all{i}-ustth_TFEM_all{i};
%     tauex_all{i} = tauRe_all{i}-tauth_TFEM_all{i};
%     tauratio_all{i} = tauex_all{i}./tauRe_all{i};
% end
%
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
%         ustex_all{i} = ustRe_all{i}-0.3474;
%         tauex_all{i} = tauRe_all{i}-0.14841;
%         tauratio_all{i} = tauex_all{i}/0.14841;
%     %apply empirical threshold - Rancho Guadalupe - all dates
%     elseif i == 2;
%         ustex_all{i} = ustRe_all{i}-0.2839;
%         tauex_all{i} = tauRe_all{i}-0.099111;
%         tauratio_all{i} = tauex_all{i}/0.099111;
%     elseif i == 3;
%         %apply empirical threshold - Oceano - May 15-19
%         date_ind = intersect(find(date_all{i}>=datetime(2015,5,15)),find(date_all{i}<=datetime(2015,5,19)));
%         ustex_all{i}(date_ind) = ustRe_all{i}(date_ind)-0.2794;
%         tauex_all{i}(date_ind) = tauRe_all{i}(date_ind)-0.096036;
%         tauratio_all{i}(date_ind) = tauex_all{i}(date_ind)/0.096036;
%         %apply empirical threshold - Oceano - May 23-31
%         date_ind = intersect(find(date_all{i}>=datetime(2015,5,23)),find(date_all{i}<=datetime(2015,5,31)));
%         ustex_all{i}(date_ind) = ustRe_all{i}(date_ind)-0.2903;
%         tauex_all{i}(date_ind) = tauRe_all{i}(date_ind)-0.10364;
%         tauratio_all{i}(date_ind) = tauex_all{i}(date_ind)/0.10364;
%         %apply empirical threshold - Oceano - June 1-4
%         date_ind = intersect(find(date_all{i}>=datetime(2015,6,1)),find(date_all{i}<=datetime(2015,6,4)));
%         ustex_all{i}(date_ind) = ustRe_all{i}(date_ind)-0.2792;
%         tauex_all{i}(date_ind) = tauRe_all{i}(date_ind)-0.095865;
%         tauratio_all{i}(date_ind) = tauex_all{i}(date_ind)/0.095865;
%     end
% end

%% CREATE BINS
%set frequency bins
fQ_bins_min = [0 0.3 0.6 0.9];
fQ_bins_max = [0.3 0.6 0.9 1];
%fQ_bins_min = [0 0.3 0.6 0.85 1];
%fQ_bins_max = [0.3 0.6 0.85 1 Inf];
fQ_bins_mid = mean([fQ_bins_min; fQ_bins_max]);
fQ_bins_legend_all = {'f_{Q} = 0-0.3','f_{Q} = 0.3-0.6', 'f_{Q} = 0.6-0.9', 'f_{Q} = 0.9-1'};
%fQ_bins_legend_all = {'f_{Q} = 0-0.3','f_{Q} = 0.3-0.6', 'f_{Q} = 0.6-0.85', 'f_{Q} = 0.85<1', 'f_{Q}=1'};
N_fQ_bins = length(fQ_bins_mid);
Markers_fQ_bins = {'rx','bv','g^','ko'};
%Markers_fQ_bins = {'rx','bv','gd','m^','ko'};

%create u* bins
ust_bins_min = 0.025:.025:0.55;
ust_bins_max = 0.05:.025:0.575;
ust_bins_mid = mean([ust_bins_min;ust_bins_max]);
N_ust_bins = length(ust_bins_mid);

%tau_bins
%tau_bins_min = 0:.025:0.375; %larger bins
%tau_bins_max = 0.025:.025:0.4; %larger bins
tau_bins_min = 0:0.02:0.38; %medium bins
tau_bins_max = 0.02:0.02:0.4; %medium bins
%tau_bins_min = 0:.0125:0.3875; %smaller bins
%tau_bins_max = 0.0125:.0125:0.4; %smaller bins
tau_bins_mid = mean([tau_bins_min;tau_bins_max]);
N_tau_bins = length(tau_bins_mid);

%tauex_bins
tauex_bins_min = -0.15:0.025:0.275; %larger bins
tauex_bins_max = -0.125:0.025:0.3; %larger bins
%tauex_bins_min = -0.1:0.0125:0.2375; %smaller bins
%tauex_bins_max = -0.0875:0.0125:0.25; %smaller bins
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
zs_fQ_bin_values = cell(N_Sites,1); %zs into frequency bins
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
    zs_fQ_bin_values{i} = cell(N_fQ_bins,1);
    ubarustd_fQ_bin_values{i} = cell(N_fQ_bins,1);
    
    %go through frequency bins, get values
    for j=1:N_fQ_bins
        bin_ind = find(fQ1_all{i}>=fQ_bins_min(j)&fQ1_all{i}<=fQ_bins_max(j)); %use 1 second frequencies for binning, soft upper limit
        %bin_ind = find(fQ1_all{i}>=fQ_bins_min(j)&fQ1_all{i}<fQ_bins_max(j)); %use 1 second frequencies for binning, hard upper limit
        %bin_ind = find(fQ_all{i}(:,1)>=fQ_bins_min(j)&fQ_all{i}(:,1)<=fQ_bins_max(j)); %use specific time interval for binning
        bin_ind = setdiff(bin_ind,ind_outlier{i}); %remove outliers specified above
        fQ_bins_legend{i}{j} = strcat(fQ_bins_legend_all{j},'; N = ',int2str(length(bin_ind)));
        if ~isempty(bin_ind)
            ust_fQ_bin_values{i}{j} = ustRe_all{i}(bin_ind);
            tau_fQ_bin_values{i}{j} = tauRe_all{i}(bin_ind);
            tauex_fQ_bin_values{i}{j} = tauex_all{i}(bin_ind);
            Q_fQ_bin_values{i}{j} = Q_all{i}(bin_ind);
            eta_fQ_bin_values{i}{j} = eta_fQ_constthr_all{i}(bin_ind); %use eta based on frequency below constant threshold method
            %eta_fQ_bin_values{i}{j} = eta_zs_constthr_all{i}(bin_ind); %use eta based on constant threshold
            %eta_fQ_bin_values{i}{j} = eta_zs_TFEMthr_all{i}(bin_ind); %use eta based on TFEM threshold 
            zs_fQ_bin_values{i}{j} = zs_all{i}(bin_ind);
            ubarustd_fQ_bin_values{i}{j} = u_bar_all{i}(bin_ind)./u_std_all{i}(bin_ind);
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
        bin_ind = find(eta_fQ_constthr_all{i}>=eta_bins_min(j)&eta_fQ_constthr_all{i}<=eta_bins_max(j));
        %bin_ind = find(eta_zs_constthr_all{i}>=eta_bins_min(j)&eta_zs_constthr_all{i}<=eta_bins_max(j));
        %bin_ind = find(eta_zs_TFEMthr_all{i}>=eta_bins_min(j)&eta_zs_TFEMthr_all{i}<=eta_bins_max(j));
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


%% FOR EACH SITE, PERFORM BINNING BY ALTERNATE (SMALLER) fQ BINS (more bins, for direct comparisons to frequency)
fQ_altbins_min = 0:0.05:0.95;
fQ_altbins_max = 0.05:0.05:1;
fQ_altbins_mid = mean([fQ_altbins_min; fQ_altbins_max]);
N_fQ_altbins = length(fQ_altbins_mid);

%separate u*th based on TFEM into frequency bins
ustth_fQ_altbin_values = cell(N_Sites,1);
ustth_fQ_altbin_avg = cell(N_Sites,1);
ustth_fQ_altbin_std = cell(N_Sites,1);
ustth_fQ_altbin_SE = cell(N_Sites,1);
tauth_fQ_altbin_values = cell(N_Sites,1);
tauth_fQ_altbin_avg = cell(N_Sites,1);
tauth_fQ_altbin_std = cell(N_Sites,1);
tauth_fQ_altbin_SE = cell(N_Sites,1);

f_uaboveft_values = cell(N_Sites,1);
f_uaboveft_avg = cell(N_Sites,1);
f_uaboveft_std = cell(N_Sites,1);
f_uaboveft_SE = cell(N_Sites,1);

f_ubelowit_values = cell(N_Sites,1);
f_ubelowit_avg = cell(N_Sites,1);
f_ubelowit_std = cell(N_Sites,1);
f_ubelowit_SE = cell(N_Sites,1);

f_ubetween_fromabove_values = cell(N_Sites,1);
f_ubetween_fromabove_avg = cell(N_Sites,1);
f_ubetween_fromabove_std = cell(N_Sites,1);
f_ubetween_fromabove_SE = cell(N_Sites,1);

f_ubetween_frombelow_values = cell(N_Sites,1);
f_ubetween_frombelow_avg = cell(N_Sites,1);
f_ubetween_frombelow_std = cell(N_Sites,1);
f_ubetween_frombelow_SE = cell(N_Sites,1);

for i = 1:N_Sites
    ustth_fQ_altbin_values{i} = cell(N_fQ_altbins,1);
    ustth_fQ_altbin_avg{i} = zeros(N_fQ_altbins,1)*NaN;
    ustth_fQ_altbin_std{i} = zeros(N_fQ_altbins,1)*NaN;
    ustth_fQ_altbin_SE{i} = zeros(N_fQ_altbins,1)*NaN;
    tauth_fQ_altbin_values{i} = cell(N_fQ_altbins,1);
    tauth_fQ_altbin_avg{i} = zeros(N_fQ_altbins,1)*NaN;
    tauth_fQ_altbin_std{i} = zeros(N_fQ_altbins,1)*NaN;
    tauth_fQ_altbin_SE{i} = zeros(N_fQ_altbins,1)*NaN;
    
    f_uaboveft_values{i} = cell(N_fQ_altbins,1);
    f_uaboveft_avg{i} = zeros(N_fQ_altbins,1)*NaN;
    f_uaboveft_std{i} = zeros(N_fQ_altbins,1)*NaN;
    f_uaboveft_SE{i} = zeros(N_fQ_altbins,1)*NaN;

    f_ubelowit_values{i} = cell(N_fQ_altbins,1);
    f_ubelowit_avg{i} = zeros(N_fQ_altbins,1)*NaN;
    f_ubelowit_std{i} = zeros(N_fQ_altbins,1)*NaN;
    f_ubelowit_SE{i} = zeros(N_fQ_altbins,1)*NaN;

    f_ubetween_fromabove_values{i} = cell(N_fQ_altbins,1);
    f_ubetween_fromabove_avg{i} = zeros(N_fQ_altbins,1)*NaN;
    f_ubetween_fromabove_std{i} = zeros(N_fQ_altbins,1)*NaN;
    f_ubetween_fromabove_SE{i} = zeros(N_fQ_altbins,1)*NaN;

    f_ubetween_frombelow_values{i} = cell(N_fQ_altbins,1);
    f_ubetween_frombelow_avg{i} = zeros(N_fQ_altbins,1)*NaN;
    f_ubetween_frombelow_std{i} = zeros(N_fQ_altbins,1)*NaN;
    f_ubetween_frombelow_SE{i} = zeros(N_fQ_altbins,1)*NaN;
    
    %get u*th's / tauth's for fQ alt bins
    for j = 1:N_fQ_altbins
        bin_ind = find(fQ1_all{i}>=fQ_altbins_min(j)&fQ1_all{i}<=fQ_altbins_max(j)); %use 1 second frequencies for binning
        ustth_fQ_altbin_values{i}{j} = ustth_TFEM_all{i}(bin_ind);
        ustth_fQ_altbin_values{i}{j} = ustth_fQ_altbin_values{i}{j}(~isnan(ustth_fQ_altbin_values{i}{j}));
        ustth_fQ_altbin_avg{i}(j) = mean(ustth_fQ_altbin_values{i}{j});
        ustth_fQ_altbin_std{i}(j) = std(ustth_fQ_altbin_values{i}{j});
        ustth_fQ_altbin_SE{i}(j) = std(ustth_fQ_altbin_values{i}{j})/sqrt(length(~isnan(ustth_fQ_altbin_values{i}{j})));
        
        tauth_fQ_altbin_values{i}{j} = tauth_TFEM_all{i}(bin_ind);
        tauth_fQ_altbin_values{i}{j} = tauth_fQ_altbin_values{i}{j}(~isnan(tauth_fQ_altbin_values{i}{j}));
        tauth_fQ_altbin_avg{i}(j) = mean(tauth_fQ_altbin_values{i}{j});
        tauth_fQ_altbin_std{i}(j) = std(tauth_fQ_altbin_values{i}{j});
        tauth_fQ_altbin_SE{i}(j) = std(tauth_fQ_altbin_values{i}{j})/sqrt(length(~isnan(tauth_fQ_altbin_values{i}{j})));
        
        f_uaboveft_values{i}{j} = f_uaboveft_all{i}(bin_ind);
        f_uaboveft_values{i}{j} = f_uaboveft_values{i}{j}(~isnan(f_uaboveft_values{i}{j}));
        f_uaboveft_avg{i}(j) = mean(f_uaboveft_values{i}{j});
        f_uaboveft_std{i}(j) = std(f_uaboveft_values{i}{j});
        f_uaboveft_SE{i}(j) = std(f_uaboveft_values{i}{j})/sqrt(length(~isnan(f_uaboveft_values{i}{j})));
                
        f_ubelowit_values{i}{j} = f_ubelowit_all{i}(bin_ind);
        f_ubelowit_values{i}{j} = f_ubelowit_values{i}{j}(~isnan(f_ubelowit_values{i}{j}));
        f_ubelowit_avg{i}(j) = mean(f_ubelowit_values{i}{j});
        f_ubelowit_std{i}(j) = std(f_ubelowit_values{i}{j});
        f_ubelowit_SE{i}(j) = std(f_ubelowit_values{i}{j})/sqrt(length(~isnan(f_ubelowit_values{i}{j})));
        
        f_ubetween_fromabove_values{i}{j} = f_ubetween_fromabove_all{i}(bin_ind);
        f_ubetween_fromabove_values{i}{j} = f_ubetween_fromabove_values{i}{j}(~isnan(f_ubetween_fromabove_values{i}{j}));
        f_ubetween_fromabove_avg{i}(j) = mean(f_ubetween_fromabove_values{i}{j});
        f_ubetween_fromabove_std{i}(j) = std(f_ubetween_fromabove_values{i}{j});
        f_ubetween_fromabove_SE{i}(j) = std(f_ubetween_fromabove_values{i}{j})/sqrt(length(~isnan(f_ubetween_fromabove_values{i}{j})));
        
        f_ubetween_frombelow_values{i}{j} = f_ubetween_frombelow_all{i}(bin_ind);
        f_ubetween_frombelow_values{i}{j} = f_ubetween_frombelow_values{i}{j}(~isnan(f_ubetween_frombelow_values{i}{j}));
        f_ubetween_frombelow_avg{i}(j) = mean(f_ubetween_frombelow_values{i}{j});
        f_ubetween_frombelow_std{i}(j) = std(f_ubetween_frombelow_values{i}{j});
        f_ubetween_frombelow_SE{i}(j) = std(f_ubetween_frombelow_values{i}{j})/sqrt(length(~isnan(f_ubetween_frombelow_values{i}{j})));
    end 
end

%% PERFORM SUB-BINNING WITHIN FREQUENCY BINS, CREATE PLOTS FOR EACH SITE
for i = 3;
%for i = 1:N_Sites
    
    %get indices of fQ bins that are not empty
    ind_fQ_bins_full = find(1-cellfun(@isempty,Q_fQ_bin_values{i}))';
    
    %fluxes and other values into u* bins within fQ bins
    Q_fQ_ust_bin_avg = cell(N_fQ_bins,1);
    Q_fQ_ust_bin_std = cell(N_fQ_bins,1);
    Q_fQ_ust_bin_SE = cell(N_fQ_bins,1);
    eta_fQ_ust_bin_avg = cell(N_fQ_bins,1);
    eta_fQ_ust_bin_std = cell(N_fQ_bins,1);
    eta_fQ_ust_bin_SE = cell(N_fQ_bins,1);
    zs_fQ_ust_bin_avg = cell(N_fQ_bins,1);
    zs_fQ_ust_bin_std = cell(N_fQ_bins,1);
    zs_fQ_ust_bin_SE = cell(N_fQ_bins,1);
    ubarustd_fQ_ust_bin_avg = cell(N_fQ_bins,1);
    ubarustd_fQ_ust_bin_std = cell(N_fQ_bins,1);
    ubarustd_fQ_ust_bin_SE = cell(N_fQ_bins,1);
    
    %fluxes into tau bins within fQ bins
    Q_fQ_tau_bin_avg = cell(N_fQ_bins,1);
    Q_fQ_tau_bin_std = cell(N_fQ_bins,1);
    Q_fQ_tau_bin_SE = cell(N_fQ_bins,1);
    
    %fluxes and other values into tauex bins within fQ bins
    Q_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
    Q_fQ_tauex_bin_std = cell(N_fQ_bins,1);
    Q_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
    eta_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
    eta_fQ_tauex_bin_std = cell(N_fQ_bins,1);
    eta_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
    zs_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
    zs_fQ_tauex_bin_std = cell(N_fQ_bins,1);
    zs_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
    
    %fluxes into eta*tau bins within fQ bins
    Q_fQ_etatau_bin_avg = cell(N_fQ_bins,1);
    Q_fQ_etatau_bin_std = cell(N_fQ_bins,1);
    Q_fQ_etatau_bin_SE = cell(N_fQ_bins,1);
    
    %go through each fQ bin
    for j = 1:N_fQ_bins
        
        %initialize u* bins for each fQ bin
        Q_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        Q_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        Q_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        eta_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        eta_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        eta_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
        zs_fQ_ust_bin_avg{j} = zeros(N_ust_bins,1)*NaN;
        zs_fQ_ust_bin_std{j} = zeros(N_ust_bins,1)*NaN;
        zs_fQ_ust_bin_SE{j} = zeros(N_ust_bins,1)*NaN;
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
        zs_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
        zs_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
        zs_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
                
        %initialize etatau bins for each fQ bin
        Q_fQ_etatau_bin_avg{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etatau_bin_std{j} = zeros(N_etatau_bins,1)*NaN;
        Q_fQ_etatau_bin_SE{j} = zeros(N_etatau_bins,1)*NaN;
              
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
                
                zs_fQ_ust_bin_values = zs_fQ_bin_values{i}{j}(bin_ind);
                zs_fQ_ust_bin_values = zs_fQ_ust_bin_values(~isnan(zs_fQ_ust_bin_values));
                zs_fQ_ust_bin_avg{j}(k) = exp(mean(log(zs_fQ_ust_bin_values))); %get mean in log space
                zs_fQ_ust_bin_SE{j}(k) = std(zs_fQ_ust_bin_values)/sqrt(length(~isnan(zs_fQ_ust_bin_values)));
                zs_fQ_ust_bin_std{j}(k) = std(zs_fQ_ust_bin_values);
                                
                ubarustd_fQ_ust_bin_values = ubarustd_fQ_bin_values{i}{j}(bin_ind);
                ubarustd_fQ_ust_bin_values = ubarustd_fQ_ust_bin_values(~isnan(ubarustd_fQ_ust_bin_values));
                ubarustd_fQ_ust_bin_avg{j}(k) = mean(ubarustd_fQ_ust_bin_values);
                ubarustd_fQ_ust_bin_SE{j}(k) = std(ubarustd_fQ_ust_bin_values)/sqrt(length(~isnan(ubarustd_fQ_ust_bin_values)));
                ubarustd_fQ_ust_bin_std{j}(k) = std(ubarustd_fQ_ust_bin_values);
            end
        end
        
        %go through each tau bin
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
                
                zs_fQ_tauex_bin_values = zs_fQ_bin_values{i}{j}(bin_ind);
                eta_fQ_tauex_bin_values = eta_fQ_tauex_bin_values(~isnan(eta_fQ_tauex_bin_values));
                zs_fQ_tauex_bin_avg{j}(k) = exp(mean(log(zs_fQ_tauex_bin_values))); %get mean in log space
                zs_fQ_tauex_bin_SE{j}(k) = std(zs_fQ_tauex_bin_values)/sqrt(length(~isnan(zs_fQ_tauex_bin_values)));
                zs_fQ_tauex_bin_std{j}(k) = std(zs_fQ_tauex_bin_values);
                
            end
        end
        
        %go through each etatau bin
        for k=1:N_etatau_bins
            bin_ind = find(eta_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}>=etatau_bins_min(k)&...
                eta_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}<=etatau_bins_max(k));
            if ~isempty(bin_ind)
                Q_fQ_etatau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
                Q_fQ_etatau_bin_values = Q_fQ_etatau_bin_values(~isnan(Q_fQ_etatau_bin_values));
                Q_fQ_etatau_bin_avg{j}(k) = mean(Q_fQ_etatau_bin_values);
                Q_fQ_etatau_bin_SE{j}(k) = std(Q_fQ_etatau_bin_values)/sqrt(length(~isnan(Q_fQ_etatau_bin_values)));
                Q_fQ_etatau_bin_std{j}(k) = std(Q_fQ_etatau_bin_values);
            end
        end    
    end

    %% PLOTS %%
        
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
    %flux versus tau
    figure(2); clf; hold on;
    for j = ind_fQ_bins_full
       errorbar(tau_bins_mid,Q_fQ_tau_bin_avg{j},Q_fQ_tau_bin_SE{j},Markers_fQ_bins{j});
       %errorbar(tau_bins_mid,Q_fQ_tau_bin_avg{j},Q_fQ_tau_bin_std{j},Markers_fQ_bins{j});
       %plot(tau_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
    
    %fit to continuous flux points
    ind_good = find(~isnan(Q_fQ_tau_bin_avg{end}));
    tau_fit = tau_bins_mid(ind_good)';
    sigma_tau_fit = ones(size(ind_good))*mode(tau_bins_max-tau_bins_min);
    Q_fit = Q_fQ_tau_bin_avg{end}(ind_good);
    sigma_Q_fit = Q_fQ_tau_bin_SE{end}(ind_good);
    [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(tau_fit, Q_fit, sigma_tau_fit, sigma_Q_fit);
    CQ_pred = b;
    sigma_CQ_pred = sigma_b;
    tau_it_pred = -a/CQ_pred;
    sigma_tau_it_pred = sigma_a/CQ_pred;
    tau_conf_min = [tau_it_pred+sigma_tau_it_pred; tau_fit];
    tau_conf_max = [tau_it_pred-sigma_tau_it_pred; tau_fit];
    Q_conf_min = [0; Q_pred-sigma_Q_pred];
    Q_conf_max = [0; Q_pred+sigma_Q_pred];
    plot([tau_it_pred max(tau_fit)],[0 CQ_pred*(max(tau_fit)-tau_it_pred)]);
    plot(tau_conf_min,Q_conf_min,'b--',tau_conf_max,Q_conf_max,'b--');
    
    ylim([0 max(Q_all{i})]);
    xlabel('\tau (Pa)');
    ylabel('Q (g m^{-1} s^{-1})');
    title(Sites{i});
    legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
    legend_items{length(legend_items)+1} = 'Q(\tau-\tau_{it}) fit';
    h_legend = legend(legend_items,'Location','NorthWest');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'Flux_Tau_fQ_',Sites{i},'.png'],'-dpng'); 

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
%     ylim([0 max(eta_zs_all{i})]);
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'eta_Ust_fQ_',Sites{i},'.png'],'-dpng');

    %eta versus tau_ex
    figure(5); clf; hold on;
    for j = ind_fQ_bins_full
       errorbar(tauex_bins_mid,eta_fQ_tauex_bin_avg{j},eta_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
       %plot(tauex_fQ_bin_values{i}{j},eta_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
    plot(tauex_bins_mid,tauratio_tauex_bin_avg{i},'k'); %plot tauex/tau versus tauex
    xlabel('\tau_{ex} (Pa)');
    ylabel('\eta = (\int u^2-u^2_{th} dt) / (\int u^2 dt)');
    tauratio_plot = [0 0.75];
    ylim([0 max(eta_fQ_constthr_all{i})]);
    %ylim([0 max(eta_zs_constthr_all{i})]);
    %ylim([0 max(eta_zs_TFEMthr_all{i})]);
    title(Sites{i});
    legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
    legend_items{length(legend_items)+1} = '\tau_{ex}/\tau';
    h_legend = legend(legend_items,'Location','SouthEast');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'eta_Tauex_fQ_',Sites{i},'.png'],'-dpng');
    
%     %zs versus u*
%     figure(6); clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,zs_fQ_ust_bin_avg{j},zs_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},zs_fQ_bin_values{i}{j},Markers_fQ_bins{j});
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
%     print([folder_Plots,'zs_Ust_fQ_',Sites{i},'.png'],'-dpng');

    %zs versus tau_ex
    figure(7); clf; hold on;
    for j = ind_fQ_bins_full
       errorbar(tauex_bins_mid,zs_fQ_tauex_bin_avg{j},zs_fQ_tauex_bin_std{j},Markers_fQ_bins{j});
       %plot(tauex_fQ_bin_values{i}{j},zs_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
    set(gca,'yscale','log');
    if i==1
        ylim([1e-6 1e-2]);
    elseif i==3
       ylim([1e-7 1e-2]);
    end
    xlabel('\tau_{ex} (Pa)');
    ylabel('z_{0} (m)');
    title(Sites{i});
    h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
    set(gca,'FontSize',16);
    set(h_legend,'FontSize',16);
    print([folder_Plots,'zs_Tauex_fQ_',Sites{i},'.png'],'-dpng');
    
    %flux versus eta*tau
    figure(8); clf; hold on; %initialize plot
    for j = ind_fQ_bins_full
       errorbar(etatau_bins_mid,Q_fQ_etatau_bin_avg{j},Q_fQ_etatau_bin_std{j},Markers_fQ_bins{j});
       %plot(eta_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j},Q_fQ_bin_values{i}{j},Markers_fQ_bins{j});
    end
    %xlim([0, min(etatau_bins_mid(etatau_bins_mid>max(eta_zs_all{i}.*tauRe_all{i})))]);
    ylim([0, max(Q_all{i})]);
    xlabel('\eta\tau (Pa)');
    ylabel('Q (g/m/s)');
    title(Sites{i});
    h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
    set(gca,'FontSize',16);
    set(gca,'FontSize',16);
    title([Sites{i}]);
    print([folder_Plots,'flux_etatau_',Sites{i},'.png'],'-dpng');

%     %ubar/ustd versus u*
%     figure(9); clf; hold on; %initialize plot
%     for j = ind_fQ_bins_full
%        errorbar(ust_bins_mid,ubarustd_fQ_ust_bin_avg{j},ubarustd_fQ_ust_bin_std{j},Markers_fQ_bins{j});
%        %plot(ust_fQ_bin_values{i}{j},ubarustd_fQ_ust_bin_values{i}{j},Markers_fQ_bins{j});
%     end
%     xlabel('u_{*} (m/s)');
%     ylabel('\mu_{u}/\sigma_{u}');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'UbarUstd_Ust_fQ_',Sites{i},'.png'],'-dpng');
    
    %tauth versus fQ
    figure(10); clf; hold on; %initialize plot
    errorbar(fQ_altbins_mid,tauth_fQ_altbin_avg{i},tauth_fQ_altbin_std{i},'ob');
    %perform fit for Oceano only
    if i==3
        ind_fit = intersect(intersect(find(fQ_altbins_mid>0.1),find(fQ_altbins_mid<0.85)),find(~isnan(tauth_fQ_altbin_avg{i})));
        [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_altbins_mid(ind_fit)', tauth_fQ_altbin_avg{i}(ind_fit), fQ_altbins_max(ind_fit)'-fQ_altbins_min(ind_fit)', tauth_fQ_altbin_SE{i}(ind_fit));
        tauft = a
        tauit = a+b
        plot([0 1],[tauft tauit],'k','LineWidth',2);
        plot(fQ_altbins_mid(ind_fit),tauth_fit+sigma_tauth_fit,'k--',fQ_altbins_mid(ind_fit),tauth_fit-sigma_tauth_fit,'k--');
    end
    xlabel('f_{Q}');
    ylabel('\tau_{th,TFEM}');
    title(Sites{i});
    set(gca,'FontSize',16);
    print([folder_Plots,'tauth_fQ_',Sites{i},'.png'],'-dpng');
    
    %frequency of u ranges versus fQ
    figure(11); clf; hold on; %initialize plot
    xlabel('f_{Q}');
    ylabel('f_{u range}');
    errorbar(fQ_altbins_mid,f_uaboveft_avg{i},f_uaboveft_std{i},'og');
    errorbar(fQ_altbins_mid,f_ubelowit_avg{i},f_ubelowit_std{i},'xr');
    errorbar(fQ_altbins_mid,f_ubetween_frombelow_avg{i},f_ubetween_frombelow_std{i},'^b');
    errorbar(fQ_altbins_mid,f_ubetween_fromabove_avg{i},f_ubetween_fromabove_std{i},'vb');
    legend('above ft','below it','from below','from above');
    title(Sites{i});
    set(gca,'FontSize',16);
    print([folder_Plots,'tauth_f_urange_',Sites{i},'.png'],'-dpng');
end

% %flux frequency versus eta
% figure(12); clf; hold on; %initialize plot
% for i=1:N_Sites
%     errorbar(eta_bins_mid,fQ_eta_bin_avg{i},fQ_eta_bin_std{i},Markers{i});
% end
% xlabel('\eta');
% ylabel('f_{Q}');
% legend(Sites,'Location','NorthWest');
% set(gca,'FontSize',16);
% print([folder_Plots,'fQ_eta.png'],'-dpng');
% 
% %ustth versus fQ
% figure(13); clf; hold on; %initialize plot
% for i=1:N_Sites
%     errorbar(fQ_bins_alt_mid,ustth_fQ_bin_avg{i},ustth_fQ_bin_std{i},Markers{i});
% end
% xlabel('f_{Q}');
% ylabel('u_{*,th,TFEM}');
% legend(Sites,'Location','SouthEast');
% set(gca,'FontSize',16);
% print([folder_Plots,'ustth_fQ.png'],'-dpng');