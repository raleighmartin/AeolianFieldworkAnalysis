%% ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;
close all;

%parameters
rho_a = 1.23; %air density (kg/m^3)
g = 9.8; %gravity (m/s^2)

%suppress figures
set(0,'DefaultFigureVisible', 'off');

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_AnalysisData = '../AnalysisData/'; %folder for storing data output
folder_Plots = '../PlotOutput/FluxFrequency/'; %folder for plots

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows_all'));

%load external data
load(strcat(folder_AnalysisData,'GreeleyNamikasData')); %literature data
Marker_Greeley96 = 'md';
Marker_Namikas03 = 'k^';
Marker_Farrell12 = 'c+';

%decide which etas to use
eta_active_all = eta_active_fQ_1s_TFEMthr_all; %eta values to use when specifying "active"
eta_inactive_all = eta_inactive_fQ_1s_TFEMthr_all; %eta values to use when specifying "inactive"
eta_all = eta_active_all; %use "active" eta values when "active" vs. "inactive" is not specified

%get info about number of sites
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};
LineColors = {'b','r','g'};

%% PARAMETERS
%tauit = [0.1206 0.0905 0.0853]; %impact threshold stress (Pa) - TFEM
tauit = [0.1484 0.1078 0.0919]; %impact threshold stress (Pa) - fit
tauft = [0.1787 0.1631 0.1299]; %fluid threshold stress (Pa)
d50_site = [0.57 0.53 0.40]; %mean grain diameter for site (mm)
z0 = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m) at threshold
theta_limit = 20; %highest angle (deg) for Oceano observations
zL_limit = -0.2; %limit for stability parameter (z/L) for all observations

%% REMOVE / FIX BAD POINTS
variable_list = who('variables','*all');
N_variables = length(variable_list);

% remove points with angle too large (Oceano only)
ind_theta = find(theta_all{3}<=theta_limit);
for j = 1:N_variables
    eval([variable_list{j},'{3}=',variable_list{j},'{3}(ind_theta);']);
end

% remove points with too much instability
for i = 1:N_Sites
    ind_zL = find(zL_all{i}>=zL_limit);
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_zL);']);
    end
end

% remove points with no flux measurement
for i = 1:N_Sites
    ind_Q_notempty = find(cellfun(@length,q_all{i})>0);
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_Q_notempty);']);
    end
end
    
% redo profiles for NaN flux points - using only qz's above 0
for i = 1:N_Sites
    ind_Q_NaN = find(isnan(Q_all{i})|isnan(sigma_Q_all{i}));
    for j = 1:length(ind_Q_NaN);
        ind_fit = find(q_all{i}{ind_Q_NaN(j)}>0);
        if length(ind_fit)>=3
            q_profile = q_all{i}{ind_Q_NaN(j)}(ind_fit);
            zW_profile = zW_all{i}{ind_Q_NaN(j)}(ind_fit);
            sigma_q_profile = sigma_q_all{i}{ind_Q_NaN(j)}(ind_fit); %use later after re-running CreateStressFluxWindows
            sigma_zW_profile = sigma_zW_all{i}{ind_Q_NaN(j)}(ind_fit); %use later after re-running CreateStressFluxWindows
            [q0,zq,sigma_q0,sigma_zq] = qz_profilefit(q_profile,zW_profile,sigma_q_profile,sigma_zW_profile);
            Q = q0*zq; %get total flux [g/m/s]
            sigma_Q = sqrt((sigma_q0*zq)^2+(sigma_zq*q0)^2); %estimate uncertainty in total flux
        else
            Q = 0;
            sigma_Q = 0;
        end
        Q_all{i}(ind_Q_NaN(j))=Q;
        sigma_Q_all{i}(ind_Q_NaN(j))=sigma_Q;
    end
end

% remove (redo) profiles with repeated Wenglor heights (Jeri only)
for i = 1:N_Sites
    N_windows = length(q_all{i});
    ind_norepeat = [];
    for j = 1:N_windows
        [n, bin] = histc(zW_all{i}{j}, unique(zW_all{i}{j}));
        if max(n)==1
            ind_norepeat = [ind_norepeat; j];
        end
    end
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_norepeat);']);
    end
end

% remove outlier points for Rancho Guadalupe (March 23, 17:29-18:29)
ind_outlier = 26:32;
ind_notoutlier = setdiff(1:length(Q_all{2}),ind_outlier);
i = 2;
for j = 1:N_variables
    eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_notoutlier);']);
end

%% COMPUTE TAUEX, TAURATIO, ETATAU, ETATAURATIO, BASED ON THRESHOLD BY SITE
tauex_all = cell(N_Sites,1);
tauratio_all = cell(N_Sites,1);
etatau_active_all = cell(N_Sites,1);
etatau_inactive_all = cell(N_Sites,1);
etatauratio_all = cell(N_Sites,1);
ustth_TFEM_all = cell(N_Sites,1);
for i=1:N_Sites
    tauex_all{i} = tauRe_all{i}-tauit(i);
    tauratio_all{i} = tauex_all{i}./tauRe_all{i};
    etatau_active_all{i} = eta_active_all{i}.*tauRe_all{i};
    etatau_inactive_all{i} = eta_inactive_all{i}.*tauRe_all{i};
    etatauratio_all{i} = eta_all{i}./tauratio_all{i};
    ustth_TFEM_all{i} = sqrt(tauth_TFEM_all{i}/rho_a);
end

%% COMPUTE ZQNORM BASED ON GRAIN SIZE BY SITE
zqnorm_all = cell(N_Sites,1);
zq_bar = zeros(N_Sites,1);
zqnorm_bar = zeros(N_Sites,1);
for i=1:N_Sites
    zq_bar(i) = mean(zq_all{i}(zq_all{i}>0));
    zqnorm_all{i} = 1000*zq_all{i}./d50_all{i};
    zqnorm_bar(i) = mean(zqnorm_all{i}(zqnorm_all{i}>0));
end

%% COMPUTE Q_NORM
Qnorm_all = cell(N_Sites,1);
for i=1:N_Sites
    Qnorm_all{i} = Q_all{i}./((1/g)*sqrt(tauit(i)/rho_a).*1e3*tauex_all{i});
end

%% CREATE BINS
%set frequency bins
fQ_continuous = 0.9; %minimum mean fQ for tauex bin to be called "continuous"
fQ_transport = 0.1; %minimum mean fQ for tauex bin to be called "transport"
fQ_bins_min = [0 0.3 0.6 0.9];
fQ_bins_max = [0.3 0.6 0.9 1];
fQ_bins_legend_all = {'f_{Q} = 0-0.3','f_{Q} = 0.3-0.6', 'f_{Q} = 0.6-0.9', 'f_{Q} = 0.9-1'};
Markers_fQ_bins = {'rx','bv','g^','ko'};
fQ_bins_mid = mean([fQ_bins_min; fQ_bins_max]);
N_fQ_bins = length(fQ_bins_mid);

%set alternate fQ bins
fQ_altbins_min = 0:0.05:0.95;
fQ_altbins_max = 0.05:0.05:1;
fQ_altbins_mid = mean([fQ_altbins_min; fQ_altbins_max]);
N_fQ_altbins = length(fQ_altbins_mid);
fQ_TFEM_fit_min = [0.65,0.7,0.1]; %mininum for fitting to obtain impact/fluid thresholds
fQ_TFEM_fit_max = [0.9,0.85,0.9]; %maximum for fitting to obtain impact/fluid thresholds
% fQ_TFEM_fit_min = [0.65,0.75,0.15]; %mininum for fitting to obtain impact/fluid thresholds
% fQ_TFEM_fit_max = [0.85,0.95,0.95]; %maximum for fitting to obtain impact/fluid thresholds

%create u* bins
ust_bins_min = 0:.01:0.56;
ust_bins_max = 0.01:.01:0.57;
ust_bins_mid = mean([ust_bins_min;ust_bins_max]);
N_ust_bins = length(ust_bins_mid);

%tau_bins
%tau_bins_min = 0.00:0.02:0.38; %medium bins
%tau_bins_max = 0.02:0.02:0.40; %medium bins
tau_bins_min = 0.00:0.01:0.39; %small bins
tau_bins_max = 0.01:0.01:0.40; %small bins
tau_bins_mid = mean([tau_bins_min;tau_bins_max]);
N_tau_bins = length(tau_bins_mid);

%tauex_bins
%tauex_bins_min = -0.10:0.02:0.28; %medium bins
%tauex_bins_max = -0.08:0.02:0.30; %medium bins
tauex_bins_min = -0.10:0.01:0.29; %small bins
tauex_bins_max = -0.09:0.01:0.30; %small bins
tauex_bins_mid = mean([tauex_bins_min;tauex_bins_max]);
N_tauex_bins = length(tauex_bins_mid);

%eta bins
eta_bins_min = 0:.025:0.475;
eta_bins_max = 0.025:.025:0.5;
eta_bins_mid = mean([eta_bins_min;eta_bins_max]);
N_eta_bins = length(eta_bins_mid);

%etatau_bins
etatau_bins_min = 0:0.01:0.24;
etatau_bins_max = 0.01:0.01:0.25;
etatau_bins_mid = mean([etatau_bins_min;etatau_bins_max]);
N_etatau_bins = length(etatau_bins_mid);

%etatauratio_bins
etatauratio_bins_min = [0 0.7 0.9 1.1 1.3];
etatauratio_bins_max = [0.7 0.9 1.1 1.3 2.5];
etatauratio_bins_mid = mean([etatau_bins_min;etatau_bins_max]);
etatauratio_bins_legend_all = {...
    '\eta/(\tau_{ex}/\tau) = 0-0.7',...
    '\eta/(\tau_{ex}/\tau) = 0.7-0.9',...
    '\eta/(\tau_{ex}/\tau) = 0.9-1.1',...
    '\eta/(\tau_{ex}/\tau) = 1.1-1.3',...
    '\eta/(\tau_{ex}/\tau) = 1.3-2.5'};
Markers_etatauratio_bins = {'rx','mv','bd','g^','ko'};
N_etatauratio_bins = length(etatauratio_bins_min);


%% FOR EACH SITE, PERFORM BINNING BY fQ
fQ_bins_legend = cell(N_Sites,1); %get legend values for each site
tau_fQ_bin_values = cell(N_Sites,1); %tau into frequency bins
tauex_fQ_bin_values = cell(N_Sites,1); %tauex into frequency bins
Q_fQ_bin_values = cell(N_Sites,1); %fluxes into frequency bins
eta_inactive_fQ_bin_values = cell(N_Sites,1); %etas into frequency bins
eta_active_fQ_bin_values = cell(N_Sites,1); %etas into frequency bins
zs_fQ_bin_values = cell(N_Sites,1); %zs into frequency bins
uratio_fQ_bin_values = cell(N_Sites,1); %ubar/ustd into frequency bins
zq_fQ_bin_values = cell(N_Sites,1); %saltation characteristic height

%apply values to bins
for i=1:N_Sites
    
    %initialize bins for each site
    fQ_bins_legend{i} = cell(N_Sites,1);
    tau_fQ_bin_values{i} = cell(N_fQ_bins,1);
    tauex_fQ_bin_values{i} = cell(N_fQ_bins,1);
    Q_fQ_bin_values{i} = cell(N_fQ_bins,1);
    eta_inactive_fQ_bin_values{i} = cell(N_fQ_bins,1);
    eta_active_fQ_bin_values{i} = cell(N_fQ_bins,1);
    zs_fQ_bin_values{i} = cell(N_fQ_bins,1);
    uratio_fQ_bin_values{i} = cell(N_fQ_bins,1);
    zq_fQ_bin_values{i} = cell(N_fQ_bins,1);
    
    %go through frequency bins, get values
    for j=1:N_fQ_bins
        bin_ind = find(fQ1_all{i}>=fQ_bins_min(j)&fQ1_all{i}<=fQ_bins_max(j)); %use 1 second frequencies for binning, soft upper limit
        fQ_bins_legend{i}{j} = strcat(fQ_bins_legend_all{j},'; N = ',int2str(length(bin_ind)));
        if ~isempty(bin_ind)
            tau_fQ_bin_values{i}{j} = tauRe_all{i}(bin_ind);
            tauex_fQ_bin_values{i}{j} = tauex_all{i}(bin_ind);
            Q_fQ_bin_values{i}{j} = Q_all{i}(bin_ind);
            eta_inactive_fQ_bin_values{i}{j} = eta_inactive_all{i}(bin_ind); %use eta based on inactive method
            eta_active_fQ_bin_values{i}{j} = eta_active_all{i}(bin_ind); %use eta based on active method
            zs_fQ_bin_values{i}{j} = zs_all{i}(bin_ind);
            uratio_fQ_bin_values{i}{j} = u_bar_all{i}(bin_ind)./u_std_all{i}(bin_ind);
            zq_fQ_bin_values{i}{j} = zq_all{i}(bin_ind);
        end
    end
end

%% FOR EACH SITE, PERFORM BINNING BY ETATAURATIO
etatauratio_bins_legend = cell(N_Sites,1); %get legend values for each site
tau_etatauratio_bin_values = cell(N_Sites,1); %tau into etatauratio bins
tauex_etatauratio_bin_values = cell(N_Sites,1); %tauex into etatauratio bins
Q_etatauratio_bin_values = cell(N_Sites,1); %fluxes into etatauratio bins
fQ_etatauratio_bin_values = cell(N_Sites,1); %fQs into etatauratio bins
zs_etatauratio_bin_values = cell(N_Sites,1); %zs into etatauratio bins
ubar_etatauratio_bin_values = cell(N_Sites,1); %ubar into etatauratio bins
ustd_etatauratio_bin_values{i} = cell(N_Sites,1); %ustd into etatauratio bins
uratio_etatauratio_bin_values{i} = cell(N_Sites,1); %uubar/ustd into etatauratio bins
zq_etatauratio_bin_values = cell(N_Sites,1); %saltation characteristic height

%apply values to bins
for i=1:N_Sites
    
    %initialize bins for each site
    etatauratio_bins_legend{i} = cell(N_Sites,1);
    tau_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    tauex_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    Q_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    fQ_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    zs_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    ubar_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    ustd_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    uratio_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    zq_etatauratio_bin_values{i} = cell(N_etatauratio_bins,1);
    
    %go through etatauratio bins, get values
    for j=1:N_etatauratio_bins
        bin_ind = find(etatauratio_all{i}>etatauratio_bins_min(j)&etatauratio_all{i}<=etatauratio_bins_max(j)); %bin by etatauratio
        etatauratio_bins_legend{i}{j} = strcat(etatauratio_bins_legend_all{j},'; N = ',int2str(length(bin_ind)));
        if ~isempty(bin_ind)
            tau_etatauratio_bin_values{i}{j} = tauRe_all{i}(bin_ind);
            tauex_etatauratio_bin_values{i}{j} = tauex_all{i}(bin_ind);
            Q_etatauratio_bin_values{i}{j} = Q_all{i}(bin_ind);
            fQ_etatauratio_bin_values{i}{j} = eta_inactive_all{i}(bin_ind);
            zs_etatauratio_bin_values{i}{j} = zs_all{i}(bin_ind);
            ubar_etatauratio_bin_values{i}{j} = u_bar_all{i}(bin_ind);
            ustd_etatauratio_bin_values{i}{j} = u_std_all{i}(bin_ind);
            uratio_etatauratio_bin_values{i}{j} = u_bar_all{i}(bin_ind)./u_std_all{i}(bin_ind);
            zq_etatauratio_bin_values{i}{j} = zq_all{i}(bin_ind);
        end
    end
end

%% FOR EACH SITE, PERFORM BINNING BY ETA

%separate fQ's into eta bins - inactive definition
fQ_eta_inactive_bin_values = cell(N_Sites,1); %fQ's into eta bins - inactive definition
fQ_eta_inactive_bin_avg = cell(N_Sites,1); %get average fQ for eta bin
fQ_eta_inactive_bin_std = cell(N_Sites,1); %get std fQ for eta bin
fQ_eta_inactive_bin_SE = cell(N_Sites,1); %get SE fQ for eta bin

%separate fQ's into eta bins - active definition
fQ_eta_active_bin_values = cell(N_Sites,1); %fQ's into eta bins - active definition
fQ_eta_active_bin_avg = cell(N_Sites,1); %get average fQ for eta bin
fQ_eta_active_bin_std = cell(N_Sites,1); %get std fQ for eta bin
fQ_eta_active_bin_SE = cell(N_Sites,1); %get SE fQ for eta bin

for i = 1:N_Sites
    %initialize fQ's into eta bins - fQ definition
    fQ_eta_inactive_bin_values{i} = cell(N_eta_bins,1); 
    fQ_eta_inactive_bin_avg{i} = zeros(N_eta_bins,1)*NaN;
    fQ_eta_inactive_bin_std{i} = zeros(N_eta_bins,1)*NaN;
    fQ_eta_inactive_bin_SE{i} = zeros(N_eta_bins,1)*NaN;
    
    %initialize fQ's into eta bins - zs definition
    fQ_eta_active_bin_values{i} = cell(N_eta_bins,1); 
    fQ_eta_active_bin_avg{i} = zeros(N_eta_bins,1)*NaN;
    fQ_eta_active_bin_std{i} = zeros(N_eta_bins,1)*NaN;
    fQ_eta_active_bin_SE{i} = zeros(N_eta_bins,1)*NaN;
    
    %get fQs for eta bins
    for j = 1:N_eta_bins
        %using inactive definition of eta
        bin_ind = find(eta_inactive_all{i}>=eta_bins_min(j)&eta_inactive_all{i}<=eta_bins_max(j));
        fQ_eta_inactive_bin_values{i}{j} = fQ1_all{i}(bin_ind);
        fQ_eta_inactive_bin_values{i}{j} = fQ_eta_inactive_bin_values{i}{j}(~isnan(fQ_eta_inactive_bin_values{i}{j}));
        fQ_eta_inactive_bin_avg{i}(j) = mean(fQ_eta_inactive_bin_values{i}{j});
        fQ_eta_inactive_bin_std{i}(j) = std(fQ_eta_inactive_bin_values{i}{j});
        fQ_eta_inactive_bin_SE{i}(j) = std(fQ_eta_inactive_bin_values{i}{j})./sqrt(length(fQ_eta_inactive_bin_values{i}{j}));

        %using active definition of eta
        bin_ind = find(eta_active_all{i}>=eta_bins_min(j)&eta_active_all{i}<=eta_bins_max(j));
        fQ_eta_active_bin_values{i}{j} = fQ1_all{i}(bin_ind);
        fQ_eta_active_bin_values{i}{j} = fQ_eta_active_bin_values{i}{j}(~isnan(fQ_eta_active_bin_values{i}{j}));
        fQ_eta_active_bin_avg{i}(j) = mean(fQ_eta_active_bin_values{i}{j});
        fQ_eta_active_bin_std{i}(j) = std(fQ_eta_active_bin_values{i}{j});
        fQ_eta_active_bin_SE{i}(j) = std(fQ_eta_active_bin_values{i}{j})./sqrt(length(fQ_eta_active_bin_values{i}{j}));
    end
end

%% FOR EACH SITE, PERFORM BINNING BY ETATAU

%separate Q's into etatau bins - inactive method
Q_etatau_inactive_bin_values = cell(N_Sites,1); %Q's into etatau bins - inactive definition
Q_etatau_inactive_bin_avg = cell(N_Sites,1); %get average Q for etatau bin
Q_etatau_inactive_bin_std = cell(N_Sites,1); %get std Q for etatau bin
Q_etatau_inactive_bin_SE = cell(N_Sites,1); %get SE Q for etatau bin

%separate Q's into etatau bins - active method
Q_etatau_active_bin_values = cell(N_Sites,1); %Q's into etatau bins - active definition
Q_etatau_active_bin_avg = cell(N_Sites,1); %get average Q for etatau bin
Q_etatau_active_bin_std = cell(N_Sites,1); %get std Q for etatau bin
Q_etatau_active_bin_SE = cell(N_Sites,1); %get SE Q for etatau bin

for i = 1:N_Sites
    %initialize Q's into etatau bins - inactive method
    Q_etatau_inactive_bin_values{i} = cell(N_etatau_bins,1); 
    Q_etatau_inactive_bin_avg{i} = zeros(N_etatau_bins,1)*NaN;
    Q_etatau_inactive_bin_std{i} = zeros(N_etatau_bins,1)*NaN;
    Q_etatau_inactive_bin_SE{i} = zeros(N_etatau_bins,1)*NaN;
    
    %initialize Q's into etatau bins - active method
    Q_etatau_active_bin_values{i} = cell(N_etatau_bins,1); 
    Q_etatau_active_bin_avg{i} = zeros(N_etatau_bins,1)*NaN;
    Q_etatau_active_bin_std{i} = zeros(N_etatau_bins,1)*NaN;
    Q_etatau_active_bin_SE{i} = zeros(N_etatau_bins,1)*NaN;
    
    %get Qs for etatau bins
    for j = 1:N_etatau_bins
        %using inactive definition of eta
        bin_ind = find(etatau_inactive_all{i}>=etatau_bins_min(j)&etatau_inactive_all{i}<=etatau_bins_max(j));
        Q_etatau_inactive_bin_values{i}{j} = Q_all{i}(bin_ind);
        ind_good = find(~isnan(Q_etatau_inactive_bin_values{i}{j}));
        Q_etatau_inactive_bin_values{i}{j} = Q_etatau_inactive_bin_values{i}{j}(~isnan(Q_etatau_inactive_bin_values{i}{j}));
        if ~isempty(ind_good)
            sigma_Q_etatau_values = sigma_Q_all{i}(bin_ind);
            sigma_Q_etatau_values = sigma_Q_etatau_values(ind_good);
            Q_etatau_inactive_bin_avg{i}(j) = mean(Q_etatau_inactive_bin_values{i}{j});
            if length(Q_etatau_inactive_bin_values{i}{j})==1 %use sample uncertainty if there is only one sample
                Q_etatau_inactive_bin_std{i}(j) = sigma_Q_etatau_values; 
                Q_etatau_inactive_bin_SE{i}(j) = sigma_Q_etatau_values;
            else %otherwise compute standard error
                Q_etatau_inactive_bin_std{i}(j) = std(Q_etatau_inactive_bin_values{i}{j});
                Q_etatau_inactive_bin_SE{i}(j) = std(Q_etatau_inactive_bin_values{i}{j})./sqrt(length(Q_etatau_inactive_bin_values{i}{j}));
            end
        end
            
        %using active definition of eta
        bin_ind = find(etatau_active_all{i}>=etatau_bins_min(j)&etatau_active_all{i}<=etatau_bins_max(j));
        Q_etatau_active_bin_values{i}{j} = Q_all{i}(bin_ind);
        ind_good = find(~isnan(Q_etatau_active_bin_values{i}{j}));
        Q_etatau_active_bin_values{i}{j} = Q_etatau_active_bin_values{i}{j}(~isnan(Q_etatau_active_bin_values{i}{j}));
        if ~isempty(ind_good)
            sigma_Q_etatau_values = sigma_Q_all{i}(bin_ind);
            sigma_Q_etatau_values = sigma_Q_etatau_values(ind_good);
            Q_etatau_active_bin_avg{i}(j) = mean(Q_etatau_active_bin_values{i}{j});
            if length(Q_etatau_active_bin_values{i}{j})==1 %use sample uncertainty if there is only one sample
                Q_etatau_active_bin_std{i}(j) = sigma_Q_etatau_values; 
                Q_etatau_active_bin_SE{i}(j) = sigma_Q_etatau_values;
            else %otherwise compute standard error
                Q_etatau_active_bin_std{i}(j) = std(Q_etatau_active_bin_values{i}{j});
                Q_etatau_active_bin_SE{i}(j) = std(Q_etatau_active_bin_values{i}{j})./sqrt(length(Q_etatau_active_bin_values{i}{j}));
            end
        end
    end
end

%% FOR EACH SITE, PERFORM BINNING BY USTAR

%separate fQs into ust bins
fQ_ust_bin_values = cell(N_Sites,1); %fQs into ust bins
fQ_ust_bin_avg = cell(N_Sites,1); %get average fQ for ust bins
fQ_ust_bin_std = cell(N_Sites,1); %get std fQ for ust bins
fQ_ust_bin_SE = cell(N_Sites,1); %get SE fQ for ust bins

%separate zq's into ust bins
zq_ust_bin_values = cell(N_Sites,1); %zq into ust bins
zq_ust_bin_avg = cell(N_Sites,1); %get average zq for ust bins
zq_ust_bin_std = cell(N_Sites,1); %get std zq for ust bins
zq_ust_bin_SE = cell(N_Sites,1); %get SE zq for ust bins

%separate zqnorm's (zq/d50) into ust bins
zqnorm_ust_bin_values = cell(N_Sites,1); %zq/d50s into ust bins
zqnorm_ust_bin_avg = cell(N_Sites,1); %get average zq/d50 for ust bins
zqnorm_ust_bin_std = cell(N_Sites,1); %get std zq/d50 for ust bins
zqnorm_ust_bin_SE = cell(N_Sites,1); %get SE zq/d50 for ust bins

%separate Qs into ust bins
Q_ust_bin_values = cell(N_Sites,1); %Qs into ust bins
Q_ust_bin_avg = cell(N_Sites,1); %get average Q for ust bins
Q_ust_bin_std = cell(N_Sites,1); %get std Q for ust bins
Q_ust_bin_SE = cell(N_Sites,1); %get SE Q for ust bins

%separate Q_norms into ust bins
Qnorm_ust_bin_values = cell(N_Sites,1); %Qs into ust bins
Qnorm_ust_bin_avg = cell(N_Sites,1); %get average Q for ust bins
Qnorm_ust_bin_std = cell(N_Sites,1); %get std Q for ust bins
Qnorm_ust_bin_SE = cell(N_Sites,1); %get SE Q for ust bins

%go through all sites
for i = 1:N_Sites

    %initialize fQ's into ust bins
    fQ_ust_bin_values{i} = cell(N_ust_bins,1); 
    fQ_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    fQ_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    fQ_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %initialize zq's into ust bins
    zq_ust_bin_values{i} = cell(N_ust_bins,1); 
    zq_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    zq_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    zq_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %initialize zqnorm's into ust bins
    zqnorm_ust_bin_values{i} = cell(N_ust_bins,1); 
    zqnorm_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    zqnorm_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    zqnorm_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %initialize Q's into ust bins
    Q_ust_bin_values{i} = cell(N_ust_bins,1); 
    Q_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    Q_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    Q_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %separate Q_norms into ust bins
    Qnorm_ust_bin_values{i} = cell(N_ust_bins,1); 
    Qnorm_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    Qnorm_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    Qnorm_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %go through each ust bin
    for j = 1:N_ust_bins
        %get indices of ust bins
        bin_ind = find(ustRe_all{i}>=ust_bins_min(j)&ustRe_all{i}<=ust_bins_max(j));

        %get fQs for ust bins
        fQ_ust_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_ust_bin_values{i}{j} = fQ_ust_bin_values{i}{j}(~isnan(fQ_ust_bin_values{i}{j}));
        fQ_ust_bin_avg{i}(j) = mean(fQ_ust_bin_values{i}{j});
        fQ_ust_bin_std{i}(j) = std(fQ_ust_bin_values{i}{j});
        fQ_ust_bin_SE{i}(j) = std(fQ_ust_bin_values{i}{j})./sqrt(length(fQ_ust_bin_values{i}{j}));
        
        %get zqs for ust bins
        zq_ust_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_ust_bin_values{i}{j} = zq_ust_bin_values{i}{j};
        ind_good = find(~isnan(zq_ust_bin_values{i}{j}));
        zq_ust_bin_values{i}{j} = zq_ust_bin_values{i}{j}(ind_good);
        if ~isempty(ind_good)
            zq_ust_bin_avg{i}(j) = mean(zq_ust_bin_values{i}{j});
            sigma_zq_ust_values = sigma_zq_all{i}(bin_ind);
            sigma_zq_ust_values = sigma_zq_ust_values(ind_good);
            if length(zq_ust_bin_values{i}{j})==1 %use sample uncertainty if there is only one sample
                zq_ust_bin_std{i}(j) = sigma_zq_ust_values; 
                zq_ust_bin_SE{i}(j) = sigma_zq_ust_values;
            else %otherwise compute standard error
                zq_ust_bin_std{i}(j) = std(zq_ust_bin_values{i}{j});
                zq_ust_bin_SE{i}(j) = std(zq_ust_bin_values{i}{j})./sqrt(length(zq_ust_bin_values{i}{j}));
            end
        end
            
        %get zqnorms for ust bins
        zqnorm_ust_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_ust_bin_values{i}{j} = zqnorm_ust_bin_values{i}{j}(~isnan(zqnorm_ust_bin_values{i}{j}));
        zqnorm_ust_bin_avg{i}(j) = mean(zqnorm_ust_bin_values{i}{j});
        zqnorm_ust_bin_std{i}(j) = std(zqnorm_ust_bin_values{i}{j});
        zqnorm_ust_bin_SE{i}(j) = std(zqnorm_ust_bin_values{i}{j})./sqrt(length(zqnorm_ust_bin_values{i}{j}));
        
        %get Qs for ust bins
        Q_ust_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_ust_bin_values{i}{j} = Q_ust_bin_values{i}{j}(~isnan(Q_ust_bin_values{i}{j}));
        Q_ust_bin_avg{i}(j) = mean(Q_ust_bin_values{i}{j});
        Q_ust_bin_std{i}(j) = std(Q_ust_bin_values{i}{j});
        Q_ust_bin_SE{i}(j) = std(Q_ust_bin_values{i}{j})./sqrt(length(Q_ust_bin_values{i}{j}));
        
        %get Qnorms for ust bins
        Qnorm_ust_bin_values{i}{j} = Qnorm_all{i}(bin_ind);
        Qnorm_ust_bin_values{i}{j} = Qnorm_ust_bin_values{i}{j}(~isnan(Qnorm_ust_bin_values{i}{j}));
        Qnorm_ust_bin_avg{i}(j) = mean(Qnorm_ust_bin_values{i}{j});
        Qnorm_ust_bin_std{i}(j) = std(Qnorm_ust_bin_values{i}{j});
        Qnorm_ust_bin_SE{i}(j) = std(Qnorm_ust_bin_values{i}{j})./sqrt(length(Qnorm_ust_bin_values{i}{j}));
        
    end
end


%% FOR EACH SITE, PERFORM BINNING BY TAU
%separate tauratios into tau bins
tauratio_tau_bin_values = cell(N_Sites,1); %tauratios into tau bins
tauratio_tau_bin_avg = cell(N_Sites,1); %get average tauratio for tau bins
tauratio_tau_bin_std = cell(N_Sites,1); %get std tauratio for tau bins
tauratio_tau_bin_SE = cell(N_Sites,1); %get SE tauratio for tau bins

%separate fQs into tau bins
fQ_tau_bin_values = cell(N_Sites,1); %fQs into tau bins
fQ_tau_bin_avg = cell(N_Sites,1); %get average fQ for tau bins
fQ_tau_bin_std = cell(N_Sites,1); %get std fQ for tau bins
fQ_tau_bin_SE = cell(N_Sites,1); %get SE fQ for tau bins

%separate Qs into tau bins
Q_tau_bin_values = cell(N_Sites,1); %Qs into tau bins
Q_tau_bin_avg = cell(N_Sites,1); %get average Q for tau bins
Q_tau_bin_std = cell(N_Sites,1); %get std Q for tau bins
Q_tau_bin_SE = cell(N_Sites,1); %get SE Q for tau bins

for i = 1:N_Sites
    %initialize tauratio's into tau bins
    tauratio_tau_bin_values{i} = cell(N_tau_bins,1); 
    tauratio_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    tauratio_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    tauratio_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;
    
    %initialize fQ's into tau bins
    fQ_tau_bin_values{i} = cell(N_tau_bins,1); 
    fQ_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    fQ_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    fQ_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;
    
    %initialize Q's into tau bins
    Q_tau_bin_values{i} = cell(N_tau_bins,1); 
    Q_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    Q_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    Q_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;
    
    %go through each tau bin
    for j = 1:N_tau_bins
        %get indices of tau bins
        bin_ind = find(tauRe_all{i}>=tau_bins_min(j)&tauRe_all{i}<=tau_bins_max(j));

        %get tauratios for tau bins
        tauratio_tau_bin_values{i}{j} = tauratio_all{i}(bin_ind);
        tauratio_tau_bin_values{i}{j} = tauratio_tau_bin_values{i}{j}(~isnan(tauratio_tau_bin_values{i}{j}));
        tauratio_tau_bin_avg{i}(j) = mean(tauratio_tau_bin_values{i}{j});
        tauratio_tau_bin_std{i}(j) = std(tauratio_tau_bin_values{i}{j});
        tauratio_tau_bin_SE{i}(j) = std(tauratio_tau_bin_values{i}{j})./sqrt(length(tauratio_tau_bin_values{i}{j}));
        
        %get fQs for tau bins
        fQ_tau_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_tau_bin_values{i}{j} = fQ_tau_bin_values{i}{j}(~isnan(fQ_tau_bin_values{i}{j}));
        fQ_tau_bin_avg{i}(j) = mean(fQ_tau_bin_values{i}{j});
        fQ_tau_bin_std{i}(j) = std(fQ_tau_bin_values{i}{j});
        fQ_tau_bin_SE{i}(j) = std(fQ_tau_bin_values{i}{j})./sqrt(length(fQ_tau_bin_values{i}{j}));
        
        %get Qs for tau bins
        Q_tau_bin_values{i}{j} = Q_all{i}(bin_ind);
        ind_good = find(~isnan(Q_tau_bin_values{i}{j}));
        Q_tau_bin_values{i}{j} = Q_tau_bin_values{i}{j}(~isnan(Q_tau_bin_values{i}{j}));
        Q_tau_bin_avg{i}(j) = mean(Q_tau_bin_values{i}{j});
        sigma_Q_tau_values = sigma_Q_all{i}(bin_ind);
        sigma_Q_tau_values = sigma_Q_tau_values(ind_good);
        if length(Q_tau_bin_values{i}{j})==1 %use sample uncertainty if there is only one sample
            Q_tau_bin_std{i}(j) = sigma_Q_tau_values; 
            Q_tau_bin_SE{i}(j) = sigma_Q_tau_values;
        else %otherwise compute standard error
            Q_tau_bin_std{i}(j) = std(Q_tau_bin_values{i}{j});
            Q_tau_bin_SE{i}(j) = std(Q_tau_bin_values{i}{j})./sqrt(length(Q_tau_bin_values{i}{j}));
        end
    end
end

%% FOR EACH SITE, PERFORM BINNING BY TAUEX
%separate tauratios into tau_ex bins
tauratio_tauex_bin_values = cell(N_Sites,1); %tauratios into tau_ex bins
tauratio_tauex_bin_avg = cell(N_Sites,1); %get average tauratio for tau_ex bins
tauratio_tauex_bin_std = cell(N_Sites,1); %get std tauratio for tau_ex bins
tauratio_tauex_bin_SE = cell(N_Sites,1); %get SE tauratio for tau_ex bins

%separate fQs into tau_ex bins
fQ_tauex_bin_values = cell(N_Sites,1); %fQs into tau_ex bins
fQ_tauex_bin_avg = cell(N_Sites,1); %get average fQ for tau_ex bins
fQ_tauex_bin_std = cell(N_Sites,1); %get std fQ for tau_ex bins
fQ_tauex_bin_SE = cell(N_Sites,1); %get SE fQ for tau_ex bins

%separate zq's into tau_ex bins
zq_tauex_bin_values = cell(N_Sites,1); %zq into tau_ex bins
zq_tauex_bin_avg = cell(N_Sites,1); %get average zq for tau_ex bins
zq_tauex_bin_std = cell(N_Sites,1); %get std zq for tau_ex bins
zq_tauex_bin_SE = cell(N_Sites,1); %get SE zq for tau_ex bins

%separate zqnorm's (zq/d50) into tau_ex bins
zqnorm_tauex_bin_values = cell(N_Sites,1); %zq/d50s into tau_ex bins
zqnorm_tauex_bin_avg = cell(N_Sites,1); %get average zq/d50 for tau_ex bins
zqnorm_tauex_bin_std = cell(N_Sites,1); %get std zq/d50 for tau_ex bins
zqnorm_tauex_bin_SE = cell(N_Sites,1); %get SE zq/d50 for tau_ex bins

%separate Qs into tau_ex bins
Q_tauex_bin_values = cell(N_Sites,1); %Qs into tau_ex bins
Q_tauex_bin_avg = cell(N_Sites,1); %get average Q for tau_ex bins
Q_tauex_bin_std = cell(N_Sites,1); %get std Q for tau_ex bins
Q_tauex_bin_SE = cell(N_Sites,1); %get SE Q for tau_ex bins

%separate Q_norms into tau_ex bins
Qnorm_tauex_bin_values = cell(N_Sites,1); %Qs into tau_ex bins
Qnorm_tauex_bin_avg = cell(N_Sites,1); %get average Q for tau_ex bins
Qnorm_tauex_bin_std = cell(N_Sites,1); %get std Q for tau_ex bins
Qnorm_tauex_bin_SE = cell(N_Sites,1); %get SE Q for tau_ex bins

%separate etas into tau_ex bins
eta_tauex_bin_values = cell(N_Sites,1); %etas into tau_ex bins
eta_tauex_bin_avg = cell(N_Sites,1); %get average eta for tau_ex bins
eta_tauex_bin_std = cell(N_Sites,1); %get std eta for tau_ex bins
eta_tauex_bin_SE = cell(N_Sites,1); %get SE eta for tau_ex bins

for i = 1:N_Sites
    %initialize tauratio's into tauex bins
    tauratio_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    tauratio_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    tauratio_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    tauratio_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize fQ's into tauex bins
    fQ_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    fQ_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    fQ_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    fQ_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;

    %initialize zq's into tauex bins
    zq_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    zq_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    zq_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    zq_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize zqnorm's into tauex bins
    zqnorm_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    zqnorm_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    zqnorm_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    zqnorm_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize Q's into tauex bins
    Q_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    Q_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    Q_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    Q_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %separate Q_norms into tau_ex bins
    Qnorm_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    Qnorm_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    Qnorm_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    Qnorm_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize eta's into tauex bins
    eta_tauex_bin_values{i} = cell(N_tauex_bins,1); 
    eta_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    eta_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    eta_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %go through each tauex bin
    for j = 1:N_tauex_bins
        %get indices of tauex bins
        bin_ind = find(tauex_all{i}>=tauex_bins_min(j)&tauex_all{i}<=tauex_bins_max(j));

        %get tauratios for tauex bins
        tauratio_tauex_bin_values{i}{j} = tauratio_all{i}(bin_ind);
        tauratio_tauex_bin_values{i}{j} = tauratio_tauex_bin_values{i}{j}(~isnan(tauratio_tauex_bin_values{i}{j}));
        tauratio_tauex_bin_avg{i}(j) = mean(tauratio_tauex_bin_values{i}{j});
        tauratio_tauex_bin_std{i}(j) = std(tauratio_tauex_bin_values{i}{j});
        tauratio_tauex_bin_SE{i}(j) = std(tauratio_tauex_bin_values{i}{j})./sqrt(length(tauratio_tauex_bin_values{i}{j}));
        
        %get fQs for tauex bins
        fQ_tauex_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_tauex_bin_values{i}{j} = fQ_tauex_bin_values{i}{j}(~isnan(fQ_tauex_bin_values{i}{j}));
        fQ_tauex_bin_avg{i}(j) = mean(fQ_tauex_bin_values{i}{j});
        fQ_tauex_bin_std{i}(j) = std(fQ_tauex_bin_values{i}{j});
        fQ_tauex_bin_SE{i}(j) = std(fQ_tauex_bin_values{i}{j})./sqrt(length(fQ_tauex_bin_values{i}{j}));

        %get zqs for tauex bins
        zq_tauex_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_tauex_bin_values{i}{j} = zq_tauex_bin_values{i}{j}(~isnan(zq_tauex_bin_values{i}{j}));
        zq_tauex_bin_avg{i}(j) = mean(zq_tauex_bin_values{i}{j});
        zq_tauex_bin_std{i}(j) = std(zq_tauex_bin_values{i}{j});
        zq_tauex_bin_SE{i}(j) = std(zq_tauex_bin_values{i}{j})./sqrt(length(zq_tauex_bin_values{i}{j}));

        %get zqnorms for tauex bins
        zqnorm_tauex_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_tauex_bin_values{i}{j} = zqnorm_tauex_bin_values{i}{j}(~isnan(zqnorm_tauex_bin_values{i}{j}));
        zqnorm_tauex_bin_avg{i}(j) = mean(zqnorm_tauex_bin_values{i}{j});
        zqnorm_tauex_bin_std{i}(j) = std(zqnorm_tauex_bin_values{i}{j});
        zqnorm_tauex_bin_SE{i}(j) = std(zqnorm_tauex_bin_values{i}{j})./sqrt(length(zqnorm_tauex_bin_values{i}{j}));
        
        %get Qs for tauex bins
        Q_tauex_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_tauex_bin_values{i}{j} = Q_tauex_bin_values{i}{j}(~isnan(Q_tauex_bin_values{i}{j}));
        Q_tauex_bin_avg{i}(j) = mean(Q_tauex_bin_values{i}{j});
        Q_tauex_bin_std{i}(j) = std(Q_tauex_bin_values{i}{j});
        Q_tauex_bin_SE{i}(j) = std(Q_tauex_bin_values{i}{j})./sqrt(length(Q_tauex_bin_values{i}{j}));
        
        %get Qnorms for tauex bins
        Qnorm_tauex_bin_values{i}{j} = Qnorm_all{i}(bin_ind);
        Qnorm_tauex_bin_values{i}{j} = Qnorm_tauex_bin_values{i}{j}(~isnan(Qnorm_tauex_bin_values{i}{j}));
        Qnorm_tauex_bin_avg{i}(j) = mean(Qnorm_tauex_bin_values{i}{j});
        Qnorm_tauex_bin_std{i}(j) = std(Qnorm_tauex_bin_values{i}{j});
        Qnorm_tauex_bin_SE{i}(j) = std(Qnorm_tauex_bin_values{i}{j})./sqrt(length(Qnorm_tauex_bin_values{i}{j}));
        
        %get etas for tauex bins
        eta_tauex_bin_values{i}{j} = eta_all{i}(bin_ind);
        eta_tauex_bin_values{i}{j} = eta_tauex_bin_values{i}{j}(~isnan(eta_tauex_bin_values{i}{j}));
        eta_tauex_bin_avg{i}(j) = mean(eta_tauex_bin_values{i}{j});
        eta_tauex_bin_std{i}(j) = std(eta_tauex_bin_values{i}{j});
        eta_tauex_bin_SE{i}(j) = std(eta_tauex_bin_values{i}{j})./sqrt(length(eta_tauex_bin_values{i}{j}));
    end
end

%% FOR EACH SITE, PERFORM BINNING BY ALTERNATE (SMALLER) fQ BINS (more bins, for direct comparisons to frequency)

%separate tauth based on TFEM into frequency bins
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
        %get indices for bins
        bin_ind = find(fQ1_all{i}>=fQ_altbins_min(j)&fQ1_all{i}<=fQ_altbins_max(j)); %use 1 second frequencies for binning
                
        %tau_th
        tauth_fQ_altbin_values{i}{j} = tauth_TFEM_all{i}(bin_ind);
        tauth_fQ_altbin_values{i}{j} = tauth_fQ_altbin_values{i}{j}(~isnan(tauth_fQ_altbin_values{i}{j}));
        tauth_fQ_altbin_avg{i}(j) = mean(tauth_fQ_altbin_values{i}{j});
        tauth_fQ_altbin_std{i}(j) = std(tauth_fQ_altbin_values{i}{j});
        tauth_fQ_altbin_SE{i}(j) = std(tauth_fQ_altbin_values{i}{j})/sqrt(length(tauth_fQ_altbin_values{i}{j}));
            
        %fraction above fluid threshold
        f_uaboveft_values{i}{j} = f_uaboveft_all{i}(bin_ind);
        f_uaboveft_values{i}{j} = f_uaboveft_values{i}{j}(~isnan(f_uaboveft_values{i}{j}));
        f_uaboveft_avg{i}(j) = mean(f_uaboveft_values{i}{j});
        f_uaboveft_std{i}(j) = std(f_uaboveft_values{i}{j});
        f_uaboveft_SE{i}(j) = std(f_uaboveft_values{i}{j})/sqrt(length(f_uaboveft_values{i}{j}));
        
        %fraction below impact threshold
        f_ubelowit_values{i}{j} = f_ubelowit_all{i}(bin_ind);
        f_ubelowit_values{i}{j} = f_ubelowit_values{i}{j}(~isnan(f_ubelowit_values{i}{j}));
        f_ubelowit_avg{i}(j) = mean(f_ubelowit_values{i}{j});
        f_ubelowit_std{i}(j) = std(f_ubelowit_values{i}{j});
        f_ubelowit_SE{i}(j) = std(f_ubelowit_values{i}{j})/sqrt(length(f_ubelowit_values{i}{j}));
        
        %fraction in hysteresis range, from above
        f_ubetween_fromabove_values{i}{j} = f_ubetween_fromabove_all{i}(bin_ind);
        f_ubetween_fromabove_values{i}{j} = f_ubetween_fromabove_values{i}{j}(~isnan(f_ubetween_fromabove_values{i}{j}));
        f_ubetween_fromabove_avg{i}(j) = mean(f_ubetween_fromabove_values{i}{j});
        f_ubetween_fromabove_std{i}(j) = std(f_ubetween_fromabove_values{i}{j});
        f_ubetween_fromabove_SE{i}(j) = std(f_ubetween_fromabove_values{i}{j})/sqrt(length(f_ubetween_fromabove_values{i}{j}));
        
        %fraction in hysteresis range, from below
        f_ubetween_frombelow_values{i}{j} = f_ubetween_frombelow_all{i}(bin_ind);
        f_ubetween_frombelow_values{i}{j} = f_ubetween_frombelow_values{i}{j}(~isnan(f_ubetween_frombelow_values{i}{j}));
        f_ubetween_frombelow_avg{i}(j) = mean(f_ubetween_frombelow_values{i}{j});
        f_ubetween_frombelow_std{i}(j) = std(f_ubetween_frombelow_values{i}{j});
        f_ubetween_frombelow_SE{i}(j) = std(f_ubetween_frombelow_values{i}{j})/sqrt(length(f_ubetween_frombelow_values{i}{j}));
    end 
end

%% FOR EACH SITE, PERFORM FITS

%linear fit to flux versus tau during transport
tau_transport = cell(N_Sites,1);
sigma_tau_transport = cell(N_Sites,1);
Q_transport = cell(N_Sites,1);
sigma_Q_transport = cell(N_Sites,1);
CQ_pred_transport = zeros(N_Sites,1);
sigma_CQ_pred_transport = zeros(N_Sites,1);
tauit_pred_transport = zeros(N_Sites,1);
sigma_tauit_pred_transport = zeros(N_Sites,1);
Q_tau_RMSE_transport = zeros(N_Sites,1); %RMSE in fit
Q_tau_r2_transport = zeros(N_Sites,1); %r2 in fit
tau_conf_min_transport = cell(N_Sites,1);
tau_conf_max_transport = cell(N_Sites,1);
Q_conf_min_transport = cell(N_Sites,1);
Q_conf_max_transport = cell(N_Sites,1);
tauex_transport = cell(N_Sites,1);
Q_tauex_transport = cell(N_Sites,1);
%unbinned values
Q_tau_RMSE_transport_unbinned = zeros(N_Sites,1); %RMSE in fit - unbinned
Q_tau_r2_transport_unbinned = zeros(N_Sites,1); %r2 in fit - unbinned

%linear fit to flux versus tau during continuous transport
tau_continuous = cell(N_Sites,1);
sigma_tau_continuous = cell(N_Sites,1);
Q_continuous = cell(N_Sites,1);
sigma_Q_continuous = cell(N_Sites,1);
CQ_pred_continuous = zeros(N_Sites,1);
sigma_CQ_pred_continuous = zeros(N_Sites,1);
tauit_pred_continuous = zeros(N_Sites,1);
sigma_tauit_pred_continuous = zeros(N_Sites,1);
tau_conf_min_continuous = cell(N_Sites,1);
tau_conf_max_continuous = cell(N_Sites,1);
Q_conf_min_continuous = cell(N_Sites,1);
Q_conf_max_continuous = cell(N_Sites,1);
tauex_continuous = cell(N_Sites,1);
Q_tauex_continuous = cell(N_Sites,1);

%determination of flux versus tau with known threshold
CQ_pred_fixedtauit = zeros(N_Sites,1);
sigma_CQ_pred_fixedtauit = zeros(N_Sites,1);
fQ_TFEM_ind = cell(N_Sites,1);
tauft_TFEM = zeros(N_Sites,1);
tauit_TFEM = zeros(N_Sites,1);
tauth_fit_TFEM = cell(N_Sites,1);
sigma_tauth_fit_TFEM = cell(N_Sites,1);

%fit values for etatau_inactive
CQ_pred_etatau_inactive = zeros(N_Sites,1);
sigma_CQ_pred_etatau_inactive = zeros(N_Sites,1);
Q_etatau_RMSE_inactive = zeros(N_Sites,1); %RMSE in fit
Q_etatau_r2_inactive = zeros(N_Sites,1); %r2 in fit
%unbinned values
Q_etatau_RMSE_inactive_unbinned = zeros(N_Sites,1); %RMSE in fit - unbinned
Q_etatau_r2_inactive_unbinned = zeros(N_Sites,1); %r2 in fit - unbinned

%fit values for etatau_active
CQ_pred_etatau_active = zeros(N_Sites,1);
sigma_CQ_pred_etatau_active = zeros(N_Sites,1);
Q_etatau_RMSE_active = zeros(N_Sites,1); %RMSE in fit
Q_etatau_r2_active = zeros(N_Sites,1); %r2 in fit
%unbinned values
Q_etatau_RMSE_active_unbinned = zeros(N_Sites,1); %RMSE in fit - unbinned
Q_etatau_r2_active_unbinned = zeros(N_Sites,1); %r2 in fit - unbinned

%fit values for zq versus ust and tau_ex
linearfit_zq_ust = cell(N_Sites,1); %linear fit for zq versus ustar
intercept_zq_ust = cell(N_Sites,1); %intercept of fit
slope_zq_ust = cell(N_Sites,1); %slope of fit
sigmaslope_zq_ust = cell(N_Sites,1); %uncertainty in slope for zq versus ustar
linearfit_zq_tauex = cell(N_Sites,1); %linear fit for zq versus tau_ex

for i = 1:N_Sites
    %fit for transport bins
    tau_bin_ind = find(fQ_tau_bin_avg{i}>=fQ_transport); %get tau bins with transport
    tau_transport{i} = tau_bins_mid(tau_bin_ind)';
    sigma_tau_transport{i} = ones(size(tau_bin_ind))*mode(tau_bins_max-tau_bins_min);
    Q_transport{i} = Q_tau_bin_avg{i}(tau_bin_ind);
    sigma_Q_transport{i} = Q_tau_bin_SE{i}(tau_bin_ind);
    [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(tau_transport{i}, Q_transport{i}, sigma_tau_transport{i}, sigma_Q_transport{i});
    CQ_pred_transport(i) = b;
    sigma_CQ_pred_transport(i) = sigma_b;
    tauit_pred_transport(i) = -a/CQ_pred_transport(i);
    sigma_tauit_pred_transport(i) = sigma_a/CQ_pred_transport(i);
    Q_tau_RMSE_transport(i) = sqrt(mean((Q_pred-Q_transport{i}).^2/length(Q_pred))); %RMSE in fit
    R = corrcoef(tau_transport{i},Q_transport{i}); %r^2 in fit
    Q_tau_r2_transport(i) = R(2).^2; %r^2 in fit
    tau_conf_min_transport{i} = [tauit_pred_transport(i)+sigma_tauit_pred_transport(i); tau_transport{i}];
    tau_conf_max_transport{i} = [tauit_pred_transport(i)-sigma_tauit_pred_transport(i); tau_transport{i}];
    Q_conf_min_transport{i} = [0; Q_pred-sigma_Q_pred];
    Q_conf_max_transport{i} = [0; Q_pred+sigma_Q_pred];
    %fit to unbinned values above first full tau bin
    ind_transport_unbinned = find(tauRe_all{i}>tau_bins_min(tau_bin_ind(1)));
    tau_transport_unbinned = tauRe_all{i}(ind_transport_unbinned);
    Q_transport_unbinned = Q_all{i}(ind_transport_unbinned);
    [~, ~, ~, ~, Q_pred, ~] = linearfit(tau_transport_unbinned, Q_transport_unbinned); %fit without including error points
    Q_tau_RMSE_transport_unbinned(i) = sqrt(mean((Q_pred-Q_transport_unbinned).^2/length(Q_pred))); %RMSE in fit - unbinned
    R = corrcoef(tau_transport_unbinned,Q_transport_unbinned); %r^2 in fit
    Q_tau_r2_transport_unbinned(i) = R(2).^2; %r2 in fit - unbinned
    
    %fit for continuous flux bins
    tau_bin_ind = find(fQ_tau_bin_avg{i}>=fQ_continuous); %get tau bins that are continuous
    tau_continuous{i} = tau_bins_mid(tau_bin_ind)';
    sigma_tau_continuous{i} = ones(size(tau_bin_ind))*mode(tau_bins_max-tau_bins_min);
    Q_continuous{i} = Q_tau_bin_avg{i}(tau_bin_ind);
    sigma_Q_continuous{i} = Q_tau_bin_SE{i}(tau_bin_ind);
    [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(tau_continuous{i}, Q_continuous{i}, sigma_tau_continuous{i}, sigma_Q_continuous{i});
    CQ_pred_continuous(i) = b;
    sigma_CQ_pred_continuous(i) = sigma_b;
    tauit_pred_continuous(i) = -a/CQ_pred_continuous(i);
    sigma_tauit_pred_continuous(i) = sigma_a/CQ_pred_continuous(i);
    tau_conf_min_continuous{i} = [tauit_pred_continuous(i)+sigma_tauit_pred_continuous(i); tau_continuous{i}];
    tau_conf_max_continuous{i} = [tauit_pred_continuous(i)-sigma_tauit_pred_continuous(i); tau_continuous{i}];
    Q_conf_min_continuous{i} = [0; Q_pred-sigma_Q_pred];
    Q_conf_max_continuous{i} = [0; Q_pred+sigma_Q_pred];
    
    %fit with known threshold
    tauex_bin_ind = find(fQ_tauex_bin_avg{i}>=fQ_continuous); %get tau bins that are continuous
    tauex_continuous{i} = tauex_bins_mid(tauex_bin_ind)';
    Q_tauex_continuous{i} = Q_tauex_bin_avg{i}(tauex_bin_ind);
    CQ_pred_fixedtauit(i) = mean(Q_tauex_continuous{i}./tauex_continuous{i});
    sigma_CQ_pred_fixedtauit(i) = std(Q_tauex_continuous{i}./tauex_continuous{i})/length(tauex_continuous{i});
%     tauex_bin_ind = find(fQ_tauex_bin_avg{i}>=fQ_transport); %get tau bins that have transport
%     tauex_transport{i} = tauex_bins_mid(tauex_bin_ind)';
%     Q_tauex_transport{i} = Q_tauex_bin_avg{i}(tauex_bin_ind);
%     CQ_pred_fixedtauit(i) = mean(Q_tauex_transport{i}./tauex_transport{i}); %use transport bins
%     sigma_CQ_pred_fixedtauit(i) = std(Q_tauex_transport{i}./tauex_transport{i})/length(tauex_transport{i});
    
    %fit values for etatau_inactive
    ind_good = find(~isnan(Q_etatau_inactive_bin_avg{i}));
    [~, b, ~, sigma_b, Q_pred, ~] = ...
        linearfit(etatau_bins_mid(ind_good)', Q_etatau_inactive_bin_avg{i}(ind_good),...
        etatau_bins_max(ind_good)'-etatau_bins_min(ind_good)', Q_etatau_inactive_bin_SE{i}(ind_good));
    CQ_pred_etatau_inactive(i) = b;
    sigma_CQ_pred_etatau_inactive(i) = sigma_b;
    Q_etatau_RMSE_inactive(i) = sqrt(mean((Q_pred-Q_etatau_inactive_bin_avg{i}(ind_good)).^2/length(ind_good))); %RMSE in fit
    R = corrcoef(etatau_bins_mid(ind_good),Q_etatau_inactive_bin_avg{i}(ind_good)); %r^2 in fit
    Q_etatau_r2_inactive(i) =  R(2).^2; %r^2 in fit
    %fit to unbinned values with etatau greater than 0
    ind_inactive_unbinned = ind_transport_unbinned; %use same values as for Q-tau comparison
    %ind_inactive_unbinned = find(etatau_inactive_all{i}>0); %reselect values based on eta*tau>0
    etatau_inactive_unbinned = etatau_inactive_all{i}(ind_inactive_unbinned);
    Q_inactive_unbinned = Q_all{i}(ind_inactive_unbinned);
    [~, ~, ~, ~, Q_pred, ~] = linearfit(etatau_inactive_unbinned, Q_inactive_unbinned); %fit without including error points
    Q_etatau_RMSE_inactive_unbinned(i) = sqrt(mean((Q_pred-Q_inactive_unbinned).^2/length(Q_pred))); %RMSE in fit - unbinned
    R = corrcoef(etatau_inactive_unbinned,Q_inactive_unbinned); %r^2 in fit
    Q_etatau_r2_inactive_unbinned(i) = R(2).^2; %r2 in fit - unbinned
    
    %fit values for etatau_active
    ind_good = find(~isnan(Q_etatau_active_bin_avg{i}));
    [~, b, ~, sigma_b, Q_pred, ~] = ...
        linearfit(etatau_bins_mid(ind_good)', Q_etatau_active_bin_avg{i}(ind_good),...
        etatau_bins_max(ind_good)'-etatau_bins_min(ind_good)', Q_etatau_active_bin_SE{i}(ind_good));
    CQ_pred_etatau_active(i) = b;
    sigma_CQ_pred_etatau_active(i) = sigma_b;
    Q_etatau_RMSE_active(i) = sqrt(mean((Q_pred-Q_etatau_active_bin_avg{i}(ind_good)).^2/length(ind_good))); %RMSE in fit
    R = corrcoef(etatau_bins_mid(ind_good),Q_etatau_active_bin_avg{i}(ind_good)); %r^2 in fit
    Q_etatau_r2_active(i) =  R(2).^2; %r^2 in fit
    %fit to unbinned values with etatau greater than 0
    ind_active_unbinned = ind_transport_unbinned; %use same values as for Q-tau comparison
    %ind_active_unbinned = find(etatau_active_all{i}>0); %reselect values based on eta*tau>0
    etatau_active_unbinned = etatau_active_all{i}(ind_active_unbinned);
    Q_active_unbinned = Q_all{i}(ind_active_unbinned);
    [~, ~, ~, ~, Q_pred, ~] = linearfit(etatau_active_unbinned, Q_active_unbinned); %fit without including error points
    Q_etatau_RMSE_active_unbinned(i) = sqrt(mean((Q_pred-Q_active_unbinned).^2/length(Q_pred))); %RMSE in fit - unbinned
    R = corrcoef(etatau_active_unbinned,Q_active_unbinned); %r^2 in fit
    Q_etatau_r2_active_unbinned(i) = R(2).^2; %r2 in fit - unbinned
    
    %fit to TFEM thresholds
    fQ_TFEM_ind{i} = intersect(intersect(find(fQ_altbins_mid>fQ_TFEM_fit_min(i)),find(fQ_altbins_mid<fQ_TFEM_fit_max(i))),find(~isnan(tauth_fQ_altbin_avg{i})));
    if i==3
        [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_altbins_mid(fQ_TFEM_ind{i})', tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}), fQ_altbins_max(fQ_TFEM_ind{i})'-fQ_altbins_min(fQ_TFEM_ind{i})', tauth_fQ_altbin_std{i}(fQ_TFEM_ind{i}));
    else %don't use confidence intervals for Jeri and RG, because there is insufficient data for these
        [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_altbins_mid(fQ_TFEM_ind{i})', tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}));
    end
    tauft_TFEM(i) = a;
    tauit_TFEM(i) = a+b;
    tauth_fit_TFEM{i} = tauth_fit;
    sigma_tauth_fit_TFEM{i} = sigma_tauth_fit;
    
    %fit values for zq versus ust
    ind_fit = find(fQ_ust_bin_avg{i}>0.5);
    [a, b, ~, sigma_b, ~, ~] = linearfit(ust_bins_mid(ind_fit)', zq_ust_bin_avg{i}(ind_fit), ust_bins_max(ind_fit)'-ust_bins_min(ind_fit)', zq_ust_bin_SE{i}(ind_fit));
    linearfit_zq_ust{i} = [a b]; %linear fit for zq versus ust
    intercept_zq_ust{i} = a; %intercept of fit
    slope_zq_ust{i} = b; %slope of fit
    sigmaslope_zq_ust{i} = sigma_b; %uncertainty in slope of fit
    
    %fit values for zq versus tau_ex
    ind_fit = find(fQ_tauex_bin_avg{i}>0.5);
    [a, b, ~, ~, ~, ~] = linearfit(tauex_bins_mid(ind_fit)', zq_tauex_bin_avg{i}(ind_fit));
    linearfit_zq_tauex{i} = [a b]; %linear fit for zq versus tau_ex
end
% 
% %% FOR EACH SITE, PERFORM SUB-BINNING WITHIN FREQUENCY AND ETATAURATIO BINS, CREATE SITE-SPECIFIC PLOTS 
% for i = 1:N_Sites
%     
%     %% FQ SUBBINS
%     %get indices of fQ bins that are not empty
%     ind_fQ_bins_full = find(1-cellfun(@isempty,Q_fQ_bin_values{i}))';
%     
%     %fluxes into tau bins within fQ bins
%     Q_fQ_tau_bin_avg = cell(N_fQ_bins,1);
%     Q_fQ_tau_bin_std = cell(N_fQ_bins,1);
%     Q_fQ_tau_bin_SE = cell(N_fQ_bins,1);
%     zq_fQ_tau_bin_avg = cell(N_fQ_bins,1);
%     zq_fQ_tau_bin_std = cell(N_fQ_bins,1);
%     zq_fQ_tau_bin_SE = cell(N_fQ_bins,1);
%     uratio_fQ_tau_bin_avg = cell(N_fQ_bins,1);
%     uratio_fQ_tau_bin_std = cell(N_fQ_bins,1);
%     uratio_fQ_tau_bin_SE = cell(N_fQ_bins,1);
%     
%     %fluxes and other values into tauex bins within fQ bins
%     Q_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
%     Q_fQ_tauex_bin_std = cell(N_fQ_bins,1);
%     Q_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
%     eta_inactive_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
%     eta_inactive_fQ_tauex_bin_std = cell(N_fQ_bins,1);
%     eta_inactive_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
%     eta_active_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
%     eta_active_fQ_tauex_bin_std = cell(N_fQ_bins,1);
%     eta_active_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
%     zs_fQ_tauex_bin_avg = cell(N_fQ_bins,1);
%     zs_fQ_tauex_bin_std = cell(N_fQ_bins,1);
%     zs_fQ_tauex_bin_SE = cell(N_fQ_bins,1);
%     
%     %fluxes into eta*tau bins within fQ bins
%     Q_fQ_etatau_inactive_bin_avg = cell(N_fQ_bins,1);
%     Q_fQ_etatau_inactive_bin_std = cell(N_fQ_bins,1);
%     Q_fQ_etatau_inactive_bin_SE = cell(N_fQ_bins,1);
%     Q_fQ_etatau_active_bin_avg = cell(N_fQ_bins,1);
%     Q_fQ_etatau_active_bin_std = cell(N_fQ_bins,1);
%     Q_fQ_etatau_active_bin_SE = cell(N_fQ_bins,1);
%     
%     %go through each fQ bin
%     for j = 1:N_fQ_bins
%         
%         %initialize tau bins for each fQ bin
%         Q_fQ_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         Q_fQ_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         Q_fQ_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         zq_fQ_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         zq_fQ_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         zq_fQ_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_fQ_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_fQ_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_fQ_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
% 
%         %initialize tau_ex bins for each fQ bin
%         Q_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
%         Q_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
%         Q_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
%         eta_inactive_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
%         eta_inactive_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
%         eta_inactive_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
%         eta_active_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
%         eta_active_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
%         eta_active_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
%         zs_fQ_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
%         zs_fQ_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
%         zs_fQ_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
%                 
%         %initialize etatau bins for each fQ bin
%         Q_fQ_etatau_inactive_bin_avg{j} = zeros(N_etatau_bins,1)*NaN;
%         Q_fQ_etatau_inactive_bin_std{j} = zeros(N_etatau_bins,1)*NaN;
%         Q_fQ_etatau_inactive_bin_SE{j} = zeros(N_etatau_bins,1)*NaN;
%         Q_fQ_etatau_active_bin_avg{j} = zeros(N_etatau_bins,1)*NaN;
%         Q_fQ_etatau_active_bin_std{j} = zeros(N_etatau_bins,1)*NaN;
%         Q_fQ_etatau_active_bin_SE{j} = zeros(N_etatau_bins,1)*NaN;
%         
%         %go through each tau bin
%         for k=1:N_tau_bins
%             bin_ind = find(tau_fQ_bin_values{i}{j}>=tau_bins_min(k)&tau_fQ_bin_values{i}{j}<=tau_bins_max(k));
%             if ~isempty(bin_ind)
%                 Q_fQ_tau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
%                 Q_fQ_tau_bin_values = Q_fQ_tau_bin_values(~isnan(Q_fQ_tau_bin_values));
%                 Q_fQ_tau_bin_avg{j}(k) = mean(Q_fQ_tau_bin_values);
%                 Q_fQ_tau_bin_SE{j}(k) = std(Q_fQ_tau_bin_values)/sqrt(length(Q_fQ_tau_bin_values));
%                 Q_fQ_tau_bin_std{j}(k) = std(Q_fQ_tau_bin_values);
%                 
%                 zq_fQ_tau_bin_values = zq_fQ_bin_values{i}{j}(bin_ind);
%                 zq_fQ_tau_bin_values = zq_fQ_tau_bin_values(~isnan(zq_fQ_tau_bin_values));
%                 zq_fQ_tau_bin_avg{j}(k) = mean(zq_fQ_tau_bin_values);
%                 zq_fQ_tau_bin_SE{j}(k) = std(zq_fQ_tau_bin_values)/sqrt(length(zq_fQ_tau_bin_values));
%                 zq_fQ_tau_bin_std{j}(k) = std(zq_fQ_tau_bin_values);
%                 
%                 uratio_fQ_tau_bin_values = uratio_fQ_bin_values{i}{j}(bin_ind);
%                 uratio_fQ_tau_bin_values = uratio_fQ_tau_bin_values(~isnan(uratio_fQ_tau_bin_values));
%                 uratio_fQ_tau_bin_avg{j}(k) = mean(uratio_fQ_tau_bin_values);
%                 uratio_fQ_tau_bin_SE{j}(k) = std(uratio_fQ_tau_bin_values)/sqrt(length(uratio_fQ_tau_bin_values));
%                 uratio_fQ_tau_bin_std{j}(k) = std(uratio_fQ_tau_bin_values);
%             end
%         end
%         
%         %go through each tau_ex bin
%         for k=1:N_tauex_bins
%             bin_ind = find(tauex_fQ_bin_values{i}{j}>=tauex_bins_min(k)&tauex_fQ_bin_values{i}{j}<=tauex_bins_max(k));
%             if ~isempty(bin_ind)
%                 Q_fQ_tauex_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
%                 Q_fQ_tauex_bin_values = Q_fQ_tauex_bin_values(~isnan(Q_fQ_tauex_bin_values));
%                 Q_fQ_tauex_bin_avg{j}(k) = mean(Q_fQ_tauex_bin_values);
%                 Q_fQ_tauex_bin_SE{j}(k) = std(Q_fQ_tauex_bin_values)/sqrt(length(Q_fQ_tauex_bin_values));
%                 Q_fQ_tauex_bin_std{j}(k) = std(Q_fQ_tauex_bin_values);
%                 
%                 eta_inactive_fQ_tauex_bin_values = eta_inactive_fQ_bin_values{i}{j}(bin_ind);
%                 eta_inactive_fQ_tauex_bin_values = eta_inactive_fQ_tauex_bin_values(~isnan(eta_inactive_fQ_tauex_bin_values));
%                 eta_inactive_fQ_tauex_bin_avg{j}(k) = mean(eta_inactive_fQ_tauex_bin_values);
%                 eta_inactive_fQ_tauex_bin_SE{j}(k) = std(eta_inactive_fQ_tauex_bin_values)/sqrt(length(eta_inactive_fQ_tauex_bin_values));
%                 eta_inactive_fQ_tauex_bin_std{j}(k) = std(eta_inactive_fQ_tauex_bin_values);
%                 
%                 eta_active_fQ_tauex_bin_values = eta_active_fQ_bin_values{i}{j}(bin_ind);
%                 eta_active_fQ_tauex_bin_values = eta_active_fQ_tauex_bin_values(~isnan(eta_active_fQ_tauex_bin_values));
%                 eta_active_fQ_tauex_bin_avg{j}(k) = mean(eta_active_fQ_tauex_bin_values);
%                 eta_active_fQ_tauex_bin_SE{j}(k) = std(eta_active_fQ_tauex_bin_values)/sqrt(length(eta_active_fQ_tauex_bin_values));
%                 eta_active_fQ_tauex_bin_std{j}(k) = std(eta_active_fQ_tauex_bin_values);
%                 
%                 zs_fQ_tauex_bin_values = zs_fQ_bin_values{i}{j}(bin_ind);
%                 zs_fQ_tauex_bin_values = zs_fQ_tauex_bin_values(~isnan(zs_fQ_tauex_bin_values));
%                 zs_fQ_tauex_bin_avg{j}(k) = exp(mean(log(zs_fQ_tauex_bin_values))); %get mean in log space
%                 zs_fQ_tauex_bin_SE{j}(k) = std(zs_fQ_tauex_bin_values)/sqrt(length(zs_fQ_tauex_bin_values));
%                 zs_fQ_tauex_bin_std{j}(k) = std(zs_fQ_tauex_bin_values);
%             end
%         end
%         
%         %go through each etatau bin
%         for k=1:N_etatau_bins
%             bin_ind = find(eta_inactive_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}>=etatau_bins_min(k)&...
%                 eta_inactive_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}<=etatau_bins_max(k));
%             if ~isempty(bin_ind)
%                 Q_fQ_eta_inactive_tau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
%                 Q_fQ_eta_inactive_tau_bin_values = Q_fQ_eta_inactive_tau_bin_values(~isnan(Q_fQ_eta_inactive_tau_bin_values));
%                 Q_fQ_etatau_inactive_bin_avg{j}(k) = mean(Q_fQ_eta_inactive_tau_bin_values);
%                 Q_fQ_etatau_inactive_bin_SE{j}(k) = std(Q_fQ_eta_inactive_tau_bin_values)/sqrt(length(Q_fQ_eta_inactive_tau_bin_values));
%                 Q_fQ_etatau_inactive_bin_std{j}(k) = std(Q_fQ_eta_inactive_tau_bin_values);
%             end
%             
%             bin_ind = find(eta_active_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}>=etatau_bins_min(k)&...
%                 eta_active_fQ_bin_values{i}{j}.*tau_fQ_bin_values{i}{j}<=etatau_bins_max(k));
%             if ~isempty(bin_ind)
%                 Q_fQ_eta_active_tau_bin_values = Q_fQ_bin_values{i}{j}(bin_ind);
%                 Q_fQ_eta_active_tau_bin_values = Q_fQ_eta_active_tau_bin_values(~isnan(Q_fQ_eta_active_tau_bin_values));
%                 Q_fQ_etatau_active_bin_avg{j}(k) = mean(Q_fQ_eta_active_tau_bin_values);
%                 Q_fQ_etatau_active_bin_SE{j}(k) = std(Q_fQ_eta_active_tau_bin_values)/sqrt(length(Q_fQ_eta_active_tau_bin_values));
%                 Q_fQ_etatau_active_bin_std{j}(k) = std(Q_fQ_eta_active_tau_bin_values);
%             end
%         end    
%     end
% 
%     %% ETATAURATIO SUBBINS
%     %get indices of etatauratio bins that are not empty
%     ind_etatauratio_bins_full = find(cellfun(@length,Q_etatauratio_bin_values{i})>1)';
%     
%     %fluxes into tau bins within etatauratio bins
%     Q_etatauratio_tau_bin_avg = cell(N_etatauratio_bins,1);
%     Q_etatauratio_tau_bin_std = cell(N_etatauratio_bins,1);
%     Q_etatauratio_tau_bin_SE = cell(N_etatauratio_bins,1);
%     zq_etatauratio_tau_bin_avg = cell(N_etatauratio_bins,1);
%     zq_etatauratio_tau_bin_std = cell(N_etatauratio_bins,1);
%     zq_etatauratio_tau_bin_SE = cell(N_etatauratio_bins,1);
%     ubar_etatauratio_tau_bin_avg = cell(N_etatauratio_bins,1);
%     ubar_etatauratio_tau_bin_std = cell(N_etatauratio_bins,1);
%     ubar_etatauratio_tau_bin_SE = cell(N_etatauratio_bins,1);
%     ustd_etatauratio_tau_bin_avg = cell(N_etatauratio_bins,1);
%     ustd_etatauratio_tau_bin_std = cell(N_etatauratio_bins,1);
%     ustd_etatauratio_tau_bin_SE = cell(N_etatauratio_bins,1);
%     uratio_etatauratio_tau_bin_avg = cell(N_etatauratio_bins,1);
%     uratio_etatauratio_tau_bin_std = cell(N_etatauratio_bins,1);
%     uratio_etatauratio_tau_bin_SE = cell(N_etatauratio_bins,1);
%     
%     %fluxes and other values into tauex bins within etatauratio bins
%     Q_etatauratio_tauex_bin_avg = cell(N_etatauratio_bins,1);
%     Q_etatauratio_tauex_bin_std = cell(N_etatauratio_bins,1);
%     Q_etatauratio_tauex_bin_SE = cell(N_etatauratio_bins,1);
%     zs_etatauratio_tauex_bin_avg = cell(N_etatauratio_bins,1);
%     zs_etatauratio_tauex_bin_std = cell(N_etatauratio_bins,1);
%     zs_etatauratio_tauex_bin_SE = cell(N_etatauratio_bins,1);
%     ubar_etatauratio_tauex_bin_avg = cell(N_etatauratio_bins,1);
%     ubar_etatauratio_tauex_bin_std = cell(N_etatauratio_bins,1);
%     ubar_etatauratio_tauex_bin_SE = cell(N_etatauratio_bins,1);
%     ustd_etatauratio_tauex_bin_avg = cell(N_etatauratio_bins,1);
%     ustd_etatauratio_tauex_bin_std = cell(N_etatauratio_bins,1);
%     ustd_etatauratio_tauex_bin_SE = cell(N_etatauratio_bins,1);
%     uratio_etatauratio_tauex_bin_avg = cell(N_etatauratio_bins,1);
%     uratio_etatauratio_tauex_bin_std = cell(N_etatauratio_bins,1);
%     uratio_etatauratio_tauex_bin_SE = cell(N_etatauratio_bins,1);
%     
%     %go through each etatauratio bin
%     for j = 1:N_etatauratio_bins
%         
%         %initialize tau bins for each etatauratio bin
%         Q_etatauratio_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         Q_etatauratio_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         Q_etatauratio_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         zq_etatauratio_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         zq_etatauratio_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         zq_etatauratio_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         ubar_etatauratio_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         ubar_etatauratio_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         ubar_etatauratio_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         ustd_etatauratio_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         ustd_etatauratio_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         ustd_etatauratio_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_etatauratio_tau_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_etatauratio_tau_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_etatauratio_tau_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
% 
%         %initialize tau_ex bins for each etatauratio bin
%         Q_etatauratio_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
%         Q_etatauratio_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
%         Q_etatauratio_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
%         zs_etatauratio_tauex_bin_avg{j} = zeros(N_tauex_bins,1)*NaN;
%         zs_etatauratio_tauex_bin_std{j} = zeros(N_tauex_bins,1)*NaN;
%         zs_etatauratio_tauex_bin_SE{j} = zeros(N_tauex_bins,1)*NaN;
%         ubar_etatauratio_tauex_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         ubar_etatauratio_tauex_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         ubar_etatauratio_tauex_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         ustd_etatauratio_tauex_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         ustd_etatauratio_tauex_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         ustd_etatauratio_tauex_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_etatauratio_tauex_bin_avg{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_etatauratio_tauex_bin_std{j} = zeros(N_tau_bins,1)*NaN;
%         uratio_etatauratio_tauex_bin_SE{j} = zeros(N_tau_bins,1)*NaN;
%          
%         %go through each tau bin
%         for k=1:N_tau_bins
%             bin_ind = find(tau_etatauratio_bin_values{i}{j}>=tau_bins_min(k)&tau_etatauratio_bin_values{i}{j}<=tau_bins_max(k));
%             if ~isempty(bin_ind)
%                 Q_etatauratio_tau_bin_values = Q_etatauratio_bin_values{i}{j}(bin_ind);
%                 Q_etatauratio_tau_bin_values = Q_etatauratio_tau_bin_values(~isnan(Q_etatauratio_tau_bin_values));
%                 Q_etatauratio_tau_bin_avg{j}(k) = mean(Q_etatauratio_tau_bin_values);
%                 Q_etatauratio_tau_bin_SE{j}(k) = std(Q_etatauratio_tau_bin_values)/sqrt(length(Q_etatauratio_tau_bin_values));
%                 Q_etatauratio_tau_bin_std{j}(k) = std(Q_etatauratio_tau_bin_values);
%                 
%                 zq_etatauratio_tau_bin_values = zq_etatauratio_bin_values{i}{j}(bin_ind);
%                 zq_etatauratio_tau_bin_values = zq_etatauratio_tau_bin_values(~isnan(zq_etatauratio_tau_bin_values));
%                 zq_etatauratio_tau_bin_avg{j}(k) = mean(zq_etatauratio_tau_bin_values);
%                 zq_etatauratio_tau_bin_SE{j}(k) = std(zq_etatauratio_tau_bin_values)/sqrt(length(zq_etatauratio_tau_bin_values));
%                 zq_etatauratio_tau_bin_std{j}(k) = std(zq_etatauratio_tau_bin_values);
%                 
%                 ubar_etatauratio_tau_bin_values = ubar_etatauratio_bin_values{i}{j}(bin_ind);
%                 ubar_etatauratio_tau_bin_values = ubar_etatauratio_tau_bin_values(~isnan(ubar_etatauratio_tau_bin_values));
%                 ubar_etatauratio_tau_bin_avg{j}(k) = mean(ubar_etatauratio_tau_bin_values);
%                 ubar_etatauratio_tau_bin_SE{j}(k) = std(ubar_etatauratio_tau_bin_values)/sqrt(length(ubar_etatauratio_tau_bin_values));
%                 ubar_etatauratio_tau_bin_std{j}(k) = std(ubar_etatauratio_tau_bin_values);
%                 
%                 ustd_etatauratio_tau_bin_values = ustd_etatauratio_bin_values{i}{j}(bin_ind);
%                 ustd_etatauratio_tau_bin_values = ustd_etatauratio_tau_bin_values(~isnan(ustd_etatauratio_tau_bin_values));
%                 ustd_etatauratio_tau_bin_avg{j}(k) = mean(ustd_etatauratio_tau_bin_values);
%                 ustd_etatauratio_tau_bin_SE{j}(k) = std(ustd_etatauratio_tau_bin_values)/sqrt(length(ustd_etatauratio_tau_bin_values));
%                 ustd_etatauratio_tau_bin_std{j}(k) = std(ustd_etatauratio_tau_bin_values);
%                 
%                 uratio_etatauratio_tau_bin_values = uratio_etatauratio_bin_values{i}{j}(bin_ind);
%                 uratio_etatauratio_tau_bin_values = uratio_etatauratio_tau_bin_values(~isnan(uratio_etatauratio_tau_bin_values));
%                 uratio_etatauratio_tau_bin_avg{j}(k) = mean(uratio_etatauratio_tau_bin_values);
%                 uratio_etatauratio_tau_bin_SE{j}(k) = std(uratio_etatauratio_tau_bin_values)/sqrt(length(uratio_etatauratio_tau_bin_values));
%                 uratio_etatauratio_tau_bin_std{j}(k) = std(uratio_etatauratio_tau_bin_values);
%             end
%         end
%         
%         %go through each tau_ex bin
%         for k=1:N_tauex_bins
%             bin_ind = find(tauex_etatauratio_bin_values{i}{j}>=tauex_bins_min(k)&tauex_etatauratio_bin_values{i}{j}<=tauex_bins_max(k));
%             if ~isempty(bin_ind)
%                 Q_etatauratio_tauex_bin_values = Q_etatauratio_bin_values{i}{j}(bin_ind);
%                 Q_etatauratio_tauex_bin_values = Q_etatauratio_tauex_bin_values(~isnan(Q_etatauratio_tauex_bin_values));
%                 Q_etatauratio_tauex_bin_avg{j}(k) = mean(Q_etatauratio_tauex_bin_values);
%                 Q_etatauratio_tauex_bin_SE{j}(k) = std(Q_etatauratio_tauex_bin_values)/sqrt(length(Q_etatauratio_tauex_bin_values));
%                 Q_etatauratio_tauex_bin_std{j}(k) = std(Q_etatauratio_tauex_bin_values);
%                 
%                 zs_etatauratio_tauex_bin_values = zs_etatauratio_bin_values{i}{j}(bin_ind);
%                 zs_etatauratio_tauex_bin_values = zs_etatauratio_tauex_bin_values(~isnan(zs_etatauratio_tauex_bin_values));
%                 zs_etatauratio_tauex_bin_avg{j}(k) = exp(mean(log(zs_etatauratio_tauex_bin_values))); %get mean in log space
%                 zs_etatauratio_tauex_bin_SE{j}(k) = std(zs_etatauratio_tauex_bin_values)/sqrt(length(zs_etatauratio_tauex_bin_values));
%                 zs_etatauratio_tauex_bin_std{j}(k) = std(zs_etatauratio_tauex_bin_values);
%                 
%                 ubar_etatauratio_tauex_bin_values = ubar_etatauratio_bin_values{i}{j}(bin_ind);
%                 ubar_etatauratio_tauex_bin_values = ubar_etatauratio_tauex_bin_values(~isnan(ubar_etatauratio_tauex_bin_values));
%                 ubar_etatauratio_tauex_bin_avg{j}(k) = mean(ubar_etatauratio_tauex_bin_values);
%                 ubar_etatauratio_tauex_bin_SE{j}(k) = std(ubar_etatauratio_tauex_bin_values)/sqrt(length(ubar_etatauratio_tauex_bin_values));
%                 ubar_etatauratio_tauex_bin_std{j}(k) = std(ubar_etatauratio_tauex_bin_values);
%                 
%                 ustd_etatauratio_tauex_bin_values = ustd_etatauratio_bin_values{i}{j}(bin_ind);
%                 ustd_etatauratio_tauex_bin_values = ustd_etatauratio_tauex_bin_values(~isnan(ustd_etatauratio_tauex_bin_values));
%                 ustd_etatauratio_tauex_bin_avg{j}(k) = mean(ustd_etatauratio_tauex_bin_values);
%                 ustd_etatauratio_tauex_bin_SE{j}(k) = std(ustd_etatauratio_tauex_bin_values)/sqrt(length(ustd_etatauratio_tauex_bin_values));
%                 ustd_etatauratio_tauex_bin_std{j}(k) = std(ustd_etatauratio_tauex_bin_values);
%                 
%                 uratio_etatauratio_tauex_bin_values = uratio_etatauratio_bin_values{i}{j}(bin_ind);
%                 uratio_etatauratio_tauex_bin_values = uratio_etatauratio_tauex_bin_values(~isnan(uratio_etatauratio_tauex_bin_values));
%                 uratio_etatauratio_tauex_bin_avg{j}(k) = mean(uratio_etatauratio_tauex_bin_values);
%                 uratio_etatauratio_tauex_bin_SE{j}(k) = std(uratio_etatauratio_tauex_bin_values)/sqrt(length(uratio_etatauratio_tauex_bin_values));
%                 uratio_etatauratio_tauex_bin_std{j}(k) = std(uratio_etatauratio_tauex_bin_values);
%             end
%         end    
%     end
%     
%     %% PLOTS %%
%     
%     %transport flux plot
%     figure; clf; hold on;
%     errorbar(tau_transport{i},Q_transport{i},sigma_Q_transport{i},Markers{i}); %plot values
%     plot([tauit(i) max(tau_transport{i})],[0 CQ_pred_fixedtauit(i)*(max(tau_transport{i})-tauit(i))],LineColors{i}); %plot fit for fixed tauit
%     plot([sigma_tauit_pred_transport(i) max(tau_transport{i})],[0 CQ_pred_transport(i)*(max(tau_transport{i})-sigma_tauit_pred_transport(i))],'k'); %plot fit with error bars
%     plot(tau_conf_min_transport{i},Q_conf_min_transport{i},'k--',tau_conf_max_transport{i},Q_conf_max_transport{i},'k--'); %include confidence intervals
%     xlabel('\tau (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     title(Sites{i});
%     legend_items = {'transport data','\tau_{it,TFEM} fit','transport fit','conf intervals',};
%     legend(legend_items,'Location','NorthWest');
%     set(gca,'FontSize',16);
%     print([folder_Plots,'flux_tau_transport',Sites{i},'.png'],'-dpng');
%     
%     %continuous flux plot
%     figure; clf; hold on;
%     errorbar(tau_continuous{i},Q_continuous{i},sigma_Q_continuous{i},Markers{i}); %plot values
%     plot([tauit(i) max(tau_continuous{i})],[0 CQ_pred_fixedtauit(i)*(max(tau_continuous{i})-tauit(i))],LineColors{i}); %plot fit for fixed tauit
%     plot([sigma_tauit_pred_continuous(i) max(tau_continuous{i})],[0 CQ_pred_continuous(i)*(max(tau_continuous{i})-sigma_tauit_pred_continuous(i))],'k'); %plot fit with error bars
%     plot(tau_conf_min_continuous{i},Q_conf_min_continuous{i},'k--',tau_conf_max_continuous{i},Q_conf_max_continuous{i},'k--'); %include confidence intervals
%     xlabel('\tau (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     title(Sites{i});
%     legend_items = {'continuous data','\tau_{it,TFEM} fit','continuous fit','conf intervals',};
%     legend(legend_items,'Location','NorthWest');
%     set(gca,'FontSize',16);
%     print([folder_Plots,'flux_tau_continuous',Sites{i},'.png'],'-dpng');
%     
%     %all flux plot
%     figure; clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tau_bins_mid,Q_fQ_tau_bin_avg{j},Q_fQ_tau_bin_SE{j},Markers_fQ_bins{j});
%     end
%     plot([tauit(i) max(tau_continuous{i})],[0 CQ_pred_fixedtauit(i)*(max(tau_continuous{i})-tauit(i))],'k'); %plot fit for fixed tauit
%     xlabel('\tau (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     title(Sites{i});
%     legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
%     legend_items{length(legend_items)+1} = 'Q(\tau-\tau_{it}) fit';
%     legend(legend_items,'Location','NorthWest');
%     set(gca,'FontSize',16);
%     print([folder_Plots,'Flux_Tau_fQ_',Sites{i},'.png'],'-dpng'); 
% 
%     %flux height versus tau
%     figure; clf; hold on;
%     zq_max = max(max([zq_fQ_tau_bin_avg{:}]+[zq_fQ_tau_bin_std{:}])); %get maximum zq
%     zq_max = ceil(zq_max*100+1)/100; %round up, add 1 cm
%     for j = ind_fQ_bins_full
%        errorbar(tau_bins_mid,zq_fQ_tau_bin_avg{j},zq_fQ_tau_bin_SE{j},Markers_fQ_bins{j});
%     end
%     ylim([0 zq_max]);
%     xlabel('\tau (Pa)');
%     ylabel('z_{q} (m)');
%     title(Sites{i});
%     legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
%     h_legend = legend(legend_items,'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'FluxHeight_Tau_fQ_',Sites{i},'.png'],'-dpng'); 
%     
%     %flux versus tauex
%     figure; clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tauex_bins_mid,Q_fQ_tauex_bin_avg{j},Q_fQ_tauex_bin_SE{j},Markers_fQ_bins{j});
%     end
%     ylim([0 max(Q_all{i})]);
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthWest');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'Flux_Tauex_fQ_',Sites{i},'.png'],'-dpng');
%   
%     %eta versus tauex
%     figure; clf;
%     tauex_min = min([tauex_bins_min(min([find(~isnan(Q_fQ_tauex_bin_avg{ind_fQ_bins_full(1)}));...
%         find(~isnan(Q_fQ_tauex_bin_avg{ind_fQ_bins_full(1)}))]));0]);
%     tauex_max = tauex_bins_max(max([find(~isnan(Q_fQ_tauex_bin_avg{end}));...
%         find(~isnan(Q_fQ_tauex_bin_avg{end}))]));
%     eta_max = max([eta_inactive_fQ_tauex_bin_avg{end}+eta_inactive_fQ_tauex_bin_std{end};...
%         eta_active_fQ_tauex_bin_avg{end}+eta_active_fQ_tauex_bin_std{end}]); %get maximum eta
%     %method 1 - inactive
%     subplot(1,2,1); hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tauex_bins_mid,eta_inactive_fQ_tauex_bin_avg{j},eta_inactive_fQ_tauex_bin_SE{j},Markers_fQ_bins{j});
%     end
%     plot(tauex_bins_mid,tauratio_tauex_bin_avg{i},'k'); %plot tauex/tau versus tauex
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('\eta_{inactive}');
%     xlim([tauex_min tauex_max]);
%     ylim([0 eta_max]);
%     title(Sites{i});
%     legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
%     legend_items{length(legend_items)+1} = '\tau_{ex}/\tau';
%     legend(legend_items,'Location','NorthWest');
%     set(gca,'FontSize',16);
%     %method 2 - active
%     subplot(1,2,2); hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tauex_bins_mid,eta_active_fQ_tauex_bin_avg{j},eta_active_fQ_tauex_bin_SE{j},Markers_fQ_bins{j});
%     end
%     plot(tauex_bins_mid,tauratio_tauex_bin_avg{i},'k'); %plot tauex/tau versus tauex
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('\eta_{active}');
%     xlim([tauex_min tauex_max]);
%     ylim([0 eta_max]);
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     set(gcf,'PaperPosition',[0 0 16 6]);
%     print([folder_Plots,'eta_Tauex_fQ_',Sites{i},'.png'],'-dpng');
% 
%     %zs versus tau_ex
%     figure; clf; hold on;
%     for j = ind_fQ_bins_full
%        errorbar(tauex_bins_mid,zs_fQ_tauex_bin_avg{j},zs_fQ_tauex_bin_SE{j},Markers_fQ_bins{j});
%     end
%     set(gca,'yscale','log');
%     if i==1
%        ylim([1e-6 1e-2]);
%     elseif i==3
%        ylim([1e-5 1e-2]);
%     end
%     xlim([tauex_min tauex_max]);
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('z_{s} (m)');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'zs_Tauex_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %flux versus eta*tau
%     figure; clf; %initialize plot
%     etatau_max = etatau_bins_max(max([find(~isnan(Q_fQ_etatau_inactive_bin_avg{end}));...
%         find(~isnan(Q_fQ_etatau_active_bin_avg{end}))]));
%     %inactive method
%     subplot(1,2,1); hold on;
%     for j = ind_fQ_bins_full
%        errorbar(etatau_bins_mid,Q_fQ_etatau_inactive_bin_avg{j},Q_fQ_etatau_inactive_bin_SE{j},Markers_fQ_bins{j});
%     end
%     plot([0 etatau_max],[0 etatau_max]*CQ_pred_fixedtauit(i),'k'); %fit prediction with CQ
%     xlim([0 etatau_max]);
%     ylim([0, max(Q_all{i})]);
%     xlabel('\eta_{inactive}\tau (Pa)');
%     ylabel('Q (g/m/s)');
%     title(Sites{i});
%     legend_items = fQ_bins_legend{i}(ind_fQ_bins_full);
%     legend_items{length(legend_items)+1}='C_{Q}\eta\tau';
%     legend(legend_items,'Location','SouthEast');
%     set(gca,'FontSize',16);
%     title([Sites{i}]);
%     %active method
%     subplot(1,2,2); hold on;
%     for j = ind_fQ_bins_full
%        errorbar(etatau_bins_mid,Q_fQ_etatau_active_bin_avg{j},Q_fQ_etatau_active_bin_SE{j},Markers_fQ_bins{j});
%     end
%     plot([0 etatau_max],[0 etatau_max]*CQ_pred_fixedtauit(i),'k'); %fit prediction with CQ
%     xlim([0 etatau_max]);
%     ylim([0, max(Q_all{i})]);
%     xlabel('\eta_{active}\tau (Pa)');
%     ylabel('Q (g m^{-1} s^{-1}');
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     title([Sites{i}]);
%     set(gcf,'PaperPosition',[0 0 14 6]);
%     print([folder_Plots,'flux_etatau_',Sites{i},'.png'],'-dpng');
% 
%     %flux versus tau comparison to flux versus eta-tau
%     figure; clf; hold on;
%     errorbar(tauex_bins_mid,Q_tauex_bin_avg{i},Q_tauex_bin_SE{i},Markers{1});
%     errorbar(etatau_bins_mid,Q_etatau_active_bin_avg{i},Q_etatau_active_bin_SE{i},Markers{2});
%     errorbar(etatau_bins_mid,Q_etatau_inactive_bin_avg{i},Q_etatau_inactive_bin_SE{i},Markers{3});
%     xlabel('\tau_{eff} (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     legend(['\tau_{ex}, r^2=',num2str(round(Q_tau_r2_transport(i),2)),', RMSE = ',num2str(round(Q_tau_RMSE_transport(i),2))],...
%         ['\eta_{active}\tau, r^2=',num2str(round(Q_etatau_r2_active(i),2)),', RMSE = ',num2str(round(Q_etatau_RMSE_active(i),2))],...
%         ['\eta_{inactive}\tau, r^2=',num2str(round(Q_etatau_r2_inactive(i),2)),', RMSE = ',num2str(round(Q_etatau_RMSE_inactive(i),2))],...
%         'Location','NorthWest');
%     taueff_max = max([tauex_bins_max(max(find(~isnan(Q_tauex_bin_avg{i})))); etatau_bins_max(max(find(~isnan(Q_etatau_active_bin_avg{i})))); etatau_bins_max(max(find(~isnan(Q_etatau_inactive_bin_avg{i}))))]);
%     xlim([0 taueff_max]);
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     print([folder_Plots,'flux_comparison_',Sites{i},'.png'],'-dpng');
%     
%     %plot as panels
%     Q_max = max([Q_tauex_bin_avg{i}+Q_tauex_bin_SE{i};Q_etatau_active_bin_avg{i}+Q_etatau_active_bin_SE{i};Q_etatau_inactive_bin_avg{i}+Q_etatau_inactive_bin_SE{i}]);
%     figure; clf;
%     subplot(1,3,1);
%     errorbar(tauex_bins_mid,Q_tauex_bin_avg{i},Q_tauex_bin_SE{i},Markers{1});
%     xlim([0 taueff_max]);
%     ylim([0 Q_max]);
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     set(gca,'FontSize',16);
%     subplot(1,3,2);
%     errorbar(etatau_bins_mid,Q_etatau_active_bin_avg{i},Q_etatau_active_bin_SE{i},Markers{2});
%     xlim([0 taueff_max]);
%     ylim([0 Q_max]);
%     xlabel('\eta_{active}\tau (Pa)');
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     subplot(1,3,3);
%     errorbar(etatau_bins_mid,Q_etatau_inactive_bin_avg{i},Q_etatau_inactive_bin_SE{i},Markers{3});
%     xlim([0 taueff_max]);
%     ylim([0 Q_max]);
%     xlabel('\eta_{inactive}\tau (Pa)');
%     set(gca,'FontSize',16);
%     set(gcf,'PaperPosition',[0 0 14 6]);
%     print([folder_Plots,'flux_comparison_panels_',Sites{i},'.png'],'-dpng');
%     
%     %flux versus tau comparison to flux versus eta-tau - unbinned data
%     figure; clf; hold on;
%     plot(tauex_all{i},Q_all{i},Markers{1});
%     plot(etatau_active_all{i},Q_all{i},Markers{2});
%     plot(etatau_inactive_all{i},Q_all{i},Markers{3});
%     xlabel('\tau_{eff} (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     legend(['\tau_{ex}, r^2=',num2str(round(Q_tau_r2_transport_unbinned(i),2)),', RMSE = ',num2str(round(Q_tau_RMSE_transport_unbinned(i),2))],...
%         ['\eta_{active}\tau, r^2=',num2str(round(Q_etatau_r2_active_unbinned(i),2)),', RMSE = ',num2str(round(Q_etatau_RMSE_active_unbinned(i),2))],...
%         ['\eta_{inactive}\tau, r^2=',num2str(round(Q_etatau_r2_inactive_unbinned(i),2)),', RMSE = ',num2str(round(Q_etatau_RMSE_inactive_unbinned(i),2))],...
%         'Location','SouthEast');
%     taueff_max = max([tauex_all{i}; etatau_active_all{i}; etatau_active_all{i}]);
%     xlim([0 taueff_max]);
%     ylim([0 max(Q_all{i})]);
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     title([Sites{i}]);
%     print([folder_Plots,'flux_comparison_unbinned_',Sites{i},'.png'],'-dpng');
%     
%     %plot as panels
%     Q_max = max(Q_all{i});
%     figure; clf;
%     subplot(1,3,1);
%     plot(tauex_all{i},Q_all{i},Markers{1});
%     xlim([0 taueff_max]);
%     ylim([0 Q_max]);
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('Q (g m^{-1} s^{-1})');
%     set(gca,'FontSize',16);
%     subplot(1,3,2);
%     plot(etatau_active_all{i},Q_all{i},Markers{2});
%     xlim([0 taueff_max]);
%     ylim([0 Q_max]);
%     xlabel('\eta_{active}\tau (Pa)');
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     subplot(1,3,3);
%     plot(etatau_inactive_all{i},Q_all{i},Markers{3});
%     xlim([0 taueff_max]);
%     ylim([0 Q_max]);
%     xlabel('\eta_{inactive}\tau (Pa)');
%     set(gca,'FontSize',16);
%     xlim([0 taueff_max]);
%     ylim([0 max(Q_all{i})]);
%     set(gca,'FontSize',16);
%     set(gcf,'PaperPosition',[0 0 14 6]);
%     print([folder_Plots,'flux_comparison_unbinned_panels_',Sites{i},'.png'],'-dpng');
%     
%     %ubar/ustd versus tau
%     figure; clf; hold on; %initialize plot
%     for j = ind_fQ_bins_full
%        errorbar(tau_bins_mid,uratio_fQ_tau_bin_avg{j},uratio_fQ_tau_bin_SE{j},Markers_fQ_bins{j});
%     end
%     xlabel('\tau (Pa)');
%     ylabel('\mu_{u}/\sigma_{u}');
%     title(Sites{i});
%     h_legend = legend(fQ_bins_legend{i}(ind_fQ_bins_full),'Location','NorthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'uratio_tau_fQ_',Sites{i},'.png'],'-dpng');
%     
%     %ubar versus tau
%     figure; clf; hold on; %initialize plot
%     for j = ind_etatauratio_bins_full
%         errorbar(tau_bins_mid,ubar_etatauratio_tau_bin_avg{j},ubar_etatauratio_tau_bin_SE{j},Markers_etatauratio_bins{j});
%     end
%     xlabel('\tau (Pa)');
%     ylabel('\mu_{u} (m/s)');
%     title(Sites{i});
%     h_legend = legend(etatauratio_bins_legend{i}(ind_etatauratio_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'ubar_tau_etatauratio_',Sites{i},'.png'],'-dpng');
%     
%     %ustd versus tau
%     figure; clf; hold on; %initialize plot
%     for j = ind_etatauratio_bins_full
%         errorbar(tau_bins_mid,ustd_etatauratio_tau_bin_avg{j},ustd_etatauratio_tau_bin_SE{j},Markers_etatauratio_bins{j});
%     end
%     xlabel('\tau (Pa)');
%     ylabel('\sigma_{u} (m/s)');
%     title(Sites{i});
%     h_legend = legend(etatauratio_bins_legend{i}(ind_etatauratio_bins_full),'Location','SouthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'ustd_tau_etatauratio_',Sites{i},'.png'],'-dpng');
%     
%     %uratio versus tau
%     figure; clf; hold on; %initialize plot
%     for j = ind_etatauratio_bins_full
%         errorbar(tau_bins_mid,uratio_etatauratio_tau_bin_avg{j},uratio_etatauratio_tau_bin_SE{j},Markers_etatauratio_bins{j});
%     end
%     xlabel('\tau (Pa)');
%     ylabel('\mu_{u}/\sigma_{u}');
%     title(Sites{i});
%     h_legend = legend(etatauratio_bins_legend{i}(ind_etatauratio_bins_full),'Location','NorthEast');
%     set(gca,'FontSize',16);
%     set(h_legend,'FontSize',16);
%     print([folder_Plots,'uratio_tau_etatauratio_',Sites{i},'.png'],'-dpng');
%     
%     %tauth versus fQ
%     figure; clf; hold on; %initialize plot
%     errorbar(fQ_altbins_mid,tauth_fQ_altbin_avg{i},tauth_fQ_altbin_SE{i},Markers{i});
%     plot([0 1],[tauft_TFEM(i) tauit_TFEM(i)],'k','LineWidth',2);
%     plot(fQ_altbins_mid(fQ_TFEM_ind{i}),tauth_fit_TFEM{i}+sigma_tauth_fit_TFEM{i},'k--',...
%         fQ_altbins_mid(fQ_TFEM_ind{i}),tauth_fit_TFEM{i}-sigma_tauth_fit_TFEM{i},'k--');
%     xlabel('f_{Q}');
%     ylabel('\tau_{th,TFEM}');
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     print([folder_Plots,'tauth_fQ_',Sites{i},'.png'],'-dpng');
%         
%     %frequency of u ranges versus fQ
%     figure; clf; hold on; %initialize plot
%     errorbar(fQ_altbins_mid,f_uaboveft_avg{i},f_uaboveft_SE{i},'og');
%     errorbar(fQ_altbins_mid,f_ubelowit_avg{i},f_ubelowit_SE{i},'xr');
%     errorbar(fQ_altbins_mid,f_ubetween_frombelow_avg{i},f_ubetween_frombelow_SE{i},'^b');
%     errorbar(fQ_altbins_mid,f_ubetween_fromabove_avg{i},f_ubetween_fromabove_SE{i},'vb');
%     ylim([0 1]);
%     xlabel('f_{Q}');
%     ylabel('f_{u range}');
%     legend('above ft','below it','from below','from above');
%     title(Sites{i});
%     set(gca,'FontSize',16);
%     print([folder_Plots,'tauth_f_urange_',Sites{i},'.png'],'-dpng');
% end

%% COMBINED PLOTS
%re-allow figures
set(0,'DefaultFigureVisible', 'on');

%raw flux versus stress
figure; clf; hold on;
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    errorbar(tau_transport{i},Q_transport{i},sigma_Q_transport{i},Markers{i}); %plot values
    plot([tauit_pred_transport(i) max(tau_transport{i})],[0 CQ_pred_transport(i)*(max(tau_transport{i})-tauit_pred_transport(i))],LineColors{i}); %plot fit for these values
    legend_items{2*i-1} = Sites{i}; %add to list of legend items
    legend_items{2*i} = 'fit'; %add to list of legend items
end
xlabel('\tau (Pa)');
ylabel('Q (g/m/s)');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',20);
set(gcf, 'PaperPosition',[0 0 8 5]);
print([folder_Plots,'flux_tau_transport_all.png'],'-dpng');

%flux versus stress for continuous tau bins
figure; clf; hold on;
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    errorbar(tau_continuous{i},Q_continuous{i},sigma_Q_continuous{i},Markers{i}); %plot values
    plot([tauit(i) max(tau_continuous{i})],[0 CQ_pred_fixedtauit(i)*(max(tau_continuous{i})-tauit(i))],LineColors{i}); %plot fit for fixed tauit
    legend_items{2*i-1} = Sites{i}; %add to list of legend items
    legend_items{2*i} = 'fit'; %add to list of legend items
end
plot(tau_Greeley96,Q_fit_Greeley96,Marker_Greeley96);
legend_items{length(legend_items)+1} = 'Greeley (96)';
plot(tau_Namikas03,Q_fit_Namikas03,Marker_Namikas03);
legend_items{length(legend_items)+1} = 'Namikas (03)';
plot(tau_Farrell12,Q_fit_Farrell12,Marker_Farrell12);
legend_items{length(legend_items)+1} = 'Farrell (12)';
xlabel('\tau (Pa)');
ylabel('Q (g/m/s)');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'flux_tau_continuous_all.png'],'-dpng');

%zq versus ust
figure; clf; hold on;
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>0.5);
    errorbar(ust_bins_mid(ind_plot),zq_ust_bin_avg{i}(ind_plot),zq_ust_bin_SE{i}(ind_plot),Markers{i},'MarkerSize',10);
end
%plot Greeley, Namikas, and Farrell
plot(ust_Greeley96,zbar_Greeley96,Marker_Greeley96,'MarkerSize',10);
plot(ust_Namikas03,zbar_Namikas03,Marker_Namikas03,'MarkerSize',10);
plot(ust_Farrell12,zbar_Farrell12,Marker_Farrell12,'MarkerSize',10);
%plot fits to field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>0.5);
    ust_fit = [ust_bins_mid(min(ind_plot)) ust_bins_mid(max(ind_plot))];
    plot(ust_fit,linearfit_zq_ust{i}(1)+linearfit_zq_ust{i}(2).*ust_fit,LineColors{i}); %plot fit for fixed tauit
end
%fit to Greeley, Namikas, and Farrell
[intercept slope] = linearfit(ust_Greeley96,zbar_Greeley96);
ust_fit_Greeley96 = [min(ust_Greeley96) max(ust_Greeley96)];
zbar_fit_Greeley96 = intercept+slope*ust_fit_Greeley96;
[intercept slope] = linearfit(ust_Namikas03,zbar_Namikas03);
ust_fit_Namikas03 = [min(ust_Namikas03) max(ust_Namikas03)];
zbar_fit_Namikas03 = intercept+slope*ust_fit_Namikas03;
[intercept slope] = linearfit(ust_Farrell12,zbar_Farrell12);
ust_fit_Farrell12 = [min(ust_Farrell12) max(ust_Farrell12)];
zbar_fit_Farrell12 = intercept+slope*ust_fit_Farrell12;
%plot fits
plot(ust_fit_Greeley96, zbar_fit_Greeley96,'m');
plot(ust_fit_Namikas03, zbar_fit_Namikas03,'k');
plot(ust_fit_Farrell12, zbar_fit_Farrell12,'c');
%organize plot
xlim([0.25 0.6]);
ylim([0 0.14]);
xlabel('shear velocity, u_{*} (m/s)');
ylabel('saltation height, z_q (m)');
legend_items = Sites;
legend_items{length(legend_items)+1} = 'Greeley et al. (1996)';
legend_items{length(legend_items)+1} = 'Namikas (2003)';
legend_items{length(legend_items)+1} = 'Farrell et al. (2012)';
legend(legend_items,'Location','SouthEast');
set(gca,'FontSize',20);
set(gcf,'PaperPosition',[0 0 11 7]);
print([folder_Plots,'zq_ust_all.png'],'-dpng');

%zqnorm versus ust
figure; clf; hold on;
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>0.5);
    errorbar(ust_bins_mid(ind_plot),1e3*zq_ust_bin_avg{i}(ind_plot)/d50_site(i),1e3*zq_ust_bin_SE{i}(ind_plot)/d50_site(i),Markers{i},'MarkerSize',10); %normalized by single grain size for site
    %errorbar(ust_bins_mid(ind_plot),zqnorm_ust_bin_avg{i}(ind_plot),zqnorm_ust_bin_SE{i}(ind_plot),Markers{i},'MarkerSize',10); %individually normalized values
end
%plot Greeley, Namikas, and Farrell
plot(ust_Greeley96,1e3*zbar_Greeley96/d50_Greeley96,Marker_Greeley96,'MarkerSize',10);
plot(ust_Namikas03,1e3*zbar_Namikas03/d50_Namikas03,Marker_Namikas03,'MarkerSize',10);
plot(ust_Farrell12,1e3*zbar_Farrell12/d50_Farrell12,Marker_Farrell12,'MarkerSize',10);
%plot fits to field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>0.5);
    ust_fit = [ust_bins_mid(min(ind_plot)) ust_bins_mid(max(ind_plot))];
    zqnorm_fit = 1e3*(intercept_zq_ust{i}+slope_zq_ust{i}.*ust_fit)/d50_site(i);
    plot(ust_fit,zqnorm_fit,LineColors{i}); %plot fit for fixed tauit
end
%plot fits to Greeley, Namikas, and Farrell
plot(ust_fit_Greeley96, 1e3*zbar_fit_Greeley96/d50_Greeley96,'m');
plot(ust_fit_Namikas03, 1e3*zbar_fit_Namikas03/d50_Namikas03,'k');
plot(ust_fit_Farrell12, 1e3*zbar_fit_Farrell12/d50_Farrell12,'c');
%organize plot
xlim([0.25 0.6]);
ylim([0 250]);
xlabel('shear velocity, u_{*} (m/s)');
ylabel('dimensionless salt. ht, z_q/d_{50}');
legend_items = Sites;
legend_items{length(legend_items)+1} = 'Greeley et al. (1996)';
legend_items{length(legend_items)+1} = 'Namikas (2003)';
legend_items{length(legend_items)+1} = 'Farrell et al. (2012)';
legend(legend_items,'Location','SouthEast');
set(gca,'FontSize',20);
set(gcf,'PaperPosition',[0 0 10 5.5]);
print([folder_Plots,'zqnorm_ust_all.png'],'-dpng');

% %zq versus tauex
% figure; clf; hold on;
% legend_items = cell(N_Sites*2,1);
% for i = 1:N_Sites
%     tauex_fit = [tauex_bins_mid(min(find(fQ_tauex_bin_avg{i}>0.5))) tauex_bins_mid(max(find(fQ_tauex_bin_avg{i}>0.5)))];
%     r_fit = corrcoef(tauex_bins_mid(find(fQ_tauex_bin_avg{i}>0.5)), zq_tauex_bin_avg{i}(find(fQ_tauex_bin_avg{i}>0.5)));
%     r2_fit = (r_fit(2)).^2
%     zq_bin_bar = mean(zq_tauex_bin_avg{i}(find(fQ_tauex_bin_avg{i}>0.5)))
%     errorbar(tauex_bins_mid,zq_tauex_bin_avg{i},zq_tauex_bin_SE{i},Markers{i});
%     plot(tauex_fit,linearfit_zq_tauex{i}(1)+linearfit_zq_tauex{i}(2).*tauex_fit,LineColors{i}); %plot fit for fixed tauit
%     legend_items{2*i-1} = Sites{i}; %add to list of legend items
%     legend_items{2*i} = 'fit'; %add to list of legend items
% end
% plot(tauex_Greeley96,zbar_Greeley96,Marker_Greeley96);
% legend_items{length(legend_items)+1} = 'Greeley et al (1996)';
% plot(tauex_Namikas03,zbar_Namikas03,Marker_Namikas03);
% legend_items{length(legend_items)+1} = 'Namikas (2003)';
% plot(tauex_Farrell12,zbar_Farrell12,Marker_Farrell12);
% legend_items{length(legend_items)+1} = 'Farrell et al (2012)';
% xlim([0 max(tauex_bins_max)]);
% ylim([0 0.14]);
% xlabel('\tau_{ex} (Pa)');
% ylabel('z_q (m)');
% legend(legend_items,'Location','eastoutside');
% set(gca,'FontSize',20);
% set(gcf,'PaperPosition',[0 0 10 5]);
% print([folder_Plots,'zq_tauex_all.png'],'-dpng');
% 
% %zq/d50 versus tauex
% figure; clf; hold on;
% legend_items = cell(N_Sites,1);
% for i = 1:N_Sites
% %    errorbar(tauex_bins_mid,zqnorm_tauex_bin_avg{i},zqnorm_tauex_bin_SE{i},Markers{i});
%     errorbar(tauex_bins_mid,1e3*zq_tauex_bin_avg{i}/d50_site(i),1e3*zq_tauex_bin_SE{i}/d50_site(i),Markers{i});
% %    plot([0 max(tauex_bins_max)],zqnorm_bar(i)*[1 1],LineColors{i}); %plot fit for fixed tauit
%     legend_items{i} = Sites{i}; %add to list of legend items
% end
% plot(tauex_Greeley96,1e3*zbar_Greeley96/d50_Greeley96,Marker_Greeley96);
% legend_items{length(legend_items)+1} = 'Greeley et al (1996)';
% plot(tauex_Namikas03,1e3*zbar_Namikas03/d50_Namikas03,Marker_Namikas03);
% legend_items{length(legend_items)+1} = 'Namikas (2003)';
% plot(tauex_Farrell12,1e3*zbar_Farrell12/d50_Farrell12,Marker_Farrell12);
% legend_items{length(legend_items)+1} = 'Farrell et al (2012)';
% xlabel('\tau_{ex} (Pa)');
% ylabel('z_q/d_{50}');
% xlim([0 max(tauex_bins_max)]);
% legend(legend_items,'Location','SouthEast');
% set(gca,'FontSize',20);
% set(gcf, 'PaperPosition',[0 0 10 6]);
% print([folder_Plots,'zqnorm_tauex_all.png'],'-dpng');

%flux versus tauex
figure; clf; hold on;
for i = 1:N_Sites
    errorbar(tauex_bins_mid,Q_tauex_bin_avg{i},Q_tauex_bin_SE{i},Markers{i},'MarkerSize',10);
end
for i = 1:N_Sites
    tauex_fit = [0 0.3];
    plot(tauex_fit, CQ_pred_fixedtauit(i)*tauex_fit,LineColors{i}); %plot linear fit value
    %plot(tauex_fit, CQ_pred_transport(i)*tauex_fit,LineColors{i}); %plot linear fit value
end
xlim([0 0.3]);
ylim([0 60]);
%plot(tauex_Greeley96,Q_fit_Greeley96,Marker_Greeley96);
%legend_items{length(legend_items)+1} = 'Greeley (96)';
%plot(tauex_Namikas03,Q_fit_Namikas03,Marker_Namikas03);
%legend_items{length(legend_items)+1} = 'Namikas (03)';
%plot(tauex_Farrell12,Q_fit_Farrell12,Marker_Farrell12);
%legend_items{length(legend_items)+1} = 'Farrell (12)';
xlabel('excess shear stress, \tau_{ex} (Pa)');
ylabel('saltation mass flux, Q (g m^{-1} s^{-1})');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',20);
set(gcf, 'PaperPosition',[0 0 11 7]);
print([folder_Plots,'Q_tauex_all.png'],'-dpng');

%normalized flux versus tauex
figure; clf; hold on;
for i = 1:N_Sites
    ind_plot = find(fQ_tauex_bin_avg{i}>=fQ_continuous);
    errorbar(tauex_bins_mid(ind_plot),Qnorm_tauex_bin_avg{i}(ind_plot),Qnorm_tauex_bin_SE{i}(ind_plot),Markers{i},'MarkerSize',10);
    %Qnorm_bar = mean(Qnorm_tauex_bin_avg{i}(intersect(find(tauex_bins_mid>0.02),find(~isnan(Qnorm_tauex_bin_avg{i})))))
end
tauex_fit = [0, max(tauex_bins_max)];
for i = 1:N_Sites
    %tauex_fit = [min(tauex_bins_mid(ind_plot)), max(tauex_bins_mid(ind_plot))];
    C = 1e-3*CQ_pred_fixedtauit(i)*g/sqrt(tauit(i)/rho_a)
    C_SE = 1e-3*sigma_CQ_pred_fixedtauit(i)*g/sqrt(tauit(i)/rho_a)
    plot(tauex_fit,C*ones(2,1),LineColors{i});
end
xlabel('excess shear stress, $$\tau_{ex}$$ (Pa)','Interpreter','Latex');
ylabel('dimensionless saltation flux, $$\hat{Q}$$','Interpreter','Latex');
legend(Sites,'Location','SouthEast');
xlim([0 max(tauex_bins_max)]);
ylim([0 10]);
set(gca,'FontSize',20);
set(gcf, 'PaperPosition',[0 0 7 7]);
print([folder_Plots,'Qnorm_tauex_all.png'],'-dpng');

%flux/tauit versus tauex
figure; clf; hold on;
legend_items = cell(N_Sites,1);
for i = 1:N_Sites
    errorbar(tauex_bins_mid,1e-3*Q_tauex_bin_avg{i}/tauit_pred_transport(i),1e-3*Q_tauex_bin_SE{i}/tauit_pred_transport(i),Markers{i});
    legend_items{i} = Sites{i};
end
plot(tauex_Greeley96,1e-3*Q_fit_Greeley96./tauit_Greeley96,Marker_Greeley96);
legend_items{length(legend_items)+1} = 'Greeley (96)';
plot(tauex_Namikas03,1e-3*Q_fit_Namikas03./tauit_Namikas03,Marker_Namikas03);
legend_items{length(legend_items)+1} = 'Namikas (03)';
plot(tauex_Farrell12,1e-3*Q_fit_Farrell12./tauit_Farrell12,Marker_Farrell12);
legend_items{length(legend_items)+1} = 'Farrell (12)';
xlabel('\tau_{ex} (Pa)');
ylabel('Q/\tau_{it} (s^{-1})');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'Qtauit_tauex_all.png'],'-dpng');

%flux/ustit versus tauex
figure; clf; hold on;
legend_items = cell(N_Sites,1);
for i = 1:N_Sites
    errorbar(tauex_bins_mid,Q_tauex_bin_avg{i}/(sqrt(tauit_pred_transport(i)/rho_a)),...
        Q_tauex_bin_SE{i}/(sqrt(tauit_pred_transport(i)/rho_a)),Markers{i});
    legend_items{i} = Sites{i};
end
plot(tauex_Greeley96,Q_fit_Greeley96./(sqrt(tauit_Greeley96/rho_a)),Marker_Greeley96);
legend_items{length(legend_items)+1} = 'Greeley (96)';
plot(tauex_Namikas03,Q_fit_Namikas03./(sqrt(tauit_Namikas03/rho_a)),Marker_Namikas03);
legend_items{length(legend_items)+1} = 'Namikas (03)';
xlabel('\tau_{ex} (Pa)');
ylabel('Q/u_{*,it} (s^{-1})');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'Qustit_tauex_all.png'],'-dpng');

%flux/d50 versus tauex
figure; clf; hold on;
legend_items = cell(N_Sites,1);
for i = 1:N_Sites
    errorbar(tauex_bins_mid,Q_tauex_bin_avg{i}/(mean(d50_all{i}(~isnan(d50_all{i})))/1e3),...
        Q_tauex_bin_SE{i}/(mean(d50_all{i}(~isnan(d50_all{i})))/1e3),Markers{i});
    legend_items{i} = Sites{i};
end
plot(tauex_Greeley96,Q_fit_Greeley96/(d50_Greeley96/1e3),Marker_Greeley96);
legend_items{length(legend_items)+1} = 'Greeley (96)';
plot(tauex_Namikas03,Q_fit_Namikas03/(d50_Namikas03/1e3),Marker_Namikas03);
legend_items{length(legend_items)+1} = 'Namikas (03)';
xlabel('\tau_{ex} (Pa)');
ylabel('Q/d_{50} (g m^{-2} s^{-1})');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',16);
set(gcf, 'PaperPosition',[0 0 10 6]);
print([folder_Plots,'Qd50_tauex_all.png'],'-dpng');

%flux/zq versus tauex
figure; clf; hold on;
legend_items = cell(N_Sites,1);
for i = 1:N_Sites
    errorbar(tauex_bins_mid,Q_tauex_bin_avg{i}/(mean(zq_all{i}(zq_all{i}>0))),...
        Q_tauex_bin_SE{i}/(mean(zq_all{i}(zq_all{i}>0))),Markers{i});
    legend_items{i} = Sites{i};
end
plot(tauex_Greeley96,Q_fit_Greeley96/mean(zbar_Greeley96),Marker_Greeley96);
legend_items{length(legend_items)+1} = 'Greeley (96)';
plot(tauex_Namikas03,Q_fit_Namikas03/mean(zbar_Namikas03),Marker_Namikas03);
legend_items{length(legend_items)+1} = 'Namikas (03)';
xlabel('\tau_{ex} (Pa)');
ylabel('Q/z_{Q} (g m^{-2} s^{-1})');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'Qzq_tauex_all.png'],'-dpng');

%eta versus tauex
figure; clf; hold on;
for i = 1:N_Sites
    errorbar(tauex_bins_mid,eta_tauex_bin_avg{i},eta_tauex_bin_SE{i},Markers{i});
end
xlabel('\tau_{ex} (Pa)');
ylabel('\eta');
legend(Sites,'Location','SouthEast');
set(gca,'FontSize',16);
print([folder_Plots,'eta_tauex_all.png'],'-dpng');

%flux versus etatau
figure; clf; hold on;
for i = 1:N_Sites
    errorbar(etatau_bins_mid,Q_etatau_active_bin_avg{i},Q_etatau_active_bin_SE{i},Markers{i});
    %plot([0 max(etatau_bins_max)],[0 CQ_pred_fixedtauit(i)*max(etatau_bins_max)],LineColors{i}); %plot fit for fixed tauit
    legend_items{2*i-1} = Sites{i}; %add to list of legend items
    legend_items{2*i} = 'fit'; %add to list of legend items
end
xlabel('\eta\tau (Pa)');
ylabel('Q (g/m/s)');
%legend(legend_items,'Location','NorthWest');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'Q_etatau_all.png'],'-dpng');

%flux frequency versus eta
figure; clf; hold on;
for i=1:N_Sites
    errorbar(eta_bins_mid,fQ_eta_inactive_bin_avg{i},fQ_eta_inactive_bin_SE{i},Markers{i});
end
xlabel('\eta');
ylabel('f_{Q}');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',16);
print([folder_Plots,'fQ_eta_all.png'],'-dpng');

%tauth versus fQ
figure; clf; hold on; %initialize plot
%legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    errorbar(fQ_altbins_mid(fQ_TFEM_ind{i}),tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}),tauth_fQ_altbin_SE{i}(fQ_TFEM_ind{i}),Markers{i},'MarkerSize',10);
    %legend_items{2*i-1} = Sites{i}; %add to list of legend items
    %legend_items{2*i} = 'fit'; %add to list of legend items
end
for i = 1:N_Sites
    plot([0 1],[tauft_TFEM(i) tauit_TFEM(i)],LineColors{i});
end
xlabel('frequency of transport');
ylabel('TFEM inferred threshold stress, \tau_{th}');
legend(Sites,'Location','NorthEast');
set(gca,'FontSize',20);
set(gcf, 'PaperPosition',[0 0 10 6]);
print([folder_Plots,'tauth_fQ_all.png'],'-dpng');

%re-allow figures
set(0,'DefaultFigureVisible', 'on');