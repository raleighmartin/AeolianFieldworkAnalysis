%% ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;
close all;

%parameters
rho_a = 1.23; %air density (kg/m^3)
g = 9.8; %gravity (m/s^2)

%frequency parameters
fQ_transport = 0.1; %minimum mean fQ for tauex bin to be called "transport"
fQ_moderate = 0.5; %mininum mean fQ for tauex bin to be called "moderate"
fQ_continuous = 0.9; %minimum mean fQ for tauex bin to be called "continuous"

%information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for storing data output
folder_Plots = '../../PlotOutput/StressFlux/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows'));
N_Sites = length(Sites);

%load external data
load(strcat(folder_AnalysisData,'LitData')); %literature data
LitNames = {'Greeley et al. (1996)'; 'Namikas (2003)';'Farrell et al. (2012)'};
ust_lit = {ust_Greeley96, ust_Namikas03, ust_Farrell12};
zq_lit = {zbar_Greeley96, zbar_Namikas03, zbar_Farrell12};
zqnorm_lit = {1e3*zbar_Greeley96/d50_Greeley96, 1e3*zbar_Namikas03/d50_Namikas03, 1e3*zbar_Farrell12/d50_Farrell12};
N_Lit = 3;

%set info for plotting
Markers_field = {'s','d','o'};
Colors_field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
MarkerSize_field = 6;
LineWidth_field = 0.75;
Markers_lit = {'x','+','*'};
MarkerSize_lit = 6;
Colors_lit = {[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
LineWidth_lit = 0.75;
PlotFont = 12;

%old info for plotting
Markers = {'bx','ro','gv'};
LineColors = {'b','r','g'};

%% PARAMETERS
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


%% COMPUTE ZQNORM BASED ON GRAIN SIZE BY SITE
zqnorm_all = cell(N_Sites,1);
for i=1:N_Sites
    zqnorm_all{i} = 1000*zq_all{i}./d50_site(i);
    %zqnorm_all{i} = 1000*zq_all{i}./d50_all{i};
end


%% CREATE BINS
%set fQ bins
fQ_bins_min = 0:0.05:0.95;
fQ_bins_max = 0.05:0.05:1;
fQ_bins_mid = mean([fQ_bins_min; fQ_bins_max]);
N_fQ_bins = length(fQ_bins_mid);
fQ_TFEM_fit_min = [0.65,0.7,0.1]; %mininum for fitting to obtain impact/fluid thresholds
fQ_TFEM_fit_max = [0.9,0.85,0.9]; %maximum for fitting to obtain impact/fluid thresholds


%% FOR EACH SITE, PERFORM BINNING BY USTAR
%create u* bins
ust_bins_min = 0:.01:0.56;
ust_bins_max = 0.01:.01:0.57;
ust_bins_mid = mean([ust_bins_min;ust_bins_max]);
N_ust_bins = length(ust_bins_mid);

%separate fQs into ust bins
fQ_ust_bin_values = cell(N_Sites,1); %fQs into ust bins
N_fQ_ust_bin_values = cell(N_Sites,1); %N fQ's in ust bins
fQ_ust_bin_avg = cell(N_Sites,1); %get average fQ for ust bins
fQ_ust_bin_std = cell(N_Sites,1); %get std fQ for ust bins
fQ_ust_bin_SE = cell(N_Sites,1); %get SE fQ for ust bins

%separate zq's into ust bins
zq_ust_bin_values = cell(N_Sites,1); %zq into ust bins
N_zq_ust_bin_values = cell(N_Sites,1); %N zq's into ust bins
zq_ust_bin_avg = cell(N_Sites,1); %get average zq for ust bins
zq_ust_bin_std = cell(N_Sites,1); %get std zq for ust bins
zq_ust_bin_SE = cell(N_Sites,1); %get SE zq for ust bins

%separate zqnorm's (zq/d50) into ust bins
zqnorm_ust_bin_values = cell(N_Sites,1); %zq/d50s into ust bins
N_zqnorm_ust_bin_values = cell(N_Sites,1); %N zq/d50s in ust bins
zqnorm_ust_bin_avg = cell(N_Sites,1); %get average zq/d50 for ust bins
zqnorm_ust_bin_std = cell(N_Sites,1); %get std zq/d50 for ust bins
zqnorm_ust_bin_SE = cell(N_Sites,1); %get SE zq/d50 for ust bins

%separate Qs into ust bins
Q_ust_bin_values = cell(N_Sites,1); %Qs into ust bins
N_Q_ust_bin_values = cell(N_Sites,1); %N Qs in ust bins
Q_ust_bin_avg = cell(N_Sites,1); %get average Q for ust bins
Q_ust_bin_std = cell(N_Sites,1); %get std Q for ust bins
Q_ust_bin_SE = cell(N_Sites,1); %get SE Q for ust bins

%go through all sites
for i = 1:N_Sites
  
    %initialize fQ's into ust bins
    fQ_ust_bin_values{i} = cell(N_ust_bins,1);
    N_fQ_ust_bin_values{i} = zeros(N_ust_bins,1);
    fQ_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    fQ_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    fQ_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %initialize zq's into ust bins
    zq_ust_bin_values{i} = cell(N_ust_bins,1);
    N_zq_ust_bin_values{i} = zeros(N_ust_bins,1);
    zq_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    zq_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    zq_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %initialize zqnorm's into ust bins
    zqnorm_ust_bin_values{i} = cell(N_ust_bins,1);
    N_zqnorm_ust_bin_values{i} = zeros(N_ust_bins,1);
    zqnorm_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    zqnorm_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    zqnorm_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %initialize Q's into ust bins
    Q_ust_bin_values{i} = cell(N_ust_bins,1);
    N_Q_ust_bin_values{i} = zeros(N_ust_bins,1);
    Q_ust_bin_avg{i} = zeros(N_ust_bins,1)*NaN;
    Q_ust_bin_std{i} = zeros(N_ust_bins,1)*NaN;
    Q_ust_bin_SE{i} = zeros(N_ust_bins,1)*NaN;
    
    %get values, avgs, and standard deviations for ust bins
    for j = 1:N_ust_bins
        %get indices
        bin_ind = find(ustRe_all{i}>=ust_bins_min(j)&ustRe_all{i}<=ust_bins_max(j));

        %get fQs for ust bins
        fQ_ust_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_ust_bin_values{i}{j} = fQ_ust_bin_values{i}{j}(~isnan(fQ_ust_bin_values{i}{j}));
        N_fQ_ust_bin_values{i}(j) = length(fQ_ust_bin_values{i}{j});
        fQ_ust_bin_avg{i}(j) = mean(fQ_ust_bin_values{i}{j});
        fQ_ust_bin_std{i}(j) = std(fQ_ust_bin_values{i}{j});
        
        %get zqs for ust bins
        zq_ust_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_ust_bin_values{i}{j} = zq_ust_bin_values{i}{j}(~isnan(zq_ust_bin_values{i}{j}));
        N_zq_ust_bin_values{i}(j) = length(zq_ust_bin_values{i}{j});
        zq_ust_bin_avg{i}(j) = mean(zq_ust_bin_values{i}{j});
        zq_ust_bin_std{i}(j) = std(zq_ust_bin_values{i}{j});
        
        %get zqnorms for ust bins
        zqnorm_ust_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_ust_bin_values{i}{j} = zqnorm_ust_bin_values{i}{j}(~isnan(zqnorm_ust_bin_values{i}{j}));
        N_zqnorm_ust_bin_values{i}(j) = length(zqnorm_ust_bin_values{i}{j});
        zqnorm_ust_bin_avg{i}(j) = mean(zqnorm_ust_bin_values{i}{j});
        zqnorm_ust_bin_std{i}(j) = std(zqnorm_ust_bin_values{i}{j});
        
        %get Qs for ust bins
        Q_ust_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_ust_bin_values{i}{j} = Q_ust_bin_values{i}{j}(~isnan(Q_ust_bin_values{i}{j}));
        N_Q_ust_bin_values{i}(j) = length(Q_ust_bin_values{i}{j});
        Q_ust_bin_avg{i}(j) = mean(Q_ust_bin_values{i}{j});
        Q_ust_bin_std{i}(j) = std(Q_ust_bin_values{i}{j});
    end

    %get indices for bins with multiple values
    ind_multibin_fQ_ust = find(N_fQ_ust_bin_values{i}>1);
    ind_multibin_zq_ust = find(N_zq_ust_bin_values{i}>1);
    ind_multibin_zqnorm_ust = find(N_zqnorm_ust_bin_values{i}>1);
    ind_multibin_Q_ust = find(N_Q_ust_bin_values{i}>1);
    
    %get standard errors for ust bins
    for j = 1:N_ust_bins
        %get SEs for fQ
        if N_fQ_ust_bin_values{i}(j)>1
            fQ_ust_bin_SE{i}(j) = fQ_ust_bin_std{i}(j)/sqrt(N_fQ_ust_bin_values{i}(j));
        elseif N_fQ_ust_bin_values{i}(j)==1
            ind_below = ind_multibin_fQ_ust(find(ind_multibin_fQ_ust<j,1,'last'));
            ind_above = ind_multibin_fQ_ust(find(ind_multibin_fQ_ust>j,1));
            fQ_ust_bin_SE{i}(j) = mean(fQ_ust_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for zq
        if N_zq_ust_bin_values{i}(j)>1
            zq_ust_bin_SE{i}(j) = zq_ust_bin_std{i}(j)/sqrt(N_zq_ust_bin_values{i}(j));
        elseif N_zq_ust_bin_values{i}(j)==1
            ind_below = ind_multibin_zq_ust(find(ind_multibin_zq_ust<j,1,'last'));
            ind_above = ind_multibin_zq_ust(find(ind_multibin_zq_ust>j,1));
            zq_ust_bin_SE{i}(j) = mean(zq_ust_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for zqnorm
        if N_zqnorm_ust_bin_values{i}(j)>1
            zqnorm_ust_bin_SE{i}(j) = zqnorm_ust_bin_std{i}(j)/sqrt(N_zqnorm_ust_bin_values{i}(j));
        elseif N_zqnorm_ust_bin_values{i}(j)==1
            ind_below = ind_multibin_zqnorm_ust(find(ind_multibin_zqnorm_ust<j,1,'last'));
            ind_above = ind_multibin_zqnorm_ust(find(ind_multibin_zqnorm_ust>j,1));
            zqnorm_ust_bin_SE{i}(j) = mean(zqnorm_ust_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for Q
        if N_Q_ust_bin_values{i}(j)>1
            Q_ust_bin_SE{i}(j) = Q_ust_bin_std{i}(j)/sqrt(N_Q_ust_bin_values{i}(j));
        elseif N_Q_ust_bin_values{i}(j)==1
            ind_below = ind_multibin_Q_ust(find(ind_multibin_Q_ust<j,1,'last'));
            ind_above = ind_multibin_Q_ust(find(ind_multibin_Q_ust>j,1));
            Q_ust_bin_SE{i}(j) = mean(Q_ust_bin_std{i}([ind_below ind_above]));
        end
    end
end


%% PLOT ZQ VERSUS UST - PANEL A
figure; subplot(10,1,1:3); hold on;

%plot binned field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>=fQ_moderate);
    errorbar(ust_bins_mid(ind_plot),zq_ust_bin_avg{i}(ind_plot),zq_ust_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot literature data
for i = 1:N_Lit
    plot(ust_lit{i},zq_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end

%fit values for zq versus ust
intercept_zq_ust = zeros(N_Sites,1); %intercept of fit
sigmaintercept_zq_ust = zeros(N_Sites,1); %uncertainty in intercept of fit
slope_zq_ust = zeros(N_Sites,1); %slope of fit for zq versus ustar
sigmaslope_zq_ust = zeros(N_Sites,1); %uncertainty in slope for zq versus ustar

for i = 1:N_Sites
    %fit values for zq versus ust
    ind_fit = find(fQ_ust_bin_avg{i}>=fQ_moderate);
    ust_fit = ust_bins_mid(ind_fit)';
    sigma_ust_fit = ust_bins_max(ind_fit)'-ust_bins_min(ind_fit)';
    zq_fit = zq_ust_bin_avg{i}(ind_fit);
    sigma_zq_fit = zq_ust_bin_SE{i}(ind_fit);

    [a, b, sigma_a, sigma_b, zq_pred, sigma_zq_pred] = linearfit(ust_bins_mid(ind_fit)', zq_ust_bin_avg{i}(ind_fit), ust_bins_max(ind_fit)'-ust_bins_min(ind_fit)', zq_ust_bin_SE{i}(ind_fit));
    intercept_zq_ust(i) = a; %intercept of fit
    sigmaintercept_zq_ust(i) = sigma_a; %uncertainty in intercept of fit
    slope_zq_ust(i) = b; %slope of fit
    sigmaslope_zq_ust(i) = sigma_b; %uncertainty in slope of fit
    
    plot(ust_fit,zq_pred,'Color',Colors_field{i},'LineWidth',LineWidth_field);
end

%fit to Greeley, Namikas, and Farrell
for i = 1:N_Lit
    [intercept, slope] = linearfit(ust_lit{i},zq_lit{i});
    ust_fit = [min(ust_lit{i}), max(ust_lit{i})];
    zq_fit = intercept+slope*ust_fit;
    plot(ust_fit,zq_fit,'Color',Colors_lit{i},'LineWidth',LineWidth_lit);
end

%organize plot
xlim([0.25 0.6]);
ylim([0 0.15]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('shear velocity, u_{*} (m/s)');
ylabel('saltation height, z_q (m)');
set(gca,'FontSize',PlotFont);


%% PLOT ZQNORM VERSUS UST - PANEL B
subplot(10,1,5:10); hold on;

%plot binned field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>=fQ_moderate);
    errorbar(ust_bins_mid(ind_plot),zqnorm_ust_bin_avg{i}(ind_plot),zqnorm_ust_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot literature data
for i = 1:N_Lit
    plot(ust_lit{i},zqnorm_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end

%fit values for zqnorm versus ust
slope_zqnorm_ust = zeros(N_Sites,1); %slope of fit for zqnorm versus ustar
sigmaslope_zqnorm_ust = zeros(N_Sites,1); %uncertainty in slope for zqnorm versus ustar
zqnorm_bar = zeros(N_Sites,1); %mean of zqnorm
sigma_zqnorm_bar = zeros(N_Sites,1); %uncertainty of zqnorm bar

for i = 1:N_Sites
    %fit values for zqnorm versus ust
    ind_fit = find(fQ_ust_bin_avg{i}>=fQ_moderate);
    ust_fit = ust_bins_mid(ind_fit)';
    sigma_ust_fit = ust_bins_max(ind_fit)'-ust_bins_min(ind_fit)';
    zqnorm_fit = zqnorm_ust_bin_avg{i}(ind_fit);
    sigma_zqnorm_fit = zqnorm_ust_bin_SE{i}(ind_fit);

    [~, b, ~, sigma_b, zqnorm_pred, sigma_zqnorm_pred] = linearfit(ust_bins_mid(ind_fit)', zqnorm_ust_bin_avg{i}(ind_fit), ust_bins_max(ind_fit)'-ust_bins_min(ind_fit)', zqnorm_ust_bin_SE{i}(ind_fit));
    slope_zqnorm_ust(i) = b; %slope of fit
    sigmaslope_zqnorm_ust(i) = sigma_b; %uncertainty in slope of fit
    zqnorm_bar(i) = mean(zqnorm_ust_bin_avg{i}(ind_fit));
    sigma_zqnorm_bar(i) = std(zqnorm_ust_bin_avg{i}(ind_fit))/sqrt(length(ind_fit));
    
    plot(ust_fit, zqnorm_bar(i)*ones(size(ust_fit)),'Color',Colors_field{i},'LineWidth',LineWidth_field); %plot mean value of zqnorm
    %plot(ust_fit,zqnorm_pred,'Color',Colors_field{i},'LineWidth',LineWidth_field); % plot fit value of zqnorm
end

%fit to Greeley, Namikas, and Farrell
for i = 1:N_Lit
    [intercept, slope] = linearfit(ust_lit{i},zqnorm_lit{i});
    ust_fit = [min(ust_lit{i}), max(ust_lit{i})];
    zqnorm_fit = intercept+slope*ust_fit;
    plot(ust_fit, mean(zqnorm_lit{i})*ones(size(ust_fit)),'Color',Colors_lit{i},'LineWidth',LineWidth_lit); %plot mean value of zqnorm
    %plot(ust_fit,zqnorm_fit,'Color',Colors_lit{i},'LineWidth',LineWidth_lit); % plot fit value of zqnorm
end

%organize plot
xlim([0.25 0.6]);
ylim([0 300]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('shear velocity, u_{*} (m/s)');
ylabel('dimensionless salt. ht, z_q/d_{50}');
legend_items = [Sites;LitNames];
legend(legend_items,'Location','SouthOutside');
set(gca,'FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 12]);
print([folder_Plots,'zq_ust.png'],'-dpng');

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.4 7]);
print([folder_Plots,'zq_ust.eps'],'-depsc');


%% FOR EACH SITE, PERFORM BINNING BY TAU
%tau_bins
tau_bins_min = 0.00:0.01:0.39;
tau_bins_max = 0.01:0.01:0.40;
tau_bins_mid = mean([tau_bins_min;tau_bins_max]);
N_tau_bins = length(tau_bins_mid);

%separate fQs into tau bins
fQ_tau_bin_values = cell(N_Sites,1); %fQs into tau bins
N_fQ_tau_bin_values = cell(N_Sites,1); %N fQ's in tau bins
fQ_tau_bin_avg = cell(N_Sites,1); %get average fQ for tau bins
fQ_tau_bin_std = cell(N_Sites,1); %get std fQ for tau bins
fQ_tau_bin_SE = cell(N_Sites,1); %get SE fQ for tau bins

%separate zq's into tau bins
zq_tau_bin_values = cell(N_Sites,1); %zq into tau bins
N_zq_tau_bin_values = cell(N_Sites,1); %N zq's into tau bins
zq_tau_bin_avg = cell(N_Sites,1); %get average zq for tau bins
zq_tau_bin_std = cell(N_Sites,1); %get std zq for tau bins
zq_tau_bin_SE = cell(N_Sites,1); %get SE zq for tau bins

%separate zqnorm's (zq/d50) into tau bins
zqnorm_tau_bin_values = cell(N_Sites,1); %zq/d50s into tau bins
N_zqnorm_tau_bin_values = cell(N_Sites,1); %N zq/d50s in tau bins
zqnorm_tau_bin_avg = cell(N_Sites,1); %get average zq/d50 for tau bins
zqnorm_tau_bin_std = cell(N_Sites,1); %get std zq/d50 for tau bins
zqnorm_tau_bin_SE = cell(N_Sites,1); %get SE zq/d50 for tau bins

%separate Qs into tau bins
Q_tau_bin_values = cell(N_Sites,1); %Qs into tau bins
N_Q_tau_bin_values = cell(N_Sites,1); %N Qs in tau bins
Q_tau_bin_avg = cell(N_Sites,1); %get average Q for tau bins
Q_tau_bin_std = cell(N_Sites,1); %get std Q for tau bins
Q_tau_bin_SE = cell(N_Sites,1); %get SE Q for tau bins

%go through all sites
for i = 1:N_Sites
  
    %initialize fQ's into tau bins
    fQ_tau_bin_values{i} = cell(N_tau_bins,1);
    N_fQ_tau_bin_values{i} = zeros(N_tau_bins,1);
    fQ_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    fQ_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    fQ_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;
    
    %initialize zq's into tau bins
    zq_tau_bin_values{i} = cell(N_tau_bins,1);
    N_zq_tau_bin_values{i} = zeros(N_tau_bins,1);
    zq_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    zq_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    zq_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;
    
    %initialize zqnorm's into tau bins
    zqnorm_tau_bin_values{i} = cell(N_tau_bins,1);
    N_zqnorm_tau_bin_values{i} = zeros(N_tau_bins,1);
    zqnorm_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    zqnorm_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    zqnorm_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;
    
    %initialize Q's into tau bins
    Q_tau_bin_values{i} = cell(N_tau_bins,1);
    N_Q_tau_bin_values{i} = zeros(N_tau_bins,1);
    Q_tau_bin_avg{i} = zeros(N_tau_bins,1)*NaN;
    Q_tau_bin_std{i} = zeros(N_tau_bins,1)*NaN;
    Q_tau_bin_SE{i} = zeros(N_tau_bins,1)*NaN;

    %get values, avgs, and standard deviations for tau bins
    for j = 1:N_tau_bins
        %get indices
        bin_ind = find(tauRe_all{i}>=tau_bins_min(j)&tauRe_all{i}<=tau_bins_max(j));

        %get fQs for tau bins
        fQ_tau_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_tau_bin_values{i}{j} = fQ_tau_bin_values{i}{j}(~isnan(fQ_tau_bin_values{i}{j}));
        N_fQ_tau_bin_values{i}(j) = length(fQ_tau_bin_values{i}{j});
        fQ_tau_bin_avg{i}(j) = mean(fQ_tau_bin_values{i}{j});
        fQ_tau_bin_std{i}(j) = std(fQ_tau_bin_values{i}{j});
        
        %get zqs for tau bins
        zq_tau_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_tau_bin_values{i}{j} = zq_tau_bin_values{i}{j}(~isnan(zq_tau_bin_values{i}{j}));
        N_zq_tau_bin_values{i}(j) = length(zq_tau_bin_values{i}{j});
        zq_tau_bin_avg{i}(j) = mean(zq_tau_bin_values{i}{j});
        zq_tau_bin_std{i}(j) = std(zq_tau_bin_values{i}{j});
        
        %get zqnorms for tau bins
        zqnorm_tau_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_tau_bin_values{i}{j} = zqnorm_tau_bin_values{i}{j}(~isnan(zqnorm_tau_bin_values{i}{j}));
        N_zqnorm_tau_bin_values{i}(j) = length(zqnorm_tau_bin_values{i}{j});
        zqnorm_tau_bin_avg{i}(j) = mean(zqnorm_tau_bin_values{i}{j});
        zqnorm_tau_bin_std{i}(j) = std(zqnorm_tau_bin_values{i}{j});
        
        %get Qs for tau bins
        Q_tau_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_tau_bin_values{i}{j} = Q_tau_bin_values{i}{j}(~isnan(Q_tau_bin_values{i}{j}));
        N_Q_tau_bin_values{i}(j) = length(Q_tau_bin_values{i}{j});
        Q_tau_bin_avg{i}(j) = mean(Q_tau_bin_values{i}{j});
        Q_tau_bin_std{i}(j) = std(Q_tau_bin_values{i}{j});
    end

    %get indices for bins with multiple values
    ind_multibin_fQ_tau = find(N_fQ_tau_bin_values{i}>1);
    ind_multibin_zq_tau = find(N_zq_tau_bin_values{i}>1);
    ind_multibin_zqnorm_tau = find(N_zqnorm_tau_bin_values{i}>1);
    ind_multibin_Q_tau = find(N_Q_tau_bin_values{i}>1);
    
    %get standard errors for tau bins
    for j = 1:N_tau_bins
        %get SEs for fQ
        if N_fQ_tau_bin_values{i}(j)>1
            fQ_tau_bin_SE{i}(j) = fQ_tau_bin_std{i}(j)/sqrt(N_fQ_tau_bin_values{i}(j));
        elseif N_fQ_tau_bin_values{i}(j)==1
            ind_below = ind_multibin_fQ_tau(find(ind_multibin_fQ_tau<j,1,'last'));
            ind_above = ind_multibin_fQ_tau(find(ind_multibin_fQ_tau>j,1));
            fQ_tau_bin_SE{i}(j) = mean(fQ_tau_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for zq
        if N_zq_tau_bin_values{i}(j)>1
            zq_tau_bin_SE{i}(j) = zq_tau_bin_std{i}(j)/sqrt(N_zq_tau_bin_values{i}(j));
        elseif N_zq_tau_bin_values{i}(j)==1
            ind_below = ind_multibin_zq_tau(find(ind_multibin_zq_tau<j,1,'last'));
            ind_above = ind_multibin_zq_tau(find(ind_multibin_zq_tau>j,1));
            zq_tau_bin_SE{i}(j) = mean(zq_tau_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for zqnorm
        if N_zqnorm_tau_bin_values{i}(j)>1
            zqnorm_tau_bin_SE{i}(j) = zqnorm_tau_bin_std{i}(j)/sqrt(N_zqnorm_tau_bin_values{i}(j));
        elseif N_zqnorm_tau_bin_values{i}(j)==1
            ind_below = ind_multibin_zqnorm_tau(find(ind_multibin_zqnorm_tau<j,1,'last'));
            ind_above = ind_multibin_zqnorm_tau(find(ind_multibin_zqnorm_tau>j,1));
            zqnorm_tau_bin_SE{i}(j) = mean(zqnorm_tau_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for Q
        if N_Q_tau_bin_values{i}(j)>1
            Q_tau_bin_SE{i}(j) = Q_tau_bin_std{i}(j)/sqrt(N_Q_tau_bin_values{i}(j));
        elseif N_Q_tau_bin_values{i}(j)==1
            ind_below = ind_multibin_Q_tau(find(ind_multibin_Q_tau<j,1,'last'));
            ind_above = ind_multibin_Q_tau(find(ind_multibin_Q_tau>j,1));
            Q_tau_bin_SE{i}(j) = mean(Q_tau_bin_std{i}([ind_below ind_above]));
        end
    end
end


%% PERFORM FIT OF Q VERSUS TAU
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

for i = 1:N_Sites
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
end

%% PLOT Q VERSUS TAU
figure; clf; hold on;
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    errorbar(tau_bins_mid,Q_tau_bin_avg{i},Q_tau_bin_SE{i},Markers{i},'MarkerSize',10);
    plot([tauit_pred_transport(i) max(tau_transport{i})],[0 CQ_pred_transport(i)*(max(tau_transport{i})-tauit_pred_transport(i))],LineColors{i}); %plot fit for these values
    legend_items{2*i-1} = Sites{i}; %add to list of legend items
    legend_items{2*i} = 'fit'; %add to list of legend items
end
xlabel('\tau (Pa)');
ylabel('Q (g/m/s)');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',PlotFont);
set(gcf, 'PaperPosition',[0 0 8 5]);
print([folder_Plots,'flux_tau.png'],'-dpng');

    
%% FOR EACH SITE, PERFORM BINNING BY TAUEX
%compute tauex based on fitted thresholds
tauex_all = cell(N_Sites,1);
for i=1:N_Sites
    tauex_all{i} = tauRe_all{i}-tauit_pred_transport(i);
end

%compute Qnorm
Qnorm_all = cell(N_Sites,1);
for i=1:N_Sites
    Qnorm_all{i} = Q_all{i}./((1/g)*sqrt(tauit_pred_transport(i)/rho_a).*1e3*tauex_all{i});
end

%tauex_bins
tauex_bins_min = -0.10:0.01:0.29; %small bins
tauex_bins_max = -0.09:0.01:0.30; %small bins
tauex_bins_mid = mean([tauex_bins_min;tauex_bins_max]);
N_tauex_bins = length(tauex_bins_mid);

%separate fQs into tauex bins
fQ_tauex_bin_values = cell(N_Sites,1); %fQs into tauex bins
N_fQ_tauex_bin_values = cell(N_Sites,1); %N fQ's in tauex bins
fQ_tauex_bin_avg = cell(N_Sites,1); %get average fQ for tauex bins
fQ_tauex_bin_std = cell(N_Sites,1); %get std fQ for tauex bins
fQ_tauex_bin_SE = cell(N_Sites,1); %get SE fQ for tauex bins

%separate zq's into tauex bins
zq_tauex_bin_values = cell(N_Sites,1); %zq into tauex bins
N_zq_tauex_bin_values = cell(N_Sites,1); %N zq's into tauex bins
zq_tauex_bin_avg = cell(N_Sites,1); %get average zq for tauex bins
zq_tauex_bin_std = cell(N_Sites,1); %get std zq for tauex bins
zq_tauex_bin_SE = cell(N_Sites,1); %get SE zq for tauex bins

%separate zqnorm's (zq/d50) into tauex bins
zqnorm_tauex_bin_values = cell(N_Sites,1); %zq/d50s into tauex bins
N_zqnorm_tauex_bin_values = cell(N_Sites,1); %N zq/d50s in tauex bins
zqnorm_tauex_bin_avg = cell(N_Sites,1); %get average zq/d50 for tauex bins
zqnorm_tauex_bin_std = cell(N_Sites,1); %get std zq/d50 for tauex bins
zqnorm_tauex_bin_SE = cell(N_Sites,1); %get SE zq/d50 for tauex bins

%separate Qs into tauex bins
Q_tauex_bin_values = cell(N_Sites,1); %Qs into tauex bins
N_Q_tauex_bin_values = cell(N_Sites,1); %N Qs in tauex bins
Q_tauex_bin_avg = cell(N_Sites,1); %get average Q for tauex bins
Q_tauex_bin_std = cell(N_Sites,1); %get std Q for tauex bins
Q_tauex_bin_SE = cell(N_Sites,1); %get SE Q for tauex bins

%separate Qnorms into tauex bins
Qnorm_tauex_bin_values = cell(N_Sites,1); %Qs into tauex bins
N_Qnorm_tauex_bin_values = cell(N_Sites,1); %N Qs in tauex bins
Qnorm_tauex_bin_avg = cell(N_Sites,1); %get average Q for tauex bins
Qnorm_tauex_bin_std = cell(N_Sites,1); %get std Q for tauex bins
Qnorm_tauex_bin_SE = cell(N_Sites,1); %get SE Q for tauex bins

%go through all sites
for i = 1:N_Sites
  
    %initialize fQ's into tauex bins
    fQ_tauex_bin_values{i} = cell(N_tauex_bins,1);
    N_fQ_tauex_bin_values{i} = zeros(N_tauex_bins,1);
    fQ_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    fQ_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    fQ_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize zq's into tauex bins
    zq_tauex_bin_values{i} = cell(N_tauex_bins,1);
    N_zq_tauex_bin_values{i} = zeros(N_tauex_bins,1);
    zq_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    zq_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    zq_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize zqnorm's into tauex bins
    zqnorm_tauex_bin_values{i} = cell(N_tauex_bins,1);
    N_zqnorm_tauex_bin_values{i} = zeros(N_tauex_bins,1);
    zqnorm_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    zqnorm_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    zqnorm_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %initialize Q's into tauex bins
    Q_tauex_bin_values{i} = cell(N_tauex_bins,1);
    N_Q_tauex_bin_values{i} = zeros(N_tauex_bins,1);
    Q_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    Q_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    Q_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;
    
    %separate Qnorms into tauex bins
    Qnorm_tauex_bin_values{i} = cell(N_tauex_bins,1);
    N_Qnorm_tauex_bin_values{i} = zeros(N_tauex_bins,1);
    Qnorm_tauex_bin_avg{i} = zeros(N_tauex_bins,1)*NaN;
    Qnorm_tauex_bin_std{i} = zeros(N_tauex_bins,1)*NaN;
    Qnorm_tauex_bin_SE{i} = zeros(N_tauex_bins,1)*NaN;

    %get values, avgs, and standard deviations for tauex bins
    for j = 1:N_tauex_bins
        %get indices
        bin_ind = find(tauex_all{i}>=tauex_bins_min(j)&tauex_all{i}<=tauex_bins_max(j));

        %get fQs for tauex bins
        fQ_tauex_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_tauex_bin_values{i}{j} = fQ_tauex_bin_values{i}{j}(~isnan(fQ_tauex_bin_values{i}{j}));
        N_fQ_tauex_bin_values{i}(j) = length(fQ_tauex_bin_values{i}{j});
        fQ_tauex_bin_avg{i}(j) = mean(fQ_tauex_bin_values{i}{j});
        fQ_tauex_bin_std{i}(j) = std(fQ_tauex_bin_values{i}{j});
        
        %get zqs for tauex bins
        zq_tauex_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_tauex_bin_values{i}{j} = zq_tauex_bin_values{i}{j}(~isnan(zq_tauex_bin_values{i}{j}));
        N_zq_tauex_bin_values{i}(j) = length(zq_tauex_bin_values{i}{j});
        zq_tauex_bin_avg{i}(j) = mean(zq_tauex_bin_values{i}{j});
        zq_tauex_bin_std{i}(j) = std(zq_tauex_bin_values{i}{j});
        
        %get zqnorms for tauex bins
        zqnorm_tauex_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_tauex_bin_values{i}{j} = zqnorm_tauex_bin_values{i}{j}(~isnan(zqnorm_tauex_bin_values{i}{j}));
        N_zqnorm_tauex_bin_values{i}(j) = length(zqnorm_tauex_bin_values{i}{j});
        zqnorm_tauex_bin_avg{i}(j) = mean(zqnorm_tauex_bin_values{i}{j});
        zqnorm_tauex_bin_std{i}(j) = std(zqnorm_tauex_bin_values{i}{j});
        
        %get Qs for tauex bins
        Q_tauex_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_tauex_bin_values{i}{j} = Q_tauex_bin_values{i}{j}(~isnan(Q_tauex_bin_values{i}{j}));
        N_Q_tauex_bin_values{i}(j) = length(Q_tauex_bin_values{i}{j});
        Q_tauex_bin_avg{i}(j) = mean(Q_tauex_bin_values{i}{j});
        Q_tauex_bin_std{i}(j) = std(Q_tauex_bin_values{i}{j});
        
        %get Qnorms for tauex bins
        Qnorm_tauex_bin_values{i}{j} = Qnorm_all{i}(bin_ind);
        Qnorm_tauex_bin_values{i}{j} = Qnorm_tauex_bin_values{i}{j}(~isnan(Qnorm_tauex_bin_values{i}{j}));
        N_Qnorm_tauex_bin_values{i}(j) = length(Qnorm_tauex_bin_values{i}{j});
        Qnorm_tauex_bin_avg{i}(j) = mean(Qnorm_tauex_bin_values{i}{j});
        Qnorm_tauex_bin_std{i}(j) = std(Qnorm_tauex_bin_values{i}{j});
    end

    %get indices for bins with multiple values
    ind_multibin_fQ_tauex = find(N_fQ_tauex_bin_values{i}>1);
    ind_multibin_zq_tauex = find(N_zq_tauex_bin_values{i}>1);
    ind_multibin_zqnorm_tauex = find(N_zqnorm_tauex_bin_values{i}>1);
    ind_multibin_Q_tauex = find(N_Q_tauex_bin_values{i}>1);
    ind_multibin_Qnorm_tauex = find(N_Qnorm_tauex_bin_values{i}>1);
    
    %get standard errors for tauex bins
    for j = 1:N_tauex_bins
        %get SEs for fQ
        if N_fQ_tauex_bin_values{i}(j)>1
            fQ_tauex_bin_SE{i}(j) = fQ_tauex_bin_std{i}(j)/sqrt(N_fQ_tauex_bin_values{i}(j));
        elseif N_fQ_tauex_bin_values{i}(j)==1
            ind_below = ind_multibin_fQ_tauex(find(ind_multibin_fQ_tauex<j,1,'last'));
            ind_above = ind_multibin_fQ_tauex(find(ind_multibin_fQ_tauex>j,1));
            fQ_tauex_bin_SE{i}(j) = mean(fQ_tauex_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for zq
        if N_zq_tauex_bin_values{i}(j)>1
            zq_tauex_bin_SE{i}(j) = zq_tauex_bin_std{i}(j)/sqrt(N_zq_tauex_bin_values{i}(j));
        elseif N_zq_tauex_bin_values{i}(j)==1
            ind_below = ind_multibin_zq_tauex(find(ind_multibin_zq_tauex<j,1,'last'));
            ind_above = ind_multibin_zq_tauex(find(ind_multibin_zq_tauex>j,1));
            zq_tauex_bin_SE{i}(j) = mean(zq_tauex_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for zqnorm
        if N_zqnorm_tauex_bin_values{i}(j)>1
            zqnorm_tauex_bin_SE{i}(j) = zqnorm_tauex_bin_std{i}(j)/sqrt(N_zqnorm_tauex_bin_values{i}(j));
        elseif N_zqnorm_tauex_bin_values{i}(j)==1
            ind_below = ind_multibin_zqnorm_tauex(find(ind_multibin_zqnorm_tauex<j,1,'last'));
            ind_above = ind_multibin_zqnorm_tauex(find(ind_multibin_zqnorm_tauex>j,1));
            zqnorm_tauex_bin_SE{i}(j) = mean(zqnorm_tauex_bin_std{i}([ind_below ind_above]));
        end
        
        %get SEs for Q
        if N_Q_tauex_bin_values{i}(j)>1
            Q_tauex_bin_SE{i}(j) = Q_tauex_bin_std{i}(j)/sqrt(N_Q_tauex_bin_values{i}(j));
        elseif N_Q_tauex_bin_values{i}(j)==1
            ind_below = ind_multibin_Q_tauex(find(ind_multibin_Q_tauex<j,1,'last'));
            ind_above = ind_multibin_Q_tauex(find(ind_multibin_Q_tauex>j,1));
            Q_tauex_bin_SE{i}(j) = mean(Q_tauex_bin_std{i}([ind_below ind_above]));
        end

        %get SEs for Qnorm
        if N_Qnorm_tauex_bin_values{i}(j)>1
            Qnorm_tauex_bin_SE{i}(j) = Qnorm_tauex_bin_std{i}(j)/sqrt(N_Qnorm_tauex_bin_values{i}(j));
        elseif N_Qnorm_tauex_bin_values{i}(j)==1
            ind_below = ind_multibin_Qnorm_tauex(find(ind_multibin_Qnorm_tauex<j,1,'last'));
            ind_above = ind_multibin_Qnorm_tauex(find(ind_multibin_Qnorm_tauex>j,1));
            Qnorm_tauex_bin_SE{i}(j) = mean(Qnorm_tauex_bin_std{i}([ind_below ind_above]));
        end
    end
end

%% FITS TO Q AND QNORM VERSUS TAUEX

%scaling coefficient for Q versus tauex
C_Q_pred = zeros(N_Sites,1); 
sigma_C_Q_pred = zeros(N_Sites,1);

%scaling coefficient for Qnorm versus tauex
C_Qnorm_pred = zeros(N_Sites,1);
sigma_C_Qnorm_pred = zeros(N_Sites,1);

%perform fits
for i = 1:N_Sites
    tauex_bin_ind = find(fQ_tauex_bin_avg{i}>=fQ_continuous); %get tau bins that are continuous
    tauex_continuous = tauex_bins_mid(tauex_bin_ind)';
    Q_continuous = Q_tauex_bin_avg{i}(tauex_bin_ind);
    Qnorm_continuous = Qnorm_tauex_bin_avg{i}(tauex_bin_ind);
    C_Q_pred(i) = mean(Q_continuous./tauex_continuous);
    sigma_C_Q_pred(i) = std(Q_continuous./tauex_continuous)/length(tauex_continuous);
    C_Qnorm_pred(i) = mean(Qnorm_continuous);
    sigma_C_Qnorm_pred(i) = std(Qnorm_continuous)/length(tauex_continuous);
end


%% PLOT Q VERSUS TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    errorbar(tauex_bins_mid,Q_tauex_bin_avg{i},Q_tauex_bin_SE{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit, C_Q_pred(i)*tauex_fit,'Color',Colors_field{i});
end

xlim([0 0.3]);
ylim([0 60]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('excess shear stress, \tau_{ex} (Pa)');
ylabel('saltation mass flux, Q (g m^{-1} s^{-1})');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',PlotFont);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 4]);
print([folder_Plots,'Q_tauex.png'],'-dpng');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 4]);
print([folder_Plots,'Q_tauex.eps'],'-depsc');


%% PLOT QNORM VS TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    ind_plot = find(fQ_tauex_bin_avg{i}>=fQ_continuous);
    errorbar(tauex_bins_mid(ind_plot),Qnorm_tauex_bin_avg{i}(ind_plot),Qnorm_tauex_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0, max(tauex_bins_max)];
for i = 1:N_Sites
    plot(tauex_fit,C_Qnorm_pred(i)*ones(2,1),'Color',Colors_field{i});
end
xlabel('excess shear stress, $$\tau_{ex}$$ (Pa)','Interpreter','Latex');
ylabel('dimensionless saltation flux, $$\hat{Q}$$','Interpreter','Latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
legend(Sites,'Location','SouthEast');
xlim([0 max(tauex_bins_max)]);
ylim([0 10]);
set(gca,'FontSize',PlotFont);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 7]);
print([folder_Plots,'Qnorm_tauex.png'],'-dpng');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.4 3]);
print([folder_Plots,'Qnorm_tauex.eps'],'-depsc');

% % %% FOR EACH SITE, PERFORM BINNING BY fQ BINS
% % 
% % %separate tauth based on TFEM into frequency bins
% % tauth_fQ_altbin_values = cell(N_Sites,1);
% % tauth_fQ_altbin_avg = cell(N_Sites,1);
% % tauth_fQ_altbin_std = cell(N_Sites,1);
% % tauth_fQ_altbin_SE = cell(N_Sites,1);
% % 
% % for i = 1:N_Sites
% %     tauth_fQ_altbin_values{i} = cell(N_fQ_bins,1);
% %     tauth_fQ_altbin_avg{i} = zeros(N_fQ_bins,1)*NaN;
% %     tauth_fQ_altbin_std{i} = zeros(N_fQ_bins,1)*NaN;
% %     tauth_fQ_altbin_SE{i} = zeros(N_fQ_bins,1)*NaN;
% %     
% %     %get tauth's for fQ bins
% %     for j = 1:N_fQ_bins
% %         %get indices for bins
% %         bin_ind = find(fQ1_all{i}>=fQ_bins_min(j)&fQ1_all{i}<=fQ_bins_max(j)); %use 1 second frequencies for binning
% %                 
% %         %tau_th
% %         tauth_fQ_altbin_values{i}{j} = tauth_TFEM_all{i}(bin_ind);
% %         tauth_fQ_altbin_values{i}{j} = tauth_fQ_altbin_values{i}{j}(~isnan(tauth_fQ_altbin_values{i}{j}));
% %         tauth_fQ_altbin_avg{i}(j) = mean(tauth_fQ_altbin_values{i}{j});
% %         tauth_fQ_altbin_std{i}(j) = std(tauth_fQ_altbin_values{i}{j});
% %         tauth_fQ_altbin_SE{i}(j) = std(tauth_fQ_altbin_values{i}{j})/sqrt(length(tauth_fQ_altbin_values{i}{j}));
% %     end 
% % end
% % 
% % %% PERFORM FITS FOR FQ
% %
% % %fit to TFEM thresholds
% % fQ_TFEM_ind = cell(N_Sites,1);
% % tauft_TFEM = zeros(N_Sites,1);
% % tauit_TFEM = zeros(N_Sites,1);
% % tauth_fit_TFEM = cell(N_Sites,1);
% % sigma_tauth_fit_TFEM = cell(N_Sites,1);
% % 
% % for i = 1:N_Sites
% %     fQ_TFEM_ind{i} = intersect(intersect(find(fQ_bins_mid>fQ_TFEM_fit_min(i)),find(fQ_bins_mid<fQ_TFEM_fit_max(i))),find(~isnan(tauth_fQ_altbin_avg{i})));
% %     if i==3
% %         [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_bins_mid(fQ_TFEM_ind{i})', tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}), fQ_bins_max(fQ_TFEM_ind{i})'-fQ_bins_min(fQ_TFEM_ind{i})', tauth_fQ_altbin_std{i}(fQ_TFEM_ind{i}));
% %     else %don't use confidence intervals for Jeri and RG, because there is insufficient data for these
% %         [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_bins_mid(fQ_TFEM_ind{i})', tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}));
% %     end
% %     tauft_TFEM(i) = a;
% %     tauit_TFEM(i) = a+b;
% %     tauth_fit_TFEM{i} = tauth_fit;
% %     sigma_tauth_fit_TFEM{i} = sigma_tauth_fit;
% % end 
% %
% % %% PLOT TAUTH VERSUS FQ
% % figure; clf; hold on; %initialize plot
% % %legend_items = cell(N_Sites*2,1);
% % for i = 1:N_Sites
% %     errorbar(fQ_bins_mid(fQ_TFEM_ind{i}),tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}),tauth_fQ_altbin_SE{i}(fQ_TFEM_ind{i}),Markers{i},'MarkerSize',10);
% %     %legend_items{2*i-1} = Sites{i}; %add to list of legend items
% %     %legend_items{2*i} = 'fit'; %add to list of legend items
% % end
% % for i = 1:N_Sites
% %     plot([0 1],[tauft_TFEM(i) tauit_TFEM(i)],LineColors{i});
% % end
% % xlabel('frequency of transport');
% % ylabel('TFEM inferred threshold stress, \tau_{th}');
% % legend(Sites,'Location','NorthEast');
% % set(gca,'FontSize',PlotFont);
% % set(gcf, 'PaperPosition',[0 0 10 6]);
% % print([folder_Plots,'tauth_fQ_all.png'],'-dpng');