%% ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;
close all;

%parameters
rho_a = 1.18; %air density (kg/m^3)
g = 9.8; %gravity (m/s^2)
z0 = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m) at threshold
theta_limit = 20; %highest angle (deg) for Oceano observations
zL_limit = -0.2; %limit for stability parameter (z/L) for all observations
zW_limit = 4; %limit on number of unique Wenglor heights in profile
zq_limit = 1; %upper limit of saltation height (m) for all observations

%fitting parameters
fQ_zQ_ust = 0; %mininum fQ for zQ versus u* comparison (for constant flux height)
fQ_Q_tau = 0.1; %minimum mean fQ for Q versus tau comparison (and determination of threshold)
fQ_Q_tauex = 0; %minimum fQ for Q versus tauex comparison
N_sigma_tauit = 2; %std deviations from impact threshold for plotting

%binning information
bin_N_min = 3; %minimum number of entries for bin
ust_bin_lower = 0;
ust_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
ust_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
tau_bin_lower = 0;
tau_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
tauex_bin_lower = 0;
tauex_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tauex_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin

%information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for storing data output
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
folder_Plots = '../../PlotOutput/StressFlux/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%load grain size data
load(strcat(folder_GrainSizeData,'MeanGrainSize'));
d50_site = d50_surface_site;

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows'));
N_Sites = length(Sites);

%load external data
load(strcat(folder_AnalysisData,'LitData')); %literature data
LitNames = {'Greeley et al. (1996)'; 'Namikas (2003)';'Farrell et al. (2012)'};
ust_lit = {ust_Greeley96, ust_Namikas03, ust_Farrell12};
sigma_ust_lit = {sigma_ust_Greeley96, sigma_ust_Namikas03, sigma_ust_Farrell12};
zq_lit = {zq_Greeley96, zq_Namikas03, zq_Farrell12};
sigma_zq_lit = {sigma_zq_Greeley96, sigma_zq_Namikas03, sigma_zq_Farrell12};
zqnorm_lit = {1e3*zq_Greeley96/d50_Greeley96, 1e3*zq_Namikas03/d50_Namikas03, 1e3*zq_Farrell12/d50_Farrell12};
sigma_zqnorm_lit = {1e3*sigma_zq_Greeley96/d50_Greeley96, 1e3*sigma_zq_Namikas03/d50_Namikas03, 1e3*sigma_zq_Farrell12/d50_Farrell12};
N_lit = 3;

%change zqnorm_lit so that Farrell doesn't show up 
zqnorm_lit{3} = zqnorm_lit{3}*NaN;

%set info for plotting
Markers_field = {'s','d','o','<','>'};
Colors_field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};
MarkerSize_field = 8;
LineWidth_field = 1;
Markers_lit = {'x','+','*'};
MarkerSize_lit = 8;
Colors_lit = {[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
LineWidth_lit = 1;
PlotFont = 14;


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

% remove points with too few Wenglor heights
for i = 1:N_Sites
    N_zW = zeros(size(zW_all{i}));
    for j = 1:length(zW_all{i})
        N_zW(j) = length(unique(zW_all{i}{j}));
    end
    ind_zW = find(N_zW>=zW_limit);
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_zW);']);
    end
end

% remove points with unreasonably large zq
for i = 1:N_Sites
    ind_zq = find(zq_all{i}<=zq_limit);
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_zq);']);
    end
end


% %% SEPARATE DATA BY CLUSTERS - EXPAND TO FIVE "SITES"
% %rename variables to '_Oceano'
% for j = 1:N_variables
%     eval(strcat(variable_list{j}(1:end-4),'_Oceano=',variable_list{j},'{3};'));
% end
% 
% variable_list_Oceano = who('variables','*Oceano');
% N_clusters = length(cluster_StartDate);
% %now separate by date
% for i = 1:N_clusters
%     ind_cluster = intersect(find(Date_Oceano>=cluster_StartDate(i)),find(Date_Oceano<=cluster_EndDate(i)));
%     for j = 1:N_variables
%         instruction = strcat(variable_list(j),'{2+i}=',variable_list_Oceano(j),'(ind_cluster);');
%         eval(instruction{1});
%     end
%     d50_site(2+i) = d50_surface_cluster(i);
%     Sites{2+i} = ['Oceano_',int2str(i)];
% end
% N_Sites = length(Sites);


%% COMPUTE ZQNORM BASED ON GRAIN SIZE BY SITE
zqnorm_all = cell(N_Sites,1);
sigma_zqnorm_all = cell(N_Sites,1);
for i=1:N_Sites
    zqnorm_all{i} = 1000*zq_all{i}./d50_site(i);
    sigma_zqnorm_all{i} = 1000*sigma_zq_all{i}./d50_site(i);
end


%% FOR EACH SITE, PERFORM BINNING BY USTAR
%initialize lists of u* bins by site
ust_bin_values_all = cell(N_Sites,1); %ust's in each bin
ust_bin_N_all = cell(N_Sites,1); %number of ust in bin
ust_bin_avg_all = cell(N_Sites,1); %avg ust in bin
ust_bin_SE_all = cell(N_Sites,1); %SE ust in bin

%separate fQs into ust bins
fQ_ust_bin_values_all = cell(N_Sites,1); %fQs into ust bins
fQ_ust_bin_avg_all = cell(N_Sites,1); %get average fQ for ust bins
fQ_ust_bin_std_all = cell(N_Sites,1); %get std fQ for ust bins
fQ_ust_bin_SE_all = cell(N_Sites,1); %get SE fQ for ust bins

%separate zq's into ust bins
zq_ust_bin_values_all = cell(N_Sites,1); %zq into ust bins
zq_ust_bin_avg_all = cell(N_Sites,1); %get average zq for ust bins
zq_ust_bin_std_all = cell(N_Sites,1); %get std zq for ust bins
zq_ust_bin_SE_all = cell(N_Sites,1); %get SE zq for ust bins

%separate zqnorm's (zq/d50) into ust bins
zqnorm_ust_bin_values_all = cell(N_Sites,1); %zq/d50s into ust bins
zqnorm_ust_bin_avg_all = cell(N_Sites,1); %get average zq/d50 for ust bins
zqnorm_ust_bin_std_all = cell(N_Sites,1); %get std zq/d50 for ust bins
zqnorm_ust_bin_SE_all = cell(N_Sites,1); %get SE zq/d50 for ust bins

%get binned values for all sites
for i = 1:N_Sites
    
    %use only values above minimum u* and Q
    ind_analysis = find((ustRe_all{i}>ust_bin_lower)&(Q_all{i}>0));
    ust = ustRe_all{i}(ind_analysis);
    fQ = fQ_all{i}(ind_analysis);
    zq = zq_all{i}(ind_analysis);
    zqnorm = zqnorm_all{i}(ind_analysis);
    
    %get u* bins
    [ust_bin_values, ust_bin_N, ust_bin_min, ust_bin_max, ust_bin_avg, ust_bin_SE] = ...
        Binning(ust, ust_bin_minrange, ust_bin_maxrange, bin_N_min);
    ust_bin_values_all{i} = ust_bin_values; %ust's in each bin
    ust_bin_N_all{i} = ust_bin_N; %number of ust in bin
    ust_bin_avg_all{i} = ust_bin_avg; %avg ust in bin
    ust_bin_SE_all{i} = ust_bin_SE; %SE ust in bin
    
    %get values in bins
    N_ust_bins = length(ust_bin_values); %number of bins
    fQ_ust_bin_values_all{i} = cell(N_ust_bins,1); %initialize fQ's
    zq_ust_bin_values_all{i} = cell(N_ust_bins,1); %initialize zq's
    zqnorm_ust_bin_values_all{i} = cell(N_ust_bins,1); %initialize zqnorm's
    for j = 1:N_ust_bins
        %get indices of values in bin
        bin_ind = find(ust>=ust_bin_min(j)&ust<=ust_bin_max(j));
        
        %get values in bins
        fQ_ust_bin_values_all{i}{j} = fQ(bin_ind); %fQ's
        zq_ust_bin_values_all{i}{j} = zq(bin_ind); %zq's
        zqnorm_ust_bin_values_all{i}{j} = zqnorm(bin_ind); %zqnorm's
    end
    
    %get fQ avg, std, and SE
    fQ_ust_bin_avg_all{i} = cellfun(@mean,fQ_ust_bin_values_all{i});
    fQ_ust_bin_std_all{i} = cellfun(@std,fQ_ust_bin_values_all{i});
    fQ_ust_bin_SE_all{i} = fQ_ust_bin_std_all{i}./cellfun(@length,fQ_ust_bin_values_all{i});
    
    %get zq avg, std, and SE
    zq_ust_bin_avg_all{i} = cellfun(@mean,zq_ust_bin_values_all{i});
    zq_ust_bin_std_all{i} = cellfun(@std,zq_ust_bin_values_all{i});
    zq_ust_bin_SE_all{i} = zq_ust_bin_std_all{i}./cellfun(@length,zq_ust_bin_values_all{i});
    
    %get zq avg, std, and SE
    zqnorm_ust_bin_avg_all{i} = cellfun(@mean,zqnorm_ust_bin_values_all{i});
    zqnorm_ust_bin_std_all{i} = cellfun(@std,zqnorm_ust_bin_values_all{i});
    zqnorm_ust_bin_SE_all{i} = zqnorm_ust_bin_std_all{i}./cellfun(@length,zqnorm_ust_bin_values_all{i});
end

%recompute standard errors for bins with fewer than bin_N_min values based on nearest Oceano u*
i_Oceano = find(strcmp(Sites,'Oceano'));
ind_SE_Oceano = find(ust_bin_N_all{i_Oceano}>=bin_N_min); %indices of Oceano bins with sufficient values for SE
ust_avg_Oceano = ust_bin_avg_all{i_Oceano}(ind_SE_Oceano); %ust values for these bins
fQ_CV_Oceano = fQ_ust_bin_std_all{i_Oceano}(ind_SE_Oceano)./fQ_ust_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of fQ values for these bins
zq_CV_Oceano = zq_ust_bin_std_all{i_Oceano}(ind_SE_Oceano)./zq_ust_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of zq values for these bins
zqnorm_CV_Oceano = zqnorm_ust_bin_std_all{i_Oceano}(ind_SE_Oceano)./zqnorm_ust_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of zqnorm values for these bins
for i = 1:N_Sites
    ind_recompute_SE = find(ust_bin_N_all{i}<bin_N_min);
    N_ind = length(ind_recompute_SE);
    for j = 1:N_ind
        k = ind_recompute_SE(j);
        ust_bin = ust_bin_avg_all{i}(k); %get avg ust for bin
        ust_diff_Oceano = abs(ust_bin-ust_avg_Oceano); %get distance of this ust from ust's for Oceano
        ind_Oceano = find(ust_diff_Oceano==min(ust_diff_Oceano),1); %get closest ust for Oceano
        fQ_std = fQ_CV_Oceano(ind_Oceano)*fQ_ust_bin_avg_all{i}(k); %std dev fQ is mean fQ value from site times CV for closest u* at Oceano
        fQ_ust_bin_SE_all{i}(k) = fQ_std/sqrt(ust_bin_N_all{i}(k)); %use this std dev to get SE
        zq_std = zq_CV_Oceano(ind_Oceano)*zq_ust_bin_avg_all{i}(k); %std dev zq is mean zq value from site times CV for closest u* at Oceano
        zq_ust_bin_SE_all{i}(k) = zq_std/sqrt(ust_bin_N_all{i}(k)); %use this std dev to get SE
        zqnorm_std = zqnorm_CV_Oceano(ind_Oceano)*zqnorm_ust_bin_avg_all{i}(k); %std dev zqnorm is mean zqnorm value from site times CV for closest u* at Oceano
        zqnorm_ust_bin_SE_all{i}(k) = zqnorm_std/sqrt(ust_bin_N_all{i}(k)); %use this std dev to get SE
    end
end

%% FITTING ZQ VS U*
%compute fit values for zq and zqnorm versus ust - field
ust_fit_all = cell(N_Sites,1); %ust's for fit

intercept_zq_ust_all = zeros(N_Sites,1); %intercept of fit
sigmaintercept_zq_ust_all = zeros(N_Sites,1); %uncertainty in intercept of fit
slope_zq_ust_all = zeros(N_Sites,1); %slope of fit for zq versus ustar
sigmaslope_zq_ust_all = zeros(N_Sites,1); %uncertainty in slope for zq versus ustar
zq_pred_all{i} = cell(N_Sites,1); %predicted zq

zqnorm_pred_all{i} = cell(N_Sites,1); %predicted zqnorm
zqnorm_bar_all = zeros(N_Sites,1); %mean of zqnorm
sigma_zqnorm_bar_all = zeros(N_Sites,1); %uncertainty of zqnorm bar
zqnorm_pred_relativerange_all = zeros(N_Sites,1); %relative range of zqnorm prediction

for i = 1:N_Sites
    ind_fit = find(fQ_ust_bin_avg_all{i}>fQ_zQ_ust);
    ust_fit_all{i} = ust_bin_avg_all{i}(ind_fit);
    sigma_ust_fit = ust_bin_SE_all{i}(ind_fit);
    zq_fit = zq_ust_bin_avg_all{i}(ind_fit);
    sigma_zq_fit = zq_ust_bin_SE_all{i}(ind_fit);
    zqnorm_fit = zqnorm_ust_bin_avg_all{i}(ind_fit);
    sigma_zqnorm_fit = zqnorm_ust_bin_SE_all{i}(ind_fit);

    [a, b, sigma_a, sigma_b, zq_pred, ~] = linearfit(ust_fit_all{i}, zq_fit, sigma_ust_fit, sigma_zq_fit);
    intercept_zq_ust_all(i) = a; %intercept of fit
    sigmaintercept_zq_ust_all(i) = sigma_a; %uncertainty in intercept of fit
    slope_zq_ust_all(i) = b; %slope of fit
    sigmaslope_zq_ust_all(i) = sigma_b; %uncertainty in slope of fit
    zq_pred_all{i} = zq_pred;
    
    [~, ~, ~, ~, zqnorm_pred, ~] = linearfit(ust_fit_all{i}, zqnorm_fit, sigma_ust_fit, sigma_zqnorm_fit);
    zqnorm_pred_all{i} = zqnorm_pred;
    [zqnorm_bar, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_fit, sigma_zqnorm_fit);
    zqnorm_bar_all(i) = zqnorm_bar;
    sigma_zqnorm_bar_all(i) = sigma_zqnorm_bar;
    zqnorm_pred_relativerange_all(i) = (zqnorm_pred(end)-zqnorm_pred(1))/zqnorm_bar_all(i);
    
end

%fit values for zq versus ust - literature
intercept_zq_ust_lit = zeros(N_lit,1); %intercept of fit
sigmaintercept_zq_ust_lit = zeros(N_lit,1); %uncertainty in intercept of fit
slope_zq_ust_lit = zeros(N_lit,1); %slope of fit for zq versus ustar
sigmaslope_zq_ust_lit = zeros(N_lit,1); %uncertainty in slope for zq versus ustar
zq_pred_lit = cell(N_lit,1); %predicted zq for lit data

zqnorm_pred_lit{i} = cell(N_lit,1); %predicted zqnorm
zqnorm_bar_lit = zeros(N_lit,1); %mean of zqnorm
sigma_zqnorm_bar_lit = zeros(N_lit,1); %uncertainty of zqnorm bar
zqnorm_pred_relativerange_lit = zeros(N_lit,1); %relative range of zqnorm prediction

for i = 1:N_lit
    [a, b, sigma_a, sigma_b, zq_pred, ~] = linearfit(ust_lit{i},zq_lit{i},sigma_ust_lit{i},sigma_zq_lit{i});
    intercept_zq_ust_lit(i) = a; %intercept of fit
    sigmaintercept_zq_ust_lit(i) = sigma_a; %uncertainty in intercept of fit
    slope_zq_ust_lit(i) = b; %slope of fit
    sigmaslope_zq_ust_lit(i) = sigma_b; %uncertainty in slope of fit
    zq_pred_lit{i} = zq_pred; %predicted zq for lit data
    
    [a, b, sigma_a, sigma_b, zqnorm_pred, ~] = linearfit(ust_lit{i},zqnorm_lit{i},sigma_ust_lit{i},sigma_zqnorm_lit{i});
    zqnorm_pred_lit{i} = zqnorm_pred;
    [zqnorm_bar, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_lit{i}, sigma_zqnorm_lit{i});
    zqnorm_bar_lit(i) = zqnorm_bar;
    sigma_zqnorm_bar_lit(i) = sigma_zqnorm_bar;
    zqnorm_pred_relativerange_lit(i) = (zqnorm_pred(ust_lit{i}==max(ust_lit{i}))-zqnorm_pred(ust_lit{i}==min(ust_lit{i})))/zqnorm_bar_lit(i);
end


%% PLOT ZQ VERSUS UST - PANEL A
figure; subplot(2,1,1); hold on;

%plot binned field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg_all{i}>fQ_zQ_ust);
    errorbar(ust_bin_avg_all{i}(ind_plot),zq_ust_bin_avg_all{i}(ind_plot),zq_ust_bin_SE_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot literature data
for i = 1:N_lit
    errorbar(ust_lit{i},zq_lit{i},sigma_zq_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end

%plot fit for field data
for i = 1:N_Sites;
    plot(ust_fit_all{i},zq_pred_all{i},'Color',Colors_field{i},'LineWidth',LineWidth_field);
end

%plot fit for lit data
for i = 1:N_lit
    plot(ust_lit{i},zq_pred_lit{i},'Color',Colors_lit{i},'LineWidth',LineWidth_lit);
end

%organize plot
xlim([0.25 0.6]);
ylim([0 0.15]);
text(0.255, 0.14,'(a)','FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
ylabel('\textbf{Saltation height, $$z_q$$ (m)}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);


%% PLOT ZQNORM VERSUS UST - PANEL B
subplot(2,1,2); hold on; %alternative subplots if legend is inside

%plot binned field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg_all{i}>fQ_zQ_ust);
    errorbar(ust_bin_avg_all{i}(ind_plot),zqnorm_ust_bin_avg_all{i}(ind_plot),zqnorm_ust_bin_SE_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot literature data
for i = 1:N_lit        
    errorbar(ust_lit{i},zqnorm_lit{i},sigma_zqnorm_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end

%plot fit for field data
for i = 1:N_Sites
    plot(ust_fit_all{i},zqnorm_pred_all{i},'Color',Colors_field{i},'LineWidth',LineWidth_field); 
end
    
for i = 1:2; %only plot first two entries, ignorning Farrell (2012) with no d50
    plot(ust_lit{i},zqnorm_pred_lit{i},'Color',Colors_lit{i},'LineWidth',LineWidth_lit);
end

%organize plot
xlim([0.25 0.6]);
ylim([0 250]);
text(0.255, 235,'(b)','FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
ylabel('\textbf{Dimensionless salt. ht., $$z_q/d_{50}$$}','Interpreter','Latex');
legend_items = [Sites;LitNames];
legend(legend_items,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
print([folder_Plots,'zq_ust.png'],'-dpng');


%% PLOT RAW ZQ / ZQNORM DATA
figure;
for i = 1:N_Sites
    ind_plot = find(fQ_all{i}>fQ_zQ_ust); %get all possible points to plot based on flux frequency

    %zq plot
    subplot(2,N_Sites,i);
    errorbar(ustRe_all{i}(ind_plot),zq_all{i}(ind_plot),sigma_zq_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    xlim([0.25 0.55]);
    ylim([0 0.14]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Saltation height, $$z_q$$ (m)}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);

    %zqnorm plot
    subplot(2,N_Sites,i+N_Sites);
    errorbar(ustRe_all{i}(ind_plot),zqnorm_all{i}(ind_plot),sigma_zqnorm_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    xlim([0.25 0.55]);
    ylim([0 250]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Dimensionless salt. ht., $$z_q/d_{50}$$}','Interpreter','Latex');
    end
    %title(Sites{i})
    set(gca,'FontSize',PlotFont);
end

%print plot for draft
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
print([folder_Plots,'zq_ust_raw.png'],'-dpng');


%% FOR EACH SITE, PERFORM BINNING BY TAU
%initialize lists of tau bins by site
tau_bin_values_all = cell(N_Sites,1); %tau's in each bin
tau_bin_N_all = cell(N_Sites,1); %number of tau in bin
tau_bin_avg_all = cell(N_Sites,1); %avg tau in bin
tau_bin_SE_all = cell(N_Sites,1); %SE tau in bin

%separate fQs into tau bins
fQ_tau_bin_values_all = cell(N_Sites,1); %fQs into tau bins
fQ_tau_bin_avg_all = cell(N_Sites,1); %get average fQ for tau bins
fQ_tau_bin_std_all = cell(N_Sites,1); %get std fQ for tau bins
fQ_tau_bin_SE_all = cell(N_Sites,1); %get SE fQ for tau bins

%separate Q's into tau bins
Q_tau_bin_values_all = cell(N_Sites,1); %Q into tau bins
Q_tau_bin_avg_all = cell(N_Sites,1); %get average Q for tau bins
Q_tau_bin_std_all = cell(N_Sites,1); %get std Q for tau bins
Q_tau_bin_SE_all = cell(N_Sites,1); %get SE Q for tau bins

%get binned values for all sites
for i = 1:N_Sites
    
    %use only values above minimum tau and Q
    ind_analysis = find((tauRe_all{i}>tau_bin_lower)&(Q_all{i}>0));
    tau = tauRe_all{i}(ind_analysis);
    fQ = fQ_all{i}(ind_analysis);
    Q = Q_all{i}(ind_analysis);
    
    %get tau bins
    [tau_bin_values, tau_bin_N, tau_bin_min, tau_bin_max, tau_bin_avg, tau_bin_SE] = ...
        Binning(tau, tau_bin_minrange, tau_bin_maxrange, bin_N_min);
    tau_bin_values_all{i} = tau_bin_values; %tau's in each bin
    tau_bin_N_all{i} = tau_bin_N; %number of tau in bin
    tau_bin_avg_all{i} = tau_bin_avg; %avg tau in bin
    tau_bin_SE_all{i} = tau_bin_SE; %SE tau in bin
    
    %get values in bins
    N_tau_bins = length(tau_bin_values); %number of bins
    fQ_tau_bin_values_all{i} = cell(N_tau_bins,1); %initialize fQ's
    Q_tau_bin_values_all{i} = cell(N_tau_bins,1); %initialize Q's
    for j = 1:N_tau_bins
        %get indices of values in bin
        bin_ind = find(tau>=tau_bin_min(j)&tau<=tau_bin_max(j));
        
        %get values in bins
        fQ_tau_bin_values_all{i}{j} = fQ(bin_ind); %fQ's
        Q_tau_bin_values_all{i}{j} = Q(bin_ind); %Q's
    end
    
    %get fQ avg, std, and SE
    fQ_tau_bin_avg_all{i} = cellfun(@mean,fQ_tau_bin_values_all{i});
    fQ_tau_bin_std_all{i} = cellfun(@std,fQ_tau_bin_values_all{i});
    fQ_tau_bin_SE_all{i} = fQ_tau_bin_std_all{i}./cellfun(@length,fQ_tau_bin_values_all{i});
    
    %get Q avg, std, and SE
    Q_tau_bin_avg_all{i} = cellfun(@mean,Q_tau_bin_values_all{i});
    Q_tau_bin_std_all{i} = cellfun(@std,Q_tau_bin_values_all{i});
    Q_tau_bin_SE_all{i} = Q_tau_bin_std_all{i}./cellfun(@length,Q_tau_bin_values_all{i});
end

%recompute standard errors for bins with fewer than bin_N_min values based on nearest Oceano tau
i_Oceano = find(strcmp(Sites,'Oceano'));
ind_SE_Oceano = find(tau_bin_N_all{i_Oceano}>=bin_N_min); %indices of Oceano bins with sufficient values for SE
tau_avg_Oceano = tau_bin_avg_all{i_Oceano}(ind_SE_Oceano); %tau values for these bins
fQ_CV_Oceano = fQ_tau_bin_std_all{i_Oceano}(ind_SE_Oceano)./fQ_tau_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of fQ values for these bins
Q_CV_Oceano = Q_tau_bin_std_all{i_Oceano}(ind_SE_Oceano)./Q_tau_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of Q values for these bins
for i = 1:N_Sites
    ind_recompute_SE = find(tau_bin_N_all{i}<bin_N_min);
    N_ind = length(ind_recompute_SE);
    for j = 1:N_ind
        k = ind_recompute_SE(j);
        tau_bin = tau_bin_avg_all{i}(k); %get avg tau for bin
        tau_diff_Oceano = abs(tau_bin-tau_avg_Oceano); %get distance of this tau from tau's for Oceano
        ind_Oceano = find(tau_diff_Oceano==min(tau_diff_Oceano),1); %get closest tau for Oceano
        fQ_std = fQ_CV_Oceano(ind_Oceano)*fQ_tau_bin_avg_all{i}(k); %std dev fQ is mean fQ value from site times CV for closest tau at Oceano
        fQ_tau_bin_SE_all{i}(k) = fQ_std/sqrt(tau_bin_N_all{i}(k)); %use this std dev to get SE
        Q_std = Q_CV_Oceano(ind_Oceano)*Q_tau_bin_avg_all{i}(k); %std dev Q is mean Q value from site times CV for closest tau at Oceano
        Q_tau_bin_SE_all{i}(k) = Q_std/sqrt(tau_bin_N_all{i}(k)); %use this std dev to get SE
    end
end


%% FITTING Q VS TAU
%compute fit values for Q versus tau - field
tau_fit_all = cell(N_Sites,1); %tau's for fit
slope_Q_tau_all = zeros(N_Sites,1); %slope of fit for Q versus tau
sigmaslope_Q_tau_all = zeros(N_Sites,1); %uncertainty in slope for Q versus tau
tauit_fit_all = zeros(N_Sites,1); %impact threshold stress
sigma_tauit_fit_all = zeros(N_Sites,1); %uncertainty in impact threshold stress
ustit_fit_all = zeros(N_Sites,1); %impact threshold shear velocity
sigma_ustit_fit_all = zeros(N_Sites,1); %uncertainty in impact threshold shear velocity
Q_pred_all = cell(N_Sites,1); %predicted Q
chi2_Qtau_all = zeros(N_Sites,1); %chi2 value for linear fit
p_Qtau_all = zeros(N_Sites,1); %p value for linear fit
% n_nonlinearfit_all = zeros(N_Sites,1); %u* scaling exponent for nonlinear fit
% tauit_nonlinearfit_all = zeros(N_Sites,1); %tauit for nonlinear fit
% C_nonlinearfit_all = zeros(N_Sites,1); %scaling coefficient for nonlinear fit
% tau_pred_nonlinear_all = cell(N_Sites,1); %tau for Q predictions for nonlinear fit
% Q_pred_nonlinear_all = cell(N_Sites,1); %predicted Q from nonlinear fit

%go through sites
for i = 1:N_Sites
    %get values for fitting
    ind_fit = find(fQ_tau_bin_avg_all{i}>fQ_Q_tau);
    tau_fit = tau_bin_avg_all{i}(ind_fit);
    sigma_tau_fit = tau_bin_SE_all{i}(ind_fit);
    Q_fit = Q_tau_bin_avg_all{i}(ind_fit);
    sigma_Q_fit = Q_tau_bin_SE_all{i}(ind_fit);
    
    %perform fitting
    [a, b, sigma_a, sigma_b, Q_pred, ~] = linearfit(tau_fit, Q_fit, sigma_tau_fit, sigma_Q_fit);
    slope_Q_tau_all(i) = b; %slope of fit
    sigmaslope_Q_tau_all(i) = sigma_b; %uncertainty in slope of fit
    tauit_fit_all(i) = -a/slope_Q_tau_all(i);
    sigma_tauit_fit_all(i) = sigma_a/slope_Q_tau_all(i);
    ustit_fit_all(i) = sqrt(tauit_fit_all(i)/rho_a);
    sigma_ustit_fit_all(i) = sigma_tauit_fit_all(i)*(1/(2*rho_a))*(1/ustit_fit_all(i));

    %compute chi2 and p-value
    sigma_Q_chi2 = sqrt(sigma_Q_fit.^2 + b^2*sigma_tau_fit.^2); %get sigma_Q values that combine uncertainties in Q and tau
    [chi2, p] = Chi2Calculation(Q_fit, sigma_Q_chi2, Q_pred, length(Q_pred-2));
    chi2_Qtau_all(i) = chi2;
    p_Qtau_all(i) = p;
    
%     %compute nonlinear fit values
%     [n, tauit, C] = NonlinearFluxFit(tau_fit, Q_fit, sigma_tau_fit, sigma_Q_fit, 100);
%     n_nonlinearfit_all(i) = n; %u* scaling exponent for nonlinear fit
%     tauit_nonlinearfit_all(i) = tauit; %tauit for nonlinear fit
%     C_nonlinearfit_all(i) = C; %scaling coefficient for nonlinear fit
%     
%     %get predictions of Q for nonlinear flux fit
%     tau_pred_nonlinear_all{i} = linspace(tauit,max(tau_fit),50);
%     ust_pred_nonlinear = sqrt(tau_pred_nonlinear_all{i}/rho_a);
%     Q_pred_nonlinear_all{i} = C*ust_pred_nonlinear.^n.*(tau_pred_nonlinear_all{i}-tauit);
    
    %add 0 point to prediction / fit values for plotting
    Q_pred_all{i} = [0; Q_pred]; 
    tau_fit_all{i} = [tauit_fit_all(i); tau_fit];
end


%% PLOT Q VERSUS TAU
figure; hold on;
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    errorbar(tau_bin_avg_all{i},Q_tau_bin_avg_all{i},Q_tau_bin_SE_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i});
    plot(tau_fit_all{i},Q_pred_all{i},'Color',Colors_field{i}); %plot fit for these values
    legend_items{2*i-1} = Sites{i}; %add to list of legend items
    legend_items{2*i} = 'linear fit'; %add to list of legend items
    %legend_items{3*i} = ['nonlinear, n=',num2str(n_nonlinearfit_all(i),'%.2f')]; %add to list of legend items
    %plot(tau_pred_nonlinear_all{i},Q_pred_nonlinear_all{i},'--','Color',Colors_field{i}); %plot nonlinear fit
end
xlim([0 0.35]);
ylim([0 50]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{$$Q$$ (g/m/s)}','Interpreter','Latex');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',PlotFont);
set(gcf, 'PaperPosition',[0 0 8 6]);
print([folder_Plots,'Q_tau.png'],'-dpng');


%% PLOT Q VERSUS TAU RAW DATA
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i);
    errorbar(tauRe_all{i},Q_all{i},sigma_Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{i});
    xlim([0 0.35]);
    ylim([0 50]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$Q$$ (g/m/s)}','Interpreter','Latex');
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end
set(gcf, 'PaperPosition',[0 0 11 4]);
print([folder_Plots,'Q_tau_raw.png'],'-dpng');


%% FOR EACH SITE, PERFORM BINNING BY TAUEX

%compute tauex based on fitted thresholds
tauex_all = cell(N_Sites,1);
for i=1:N_Sites
    tauex_all{i} = tauRe_all{i}-tauit_fit_all(i);
end

%compute Qhat
Qhat_all = cell(N_Sites,1);
sigma_Qhat_all = cell(N_Sites,1);
for i=1:N_Sites
    Qhat_all{i} = Q_all{i}./((1/g)*ustit_fit_all(i).*1e3*tauex_all{i});
    sigma_Qhat_all{i} = sigma_Q_all{i}./((1/g)*ustit_fit_all(i).*1e3*tauex_all{i});
end

%initialize lists of tauex bins by site
tauex_bin_values_all = cell(N_Sites,1); %tauex's in each bin
tauex_bin_N_all = cell(N_Sites,1); %number of tauex in bin
tauex_bin_avg_all = cell(N_Sites,1); %avg tauex in bin
tauex_bin_SE_all = cell(N_Sites,1); %SE tauex in bin
tauex_bin_min_all = cell(N_Sites,1); %minimum tauex in bin

%separate fQs into tauex bins
fQ_tauex_bin_values_all = cell(N_Sites,1); %fQs into tauex bins
fQ_tauex_bin_avg_all = cell(N_Sites,1); %get average fQ for tauex bins
fQ_tauex_bin_std_all = cell(N_Sites,1); %get std fQ for tauex bins
fQ_tauex_bin_SE_all = cell(N_Sites,1); %get SE fQ for tauex bins

%separate Q's into tauex bins
Q_tauex_bin_values_all = cell(N_Sites,1); %Q into tauex bins
Q_tauex_bin_avg_all = cell(N_Sites,1); %get average Q for tauex bins
Q_tauex_bin_std_all = cell(N_Sites,1); %get std Q for tauex bins
Q_tauex_bin_SE_all = cell(N_Sites,1); %get SE Q for tauex bins

%get binned values for all sites
for i = 1:N_Sites
    
    %use only values above minimum tauex and Q
    ind_analysis = find(tauex_all{i}>tauex_bin_lower);
    tauex = tauex_all{i}(ind_analysis);
    fQ = fQ_all{i}(ind_analysis);
    Q = Q_all{i}(ind_analysis);
    
    %get tauex bins
    [tauex_bin_values, tauex_bin_N, tauex_bin_min, tauex_bin_max, tauex_bin_avg, tauex_bin_SE] = ...
        Binning(tauex, tauex_bin_minrange, tauex_bin_maxrange, bin_N_min);
    tauex_bin_values_all{i} = tauex_bin_values; %tauex's in each bin
    tauex_bin_N_all{i} = tauex_bin_N; %number of tauex in bin
    tauex_bin_avg_all{i} = tauex_bin_avg; %avg tauex in bin
    tauex_bin_SE_all{i} = tauex_bin_SE; %SE tauex in bin
    tauex_bin_min_all{i} = tauex_bin_min; %min tauex in bin
    
    %get values in bins
    N_tauex_bins = length(tauex_bin_values); %number of bins
    fQ_tauex_bin_values_all{i} = cell(N_tauex_bins,1); %initialize fQ's
    Q_tauex_bin_values_all{i} = cell(N_tauex_bins,1); %initialize Q's
    for j = 1:N_tauex_bins
        %get indices of values in bin
        bin_ind = find(tauex>=tauex_bin_min(j)&tauex<=tauex_bin_max(j));
        
        %get values in bins
        fQ_tauex_bin_values_all{i}{j} = fQ(bin_ind); %fQ's
        Q_tauex_bin_values_all{i}{j} = Q(bin_ind); %Q's
    end
    
    %get fQ avg, std, and SE
    fQ_tauex_bin_avg_all{i} = cellfun(@mean,fQ_tauex_bin_values_all{i});
    fQ_tauex_bin_std_all{i} = cellfun(@std,fQ_tauex_bin_values_all{i});
    fQ_tauex_bin_SE_all{i} = fQ_tauex_bin_std_all{i}./cellfun(@length,fQ_tauex_bin_values_all{i});
    
    %get Q avg, std, and SE
    Q_tauex_bin_avg_all{i} = cellfun(@mean,Q_tauex_bin_values_all{i});
    Q_tauex_bin_std_all{i} = cellfun(@std,Q_tauex_bin_values_all{i});
    Q_tauex_bin_SE_all{i} = Q_tauex_bin_std_all{i}./cellfun(@length,Q_tauex_bin_values_all{i});
end

%recompute standard errors for bins with fewer than bin_N_min values based on nearest Oceano tauex
i_Oceano = find(strcmp(Sites,'Oceano'));
ind_SE_Oceano = find(tauex_bin_N_all{i_Oceano}>=bin_N_min); %indices of Oceano bins with sufficient values for SE
tauex_avg_Oceano = tauex_bin_avg_all{i_Oceano}(ind_SE_Oceano); %tauex values for these bins
fQ_CV_Oceano = fQ_tauex_bin_std_all{i_Oceano}(ind_SE_Oceano)./fQ_tauex_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of fQ values for these bins
Q_CV_Oceano = Q_tauex_bin_std_all{i_Oceano}(ind_SE_Oceano)./Q_tauex_bin_avg_all{i_Oceano}(ind_SE_Oceano); %std/mean of Q values for these bins
for i = 1:N_Sites
    ind_recompute_SE = find(tauex_bin_N_all{i}<bin_N_min);
    N_ind = length(ind_recompute_SE);
    for j = 1:N_ind
        k = ind_recompute_SE(j);
        tauex_bin = tauex_bin_avg_all{i}(k); %get avg tauex for bin
        tauex_diff_Oceano = abs(tauex_bin-tauex_avg_Oceano); %get distance of this tauex from tauex's for Oceano
        ind_Oceano = find(tauex_diff_Oceano==min(tauex_diff_Oceano),1); %get closest tauex for Oceano
        fQ_std = fQ_CV_Oceano(ind_Oceano)*fQ_tauex_bin_avg_all{i}(k); %std dev fQ is mean fQ value from site times CV for closest tauex at Oceano
        fQ_tauex_bin_SE_all{i}(k) = fQ_std/sqrt(tauex_bin_N_all{i}(k)); %use this std dev to get SE
        Q_std = Q_CV_Oceano(ind_Oceano)*Q_tauex_bin_avg_all{i}(k); %std dev Q is mean Q value from site times CV for closest tauex at Oceano
        Q_tauex_bin_SE_all{i}(k) = Q_std/sqrt(tauex_bin_N_all{i}(k)); %use this std dev to get SE
    end
end

%compute Qhat_bin from bin averages
Qhat_tauex_bin_avg_all = cell(N_Sites,1); %get average Qhat for tauex bins
Qhat_tauex_bin_SE_all = cell(N_Sites,1); %get SE Qhat for tauex bins

for i=1:N_Sites
    Q = Q_tauex_bin_avg_all{i}*1e-3;
    sigma_Q = Q_tauex_bin_SE_all{i}*1e-3;
    ustit = ustit_fit_all(i);
    sigma_ustit = sigma_ustit_fit_all(i);
    tauex = tauex_bin_avg_all{i};
    sigma_tauex = tauex_bin_SE_all{i};
    
    Qhat_tauex_bin_avg_all{i} = (g/ustit)*(Q./tauex); %compute average Qhat
    sigma_Qhat_Q = (g/ustit)*(sigma_Q./tauex); %compute contribution of Q to sigma_Qhat
    sigma_Qhat_ustit = (sigma_ustit/ustit)*(g/ustit)*(Q./tauex); %compute average Qhat
    sigma_Qhat_tauex = (sigma_tauex./tauex).*(g/ustit).*(Q./tauex); %compute average Qhat
    Qhat_tauex_bin_SE_all{i} = sqrt(sigma_Qhat_ustit.^2+sigma_Qhat_tauex.^2+sigma_Qhat_Q.^2); %compute total uncertainty for Qhat
end


%% FITTING Q VS TAUEX
%compute fit values for Q versus tauex - field
tauex_fit_all = cell(N_Sites,1); %tauex's for fit
Qtauratio_bar_all = zeros(N_Sites,1); %scaling coefficient for Q/tauex versus tauex 
Qhat_bar_all = zeros(N_Sites,1); %scaling coefficient for Qhat versus tauex
sigma_Qhat_bar_all = zeros(N_Sites,1); %scaling coefficient for Qhat versus tauex

%go through sites
for i = 1:N_Sites
    %get values for fitting
    ind_fit = intersect(find(fQ_tauex_bin_avg_all{i}>fQ_Q_tauex),... %must have sufficient transport
        find(tauex_bin_min_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i))); %must exceed threshold
    tauex_fit_all{i} = tauex_bin_avg_all{i}(ind_fit);
    sigma_tauex_fit = tauex_bin_SE_all{i}(ind_fit);
    Q_fit = Q_tauex_bin_avg_all{i}(ind_fit);
    sigma_Q_fit = Q_tauex_bin_SE_all{i}(ind_fit);
    Qtauratio_fit = Q_fit./tauex_fit_all{i};
    sigma_Qtauratio_fit = sigma_Q_fit./tauex_fit_all{i};
    Qhat_fit = Qhat_tauex_bin_avg_all{i}(ind_fit);
    sigma_Qhat_fit = Qhat_tauex_bin_SE_all{i}(ind_fit);
    
    %perform fitting
    [Qtauratio_bar, ~] = MeanUncertainty(Qtauratio_fit, sigma_Qtauratio_fit);
    Qtauratio_bar_all(i) = Qtauratio_bar;
    [Qhat_bar, sigma_Qhat_bar] = MeanUncertainty(Qhat_fit, sigma_Qhat_fit);
    Qhat_bar_all(i) = Qhat_bar;
    sigma_Qhat_bar_all(i) = sigma_Qhat_bar;
end


%% PLOT Q VERSUS TAUEX BINNED DATA
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    ind_plot = intersect(find(fQ_tauex_bin_avg_all{i}>fQ_Q_tauex),... %must have sufficient transport
        find(tauex_bin_min_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i))); %must exceed threshold
    errorbar(tauex_bin_avg_all{i}(ind_plot),Q_tauex_bin_avg_all{i}(ind_plot),Q_tauex_bin_SE_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit, Qtauratio_bar_all(i)*tauex_fit,'Color',Colors_field{i});
end

xlim([0 0.3]);
ylim([0 50]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%draft paper plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'Q_tauex.png'],'-dpng');


%% PLOT Q VERSUS TAUEX RAW DATA
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i);
    ind_plot = intersect(find(fQ_all{i}>fQ_Q_tauex),... %must have sufficient transport
        find(tauex_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i))); %must exceed threshold
    errorbar(tauex_all{i}(ind_plot),Q_all{i}(ind_plot),sigma_Q_all{i}(ind_plot),Markers_field{i},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{i});
    xlim([0 0.3]);
    ylim([0 50]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end
set(gcf, 'PaperPosition',[0 0 11 4]);
print([folder_Plots,'Q_tauex_raw.png'],'-dpng');


%% PLOT QHAT VS TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    ind_plot = intersect(find(fQ_tauex_bin_avg_all{i}>fQ_Q_tauex),... %must have sufficient transport
        find(tauex_bin_min_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i))); %must exceed threshold
    errorbar(tauex_bin_avg_all{i}(ind_plot),Qhat_tauex_bin_avg_all{i}(ind_plot),Qhat_tauex_bin_SE_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit,Qhat_bar_all(i)*ones(2,1),'Color',Colors_field{i});
end
xlim([0 0.3]);
ylim([0 9]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Dimensionless saltation flux, $$\hat{Q}$$}','Interpreter','Latex');
legend(Sites,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%draft plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'Qhat_tauex.png'],'-dpng');


%% PLOT QHAT VERSUS TAUEX RAW DATA
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i);
    ind_plot = intersect(find(fQ_all{i}>fQ_Q_tauex),... %must have sufficient transport
        find(tauex_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i))); %must exceed threshold
    errorbar(tauex_all{i}(ind_plot),Qhat_all{i}(ind_plot),sigma_Qhat_all{i}(ind_plot),Markers_field{i},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{i});
    xlim([0 0.3]);
    ylim([0 9]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{Dimensionless saltation flux, $$\hat{Q}$$}','Interpreter','Latex');
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end
set(gcf, 'PaperPosition',[0 0 11 4]);
print([folder_Plots,'Qhat_tauex_raw.png'],'-dpng');

% 
% %% FOR EACH SITE, PERFORM BINNING BY UST*TAUEX
% %compute ust*tauex
% usttauex_all = cell(N_Sites,1);
% for i=1:N_Sites
%     usttauex_all{i} = ustRe_all{i}.*tauex_all{i};
% end
% 
% %create initial tauex bins
% usttauex_bin_lower = 0;
% usttauex_bin_upper = 0.145;
% usttauex_bin_increment = 0.005;
% usttauex_bin_edges_init = usttauex_bin_lower:usttauex_bin_increment:usttauex_bin_upper;
% usttauex_bin_min_init = usttauex_bin_edges_init(1:end-1);
% usttauex_bin_max_init = usttauex_bin_edges_init(2:end);
% N_usttauex_bins_init = length(usttauex_bin_edges_init)-1;
% 
% %initialize lists of usttauex bins by site
% usttauex_bin_values = cell(N_Sites,1); %usttauex's in each bin
% N_usttauex_bin_values = cell(N_Sites,1); %N usttauex's in each bin
% usttauex_bin_avg = cell(N_Sites,1); %avg usttauex in bin
% usttauex_bin_std = cell(N_Sites,1); %std usttauex in bin
% usttauex_bin_SE = cell(N_Sites,1); %SE usttauex in bin
% usttauex_bin_min = cell(N_Sites,1); %minimum usttauex in bin
% 
% %separate Qs into usttauex bins
% Q_usttauex_bin_values = cell(N_Sites,1); %Qs into usttauex bins
% Q_usttauex_bin_avg = cell(N_Sites,1); %get average Q for usttauex bins
% Q_usttauex_bin_std = cell(N_Sites,1); %get std Q for usttauex bins
% Q_usttauex_bin_SE = cell(N_Sites,1); %get SE Q for usttauex bins
% 
% %go through all sites
% for i = 1:N_Sites
%     
%     %get usttauex bins
%     bin_N_init = histcounts(usttauex_all{i},usttauex_bin_edges_init); %get number of usttauex in each initial bin
%     usttauex_bin_min_Site = []; %initialize list of usttauex bin mins for site
%     usttauex_bin_max_Site = []; %initialize list of usttauex bin maxes for site
%     cum_N = 0; %initialize cumulative usttauex counts for site bin
%     usttauex_bin_min_this = usttauex_bin_min_init(1); %initialize usttauex bin min for this bin
%     for j = 1:N_usttauex_bins_init; %go through all initial usttauex bins to determine usttauex bins for site
%         cum_N = cum_N+bin_N_init(j); %add to cumulative usttauex counts
%         if cum_N>=bin_N_min; %if cumulative usttauex counts exceeds limit
%             usttauex_bin_max_this = usttauex_bin_max_init(j); %set usttauex bin max for this bin
%             usttauex_bin_min_Site = [usttauex_bin_min_Site, usttauex_bin_min_this]; %add to list of usttauex bin mins for Site
%             usttauex_bin_max_Site = [usttauex_bin_max_Site, usttauex_bin_max_this]; %add to list of usttauex bin maxes for Site
%             cum_N = 0; %reset cumulative usttauex counts for site bin
%             usttauex_bin_min_this = usttauex_bin_max_init(j); %reset usttauex bin min for next bin
%         end
%     end
%     N_usttauex_bins_Site = length(usttauex_bin_min_Site); %get number of usttauex bins for Site
%    
%     %initialize usttauex bins
%     usttauex_bin_values{i} = cell(N_usttauex_bins_Site,1);
%     N_usttauex_bin_values{i} = zeros(N_usttauex_bins_Site,1)*NaN;
%     usttauex_bin_avg{i} = zeros(N_usttauex_bins_Site,1)*NaN;
%     usttauex_bin_std{i} = zeros(N_usttauex_bins_Site,1)*NaN;
%     usttauex_bin_SE{i} = zeros(N_usttauex_bins_Site,1)*NaN;
%     usttauex_bin_min{i} = zeros(N_usttauex_bins_Site,1)*NaN;
% 
%     %initialize Q's into usttauex bins
%     Q_usttauex_bin_values{i} = cell(N_usttauex_bins_Site,1);
%     Q_usttauex_bin_avg{i} = zeros(N_usttauex_bins_Site,1)*NaN;
%     Q_usttauex_bin_std{i} = zeros(N_usttauex_bins_Site,1)*NaN;
%     Q_usttauex_bin_SE{i} = zeros(N_usttauex_bins_Site,1)*NaN;
% 
%     %get values, avgs, and standard deviations for usttauex bins
%     for j = 1:N_usttauex_bins_Site
%         %get indices
%         bin_ind = find(usttauex_all{i}>=usttauex_bin_min_Site(j)&usttauex_all{i}<=usttauex_bin_max_Site(j));
% 
%         %get usttauex in bin
%         usttauex_bin_values{i}{j} = usttauex_all{i}(bin_ind);
%         N_usttauex_bin_values{i}(j) = length(usttauex_bin_values{i}{j});
%         usttauex_bin_avg{i}(j) = mean(usttauex_bin_values{i}{j});
%         usttauex_bin_std{i}(j) = std(usttauex_bin_values{i}{j});
%         usttauex_bin_SE{i}(j) = usttauex_bin_std{i}(j)./sqrt(N_usttauex_bin_values{i}(j));
%         usttauex_bin_min{i}(j) = min(usttauex_bin_values{i}{j});
%        
%         %get Qs for usttauex bins
%         Q_usttauex_bin_values{i}{j} = Q_all{i}(bin_ind);
%         Q_usttauex_bin_avg{i}(j) = mean(Q_usttauex_bin_values{i}{j});
%         Q_usttauex_bin_std{i}(j) = std(Q_usttauex_bin_values{i}{j});
%         Q_usttauex_bin_SE{i}(j) = Q_usttauex_bin_std{i}(j)./sqrt(N_usttauex_bin_values{i}(j));
%     end
% end
% 
% 
% %% FITS TO Q UST*TAUEX
% 
% %initialize scaling coefficient for Q versus ust*tauex
% C_Qusttauex_fit = zeros(N_Sites,1);
% Qpred_usttauex = cell(N_Sites,1);
% chi2_Qusttauex = zeros(N_Sites,1);
% P_Qusttauex = zeros(N_Sites,1);
% 
% %perform fits
% for i = 1:N_Sites
%    
%     %compute scaling coefficient for Q versus ust*tauex
%     C_Qusttauex_fit(i) = mean(Q_usttauex_bin_avg{i}./usttauex_bin_avg{i});
%     
%     %predict values of Q versus ust*tauex
%     Qpred_usttauex{i} = usttauex_bin_avg{i}*C_Qusttauex_fit(i);
%     
%     %compute chi2 value based on prediction
%     chi2_Qusttauex(i) = sum((1./Q_usttauex_bin_SE{i}.^2).*(Q_usttauex_bin_avg{i}-Qpred_usttauex{i}).^2);
%     N_DF = length(usttauex_bin_avg{i})-2; %number of degrees of freedom
%     
%     P_Qusttauex(i) = chi2cdf(chi2_Qusttauex(i),N_DF);
% end
% 
% 
% %% PLOT Q VERSUS UST*TAUEX
% figure; clf; hold on;
% 
% %plot binned values
% for i = 1:N_Sites
%     errorbar(usttauex_bin_avg{i},Q_usttauex_bin_avg{i},Q_usttauex_bin_SE{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
% end
% 
% %plot fit values
% usttauex_fit = [0 0.14];
% for i = 1:N_Sites
%     plot(usttauex_bin_avg{i}, Qpred_usttauex{i},'Color',Colors_field{i});
% end
% 
% xlim([0 0.14]);
% ylim([0 45]);
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% xlabel('\textbf{Shear velocity - excess stress product, $$u_{*} \tau_{ex}$$ (kg s$$^{-3}$$)}','Interpreter','Latex');
% ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
% legend(Sites,'Location','SouthEast');
% set(gca,'FontSize',PlotFont);
% 
% %draft paper plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
% print([folder_Plots,'Q_usttauex.png'],'-dpng');
% 
% 
% 
%% save analysis data
save([folder_AnalysisData,'StressFluxWindows_Analysis.mat'],'*all','*lit');
% 
% 
% % %% DIFFERENT METRICS OF WIND
% % i = 3;
% % 
% % %mean wind velocity
% % figure; clf;
% % subplot(3,1,1); hold on;
% % plot(tauRe_all{i},Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/4);
% % xlabel('Shear stress, $$\tau_{Re}$$ (Pa)','Interpreter','Latex');
% % ylabel('Saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
% % set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% % set(gca,'FontSize',PlotFont);
% % 
% % subplot(3,1,2); hold on;
% % plot(ubar_all{i},Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/4);
% % xlabel('Mean wind velocity, $$\overline{u}$$ (m/s)','Interpreter','Latex');
% % ylabel('Saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
% % set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% % set(gca,'FontSize',PlotFont);
% % 
% % subplot(3,1,3); hold on;
% % plot(abs(uw_all{i}),Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/4);
% % xlabel('Momentum flux, $$\overline{uw^{\prime}} \textrm{ (m}^2/\textrm{s}^2)$$ ','Interpreter','Latex');
% % ylabel('Saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
% % set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% % set(gca,'FontSize',PlotFont);
% % 
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
% % print([folder_Plots,'Q_wind_raw.png'],'-dpng');
% % 
% %
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
% %
% % %% PERFORM FITS FOR FQ
% % %set fQ bins
% % fQ_bins_min = 0:0.05:0.95;
% % fQ_bins_max = 0.05:0.05:1;
% % fQ_bins_mid = mean([fQ_bins_min; fQ_bins_max]);
% % N_fQ_bins = length(fQ_bins_mid);
% % fQ_TFEM_fit_min = [0.65,0.7,0.1]; %mininum for fitting to obtain impact/fluid thresholds
% % fQ_TFEM_fit_max = [0.9,0.85,0.9]; %maximum for fitting to obtain impact/fluid thresholds
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