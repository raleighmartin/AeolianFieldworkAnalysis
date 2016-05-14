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

%analysis parameters
fQ_Q_tau = 0.1; %minimum mean fQ for Q versus tau comparison (and determination of threshold)

%binning information
bin_N_min = 3; %minimum number of entries for bin
ust_bin_lower = 0;
ust_bin_minrange = 0.005; %minimum difference between upper and lower value in bin
ust_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
tau_bin_lower = 0;
tau_bin_minrange = 0.005; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
tauex_bin_lower = 0;
tauex_bin_minrange = 0.005; %minimum difference between upper and lower value in bin
tauex_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin

%fitting parameters
fQ_zQ_ust = 0; %mininum fQ for zQ versus u* comparison (for constant flux height)
fQ_Q_tauex = 0; %minimum fQ for Q versus tauex comparison
fQ_Qhat_tauex = 0; %minimum fQ for Qhat versus tauex comparison
N_sigma_tauit = 0; %std deviations from impact threshold for plotting

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
N_Lit = 3;

%change zqnorm_lit so that Farrell doesn't show up 
zqnorm_lit{3} = zqnorm_lit{3}*1e6;

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
ust_bin_ind = cell(N_Sites,1); %indices of ustRe_all for bin
ust_bin_values = cell(N_Sites,1); %ust's in each bin
N_ust_bin_values = cell(N_Sites,1); %N ust's in each bin
ust_bin_avg = cell(N_Sites,1); %avg ust in bin
ust_bin_std = cell(N_Sites,1); %std ust in bin
ust_bin_SE = cell(N_Sites,1); %SE ust in bin

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

%go through all sites
for i = 1:N_Sites
    
    %get u* bins
    ust_analysis = ustRe_all{i}((ustRe_all{i}>ust_bin_lower)&(Q_all{i}>0));
    N_ust = length(ust_analysis); %get number of u*
    ust_sort = sort(ust_analysis); %sort u*
    ind_bin = 1; %initialize index of bin
    ust_bin_list = []; %initialize list of ust in bin
    for j=1:N_ust
        if isempty(ust_bin_list) %if no u*'s in bin, create first value for bin
            ust_bin_list = ust_sort(j); %add u* value to list
        elseif (ust_sort(j)-min(ust_bin_list))<=ust_bin_minrange; %if value of u* is within min range, add to list
            ust_bin_list = [ust_bin_list; ust_sort(j)];     
        elseif (ust_sort(j)-min(ust_bin_list))>ust_bin_maxrange; %if value of u* is outside max range, create new bin
            ust_bin_values{i}{ind_bin} = ust_bin_list; %add current u* list to cell array
            ind_bin = ind_bin+1; %increment to next u* bin
            ust_bin_list = ust_sort(j); %add u* value to new list
        elseif length(ust_bin_list)<bin_N_min; %if next u* is within min-max range, but have not reached minimum number of entries, add to list
            ust_bin_list = [ust_bin_list; ust_sort(j)]; 
        elseif (N_ust-j)<bin_N_min; %if new bin could be created, but there are few entries left at the end of the list, add to current list
            ust_bin_list = [ust_bin_list; ust_sort(j)];
        else %create new bin
            ust_bin_values{i}{ind_bin} = ust_bin_list; %add current u* list to cell array
            ind_bin = ind_bin+1; %increment to next u*bin
            ust_bin_list = ust_sort(j); %add u* value to new list
        end
        if j == N_ust %if at the end of u*s, add current u* list to cell array
            ust_bin_values{i}{ind_bin} = ust_bin_list;
        end
    end


        elseif length(ust_bin_list)<=bin_N_min;
            ust_bin_minrange
            if 
            range([ust_bin_list; ust_sort(ind)])
            
        cum_N = cum_N+1; %add to cumulative u* counts
        %ust_range = range(ustRe_all{i}((ustRe_all{i}>=ust_bin_min_this&ustRe_all{i}<=ust_bin_max_init(j))&(Q_all{i}>0))); 
        if cum_N>=bin_N_min; %if cumulative u* counts exceeds limit
            ust_bin_max_this = ust_bin_max_init(j); %set u* bin max for this bin
            ust_bin_min_Site = [ust_bin_min_Site, ust_bin_min_this]; %add to list of u* bin mins for Site
            ust_bin_max_Site = [ust_bin_max_Site, ust_bin_max_this]; %add to list of u* bin maxes for Site
            cum_N = 0; %reset cumulative u* counts for site bin
            ust_bin_min_this = ust_bin_max_init(j); %reset u* bin min for next bin
        end
    end
    
    ust_sort = ust_sort(ind_ust_sort);
    bin_N_init = histcounts(ustRe_all{i}(Q_all{i}>0),ust_bin_edges_init); %get number of u* in each initial bin, limit only to entries with Q>0
    ust_bin_min_Site = []; %initialize list of u* bin mins for site
    ust_bin_max_Site = []; %initialize list of u* bin maxes for site
    cum_N = 0; %initialize cumulative u* counts for site bin
    ust_bin_min_this = ust_bin_min_init(1); %initialize u* bin min for this bin
    for j = 1:N_ust_bins_init; %go through all initial u* bins to determine u* bins for site
        cum_N = cum_N+bin_N_init(j); %add to cumulative u* counts
        %ust_range = range(ustRe_all{i}((ustRe_all{i}>=ust_bin_min_this&ustRe_all{i}<=ust_bin_max_init(j))&(Q_all{i}>0))); 
        if cum_N>=bin_N_min; %if cumulative u* counts exceeds limit
            ust_bin_max_this = ust_bin_max_init(j); %set u* bin max for this bin
            ust_bin_min_Site = [ust_bin_min_Site, ust_bin_min_this]; %add to list of u* bin mins for Site
            ust_bin_max_Site = [ust_bin_max_Site, ust_bin_max_this]; %add to list of u* bin maxes for Site
            cum_N = 0; %reset cumulative u* counts for site bin
            ust_bin_min_this = ust_bin_max_init(j); %reset u* bin min for next bin
        end
    end
    N_ust_bins_Site = length(ust_bin_min_Site); %get number of u* bins for Site
   
    %initialize u* bins
    ust_bin_values{i} = cell(N_ust_bins_Site,1);
    N_ust_bin_values{i} = zeros(N_ust_bins_Site,1)*NaN;
    ust_bin_avg{i} = zeros(N_ust_bins_Site,1)*NaN;
    ust_bin_std{i} = zeros(N_ust_bins_Site,1)*NaN;
    ust_bin_SE{i} = zeros(N_ust_bins_Site,1)*NaN;

    %initialize fQ's into ust bins
    fQ_ust_bin_values{i} = cell(N_ust_bins_Site,1);
    fQ_ust_bin_avg{i} = zeros(N_ust_bins_Site,1)*NaN;
    fQ_ust_bin_std{i} = zeros(N_ust_bins_Site,1)*NaN;
    fQ_ust_bin_SE{i} = zeros(N_ust_bins_Site,1)*NaN;
    
    %initialize zq's into ust bins
    zq_ust_bin_values{i} = cell(N_ust_bins_Site,1);
    zq_ust_bin_avg{i} = zeros(N_ust_bins_Site,1)*NaN;
    zq_ust_bin_std{i} = zeros(N_ust_bins_Site,1)*NaN;
    zq_ust_bin_SE{i} = zeros(N_ust_bins_Site,1)*NaN;
    
    %initialize zqnorm's into ust bins
    zqnorm_ust_bin_values{i} = cell(N_ust_bins_Site,1);
    zqnorm_ust_bin_avg{i} = zeros(N_ust_bins_Site,1)*NaN;
    zqnorm_ust_bin_std{i} = zeros(N_ust_bins_Site,1)*NaN;
    zqnorm_ust_bin_SE{i} = zeros(N_ust_bins_Site,1)*NaN;
    
    %initialize Q's into ust bins
    Q_ust_bin_values{i} = cell(N_ust_bins_Site,1);
    Q_ust_bin_avg{i} = zeros(N_ust_bins_Site,1)*NaN;
    Q_ust_bin_std{i} = zeros(N_ust_bins_Site,1)*NaN;
    Q_ust_bin_SE{i} = zeros(N_ust_bins_Site,1)*NaN;
    
    %get values, avgs, and standard deviations for ust bins
    for j = 1:N_ust_bins_Site
        %get indices
        %bin_ind = find(ustRe_all{i}>=ust_bin_min_Site(j)&ustRe_all{i}<=ust_bin_max_Site(j));
        bin_ind = find((ustRe_all{i}>=ust_bin_min_Site(j)&ustRe_all{i}<=ust_bin_max_Site(j))&(Q_all{i}>0)); %limit only to entries with Q>0

        %get ust in bin
        ust_bin_values{i}{j} = ustRe_all{i}(bin_ind);
        N_ust_bin_values{i}(j) = length(ust_bin_values{i}{j});
        ust_bin_avg{i}(j) = mean(ust_bin_values{i}{j});
        ust_bin_std{i}(j) = std(ust_bin_values{i}{j});
        ust_bin_SE{i}(j) = ust_bin_std{i}(j)./sqrt(N_ust_bin_values{i}(j));
        
        %get fQs for ust bins
        fQ_ust_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_ust_bin_avg{i}(j) = mean(fQ_ust_bin_values{i}{j});
        fQ_ust_bin_std{i}(j) = std(fQ_ust_bin_values{i}{j});
        fQ_ust_bin_SE{i}(j) = fQ_ust_bin_std{i}(j)./sqrt(N_ust_bin_values{i}(j));
        
        %get zqs for ust bins
        zq_ust_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_ust_bin_avg{i}(j) = mean(zq_ust_bin_values{i}{j});
        zq_ust_bin_std{i}(j) = std(zq_ust_bin_values{i}{j});
        zq_ust_bin_SE{i}(j) = zq_ust_bin_std{i}(j)./sqrt(N_ust_bin_values{i}(j));
        
        %get zqnorms for ust bins
        zqnorm_ust_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_ust_bin_avg{i}(j) = mean(zqnorm_ust_bin_values{i}{j});
        zqnorm_ust_bin_std{i}(j) = std(zqnorm_ust_bin_values{i}{j});
        zqnorm_ust_bin_SE{i}(j) = zqnorm_ust_bin_std{i}(j)./sqrt(N_ust_bin_values{i}(j));
        
        %get Qs for ust bins
        Q_ust_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_ust_bin_avg{i}(j) = mean(Q_ust_bin_values{i}{j});
        Q_ust_bin_std{i}(j) = std(Q_ust_bin_values{i}{j});
        Q_ust_bin_SE{i}(j) = Q_ust_bin_std{i}(j)./sqrt(N_ust_bin_values{i}(j));
    end
end


%% PLOT ZQ VERSUS UST - PANEL A
figure; subplot(2,1,1); hold on; %alternative subplots if legend is inside

%plot binned field data
for i = 1:N_Sites
    ind_plot = find(fQ_ust_bin_avg{i}>fQ_zQ_ust);
    errorbar(ust_bin_avg{i}(ind_plot),zq_ust_bin_avg{i}(ind_plot),zq_ust_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot literature data
for i = 1:N_Lit
    errorbar(ust_lit{i},zq_lit{i},sigma_zq_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end

%fit values for zq versus ust - field
intercept_zq_ust_field = zeros(N_Sites,1); %intercept of fit
sigmaintercept_zq_ust_field = zeros(N_Sites,1); %uncertainty in intercept of fit
slope_zq_ust_field = zeros(N_Sites,1); %slope of fit for zq versus ustar
sigmaslope_zq_ust_field = zeros(N_Sites,1); %uncertainty in slope for zq versus ustar
zq_bar_field = zeros(N_Sites,1); %mean zq for field sites

for i = 1:N_Sites
    %fit values for zq versus ust
    ind_fit = find(fQ_ust_bin_avg{i}>fQ_zQ_ust);
    ust_fit = ust_bin_avg{i}(ind_fit);
    sigma_ust_fit = ust_bin_SE{i}(ind_fit);
    zq_fit = zq_ust_bin_avg{i}(ind_fit);
    sigma_zq_fit = zq_ust_bin_SE{i}(ind_fit);

    [a, b, sigma_a, sigma_b, zq_pred, sigma_zq_pred] = linearfit(ust_fit, zq_fit, sigma_ust_fit, sigma_zq_fit);
    intercept_zq_ust_field(i) = a; %intercept of fit
    sigmaintercept_zq_ust_field(i) = sigma_a; %uncertainty in intercept of fit
    slope_zq_ust_field(i) = b; %slope of fit
    sigmaslope_zq_ust_field(i) = sigma_b; %uncertainty in slope of fit
    zq_bar_field(i) = mean(zq_ust_bin_avg{i}(ind_fit));
    
    plot(ust_fit,zq_pred,'Color',Colors_field{i},'LineWidth',LineWidth_field);
end

%fit values for zq versus ust - literature
intercept_zq_ust_lit = zeros(N_Lit,1); %intercept of fit
sigmaintercept_zq_ust_lit = zeros(N_Lit,1); %uncertainty in intercept of fit
slope_zq_ust_lit = zeros(N_Lit,1); %slope of fit for zq versus ustar
sigmaslope_zq_ust_lit = zeros(N_Lit,1); %uncertainty in slope for zq versus ustar

for i = 1:N_Lit
    [a, b, sigma_a, sigma_b, zq_pred, sigma_zq_pred] = linearfit(ust_lit{i},zq_lit{i},sigma_ust_lit{i},sigma_zq_lit{i});
    intercept_zq_ust_lit(i) = a; %intercept of fit
    sigmaintercept_zq_ust_lit(i) = sigma_a; %uncertainty in intercept of fit
    slope_zq_ust_lit(i) = b; %slope of fit
    sigmaslope_zq_ust_lit(i) = sigma_b; %uncertainty in slope of fit
    
    plot(ust_lit{i},zq_pred,'Color',Colors_lit{i},'LineWidth',LineWidth_lit);
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
    ind_plot = find(fQ_ust_bin_avg{i}>fQ_zQ_ust);
    errorbar(ust_bin_avg{i}(ind_plot),zqnorm_ust_bin_avg{i}(ind_plot),zqnorm_ust_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot literature data
for i = 1:N_Lit        
    errorbar(ust_lit{i},zqnorm_lit{i},sigma_zqnorm_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end

%fit values for zqnorm versus ust
slope_zqnorm_ust_field_all = zeros(N_Sites,1); %slope of fit for zqnorm versus ustar
sigmaslope_zqnorm_ust_field_all = zeros(N_Sites,1); %uncertainty in slope for zqnorm versus ustar
zqnorm_pred_field{i} = cell(N_Sites,1); %predicted zqnorm
zqnorm_bar_field_all = zeros(N_Sites,1); %mean of zqnorm
sigma_zqnorm_bar_field_all = zeros(N_Sites,1); %uncertainty of zqnorm bar
zqnorm_pred_relativerange_field_all = zeros(N_Sites,1); %relative range of zqnorm prediction


for i = 1:N_Sites
    %fit values for zqnorm versus ust
    ind_fit = find(fQ_ust_bin_avg{i}>fQ_zQ_ust);
    ust_fit = ust_bin_avg{i}(ind_fit);
    sigma_ust_fit = ust_bin_SE{i}(ind_fit);
    zqnorm_fit = zqnorm_ust_bin_avg{i}(ind_fit);
    sigma_zqnorm_fit = zqnorm_ust_bin_SE{i}(ind_fit);

    [~, b, ~, sigma_b, zqnorm_pred, ~] = linearfit(ust_fit, zqnorm_fit, sigma_ust_fit, sigma_zqnorm_fit);
    slope_zqnorm_ust_field_all(i) = b; %slope of fit
    sigmaslope_zqnorm_ust_field_all(i) = sigma_b; %uncertainty in slope of fit
    zqnorm_pred_field{i} = zqnorm_pred;
    %zqnorm_bar_field_all(i) = mean(zqnorm_ust_bin_avg{i}(ind_fit));
    %sigma_zqnorm_bar_field_all(i) = std(zqnorm_ust_bin_avg{i}(ind_fit))/sqrt(length(ind_fit));
    [zqnorm_bar_field, sigma_zqnorm_bar_field] = MeanUncertainty(zqnorm_ust_bin_avg{i}(ind_fit), zqnorm_ust_bin_SE{i}(ind_fit));
    zqnorm_bar_field_all(i) = zqnorm_bar_field;
    sigma_zqnorm_bar_field_all(i) = sigma_zqnorm_bar_field;
    zqnorm_pred_relativerange_field_all(i) = (zqnorm_pred(end)-zqnorm_pred(1))/zqnorm_bar_field_all(i);
    
    % plot fit value of zqnorm
    plot(ust_fit,zqnorm_pred_field{i},'Color',Colors_field{i},'LineWidth',LineWidth_field); 
end

%fit values for zqnorm versus ust - literature
slope_zqnorm_ust_lit = zeros(N_Lit,1); %slope of fit for zqnorm versus ustar
sigmaslope_zqnorm_ust_lit = zeros(N_Lit,1); %uncertainty in slope for zqnorm versus ustar
zqnorm_pred_lit{i} = cell(N_Lit,1); %predicted zqnorm
zqnorm_bar_lit_all = zeros(N_Lit,1); %mean of zqnorm
sigma_zqnorm_bar_lit_all = zeros(N_Lit,1); %uncertainty of zqnorm bar
zqnorm_pred_relativerange_lit_all = zeros(N_Lit,1); %relative range of zqnorm prediction


for i = 1:N_Lit 
    %fit values for zqnorm versus ust
    [a, b, sigma_a, sigma_b, zqnorm_pred, sigma_zqnorm_pred] = linearfit(ust_lit{i},zqnorm_lit{i},sigma_ust_lit{i},sigma_zqnorm_lit{i});
    slope_zqnorm_ust_lit(i) = b; %slope of fit
    sigmaslope_zqnorm_ust_lit(i) = sigma_b; %uncertainty in slope of fit
    zqnorm_pred_lit{i} = zqnorm_pred;
    %zqnorm_bar_lit_all(i) = mean(zqnorm_lit{i});
    %sigma_zqnorm_bar_lit_all(i) = std(zqnorm_lit{i})/sqrt(length(zqnorm_lit));
    [zqnorm_bar_lit, sigma_zqnorm_bar_lit] = MeanUncertainty(zqnorm_lit{i}, sigma_zqnorm_lit{i});
    zqnorm_bar_lit_all(i) = zqnorm_bar_lit;
    sigma_zqnorm_bar_lit_all(i) = sigma_zqnorm_bar_lit;
    
    zqnorm_pred_relativerange_lit_all(i) = (zqnorm_pred(ust_lit{i}==max(ust_lit{i}))-zqnorm_pred(ust_lit{i}==min(ust_lit{i})))/zqnorm_bar_lit_all(i);
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

%print plot for draft
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
%print([folder_Plots,'zq_ust_',int2str(fQ_zQ_ust*100),'.png'],'-dpng');
print([folder_Plots,'zq_ust.png'],'-dpng');


%% replot raw data
subplot(2,1,1); cla;
subplot(2,1,2); cla;
for i = 1:N_Sites
    ind_plot_all = find(fQ_all{i}>fQ_zQ_ust); %get all possible points to plot based on flux frequency
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    %plot in both subplot windows
    subplot(2,1,1);
    errorbar(ustRe_all{i}(ind_plot),zq_all{i}(ind_plot),sigma_zq_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    subplot(2,1,2);
    errorbar(ustRe_all{i}(ind_plot),zqnorm_all{i}(ind_plot),sigma_zqnorm_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end

for i = 2 %replot first two sites for clarity
    ind_plot_all = find(fQ_all{i}>fQ_zQ_ust); %get all possible points to plot based on flux frequency
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    %plot in both subplot windows
    subplot(2,1,1);
    errorbar(ustRe_all{i}(ind_plot),zq_all{i}(ind_plot),sigma_zq_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    subplot(2,1,2);
    errorbar(ustRe_all{i}(ind_plot),zqnorm_all{i}(ind_plot),sigma_zqnorm_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end

%redo ylim and legend
subplot(2,1,2);
ylim([0 300]);
legend(Sites,'Location','SouthEast');

%print plot for draft
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
%print([folder_Plots,'zq_ust_raw_',int2str(fQ_zQ_ust*100),'.png'],'-dpng');
print([folder_Plots,'zq_ust_raw.png'],'-dpng');


%% FOR EACH SITE, PERFORM BINNING BY TAU
%create initial tau bins
tau_bin_edges_init = tau_bin_lower:tau_bin_increment:tau_bin_upper;
tau_bin_min_init = tau_bin_edges_init(1:end-1);
tau_bin_max_init = tau_bin_edges_init(2:end);
N_tau_bins_init = length(tau_bin_edges_init)-1;

%initialize lists of tau bins by site
tau_bin_values = cell(N_Sites,1); %tau's in each bin
N_tau_bin_values = cell(N_Sites,1); %N tau's in each bin
tau_bin_avg = cell(N_Sites,1); %avg tau in bin
tau_bin_std = cell(N_Sites,1); %std tau in bin
tau_bin_SE = cell(N_Sites,1); %SE tau in bin

%separate fQs into tau bins
fQ_tau_bin_values = cell(N_Sites,1); %fQs into tau bins
fQ_tau_bin_avg = cell(N_Sites,1); %get average fQ for tau bins
fQ_tau_bin_std = cell(N_Sites,1); %get std fQ for tau bins
fQ_tau_bin_SE = cell(N_Sites,1); %get SE fQ for tau bins

%separate zq's into tau bins
zq_tau_bin_values = cell(N_Sites,1); %zq into tau bins
zq_tau_bin_avg = cell(N_Sites,1); %get average zq for tau bins
zq_tau_bin_std = cell(N_Sites,1); %get std zq for tau bins
zq_tau_bin_SE = cell(N_Sites,1); %get SE zq for tau bins

%separate zqnorm's (zq/d50) into tau bins
zqnorm_tau_bin_values = cell(N_Sites,1); %zq/d50s into tau bins
zqnorm_tau_bin_avg = cell(N_Sites,1); %get average zq/d50 for tau bins
zqnorm_tau_bin_std = cell(N_Sites,1); %get std zq/d50 for tau bins
zqnorm_tau_bin_SE = cell(N_Sites,1); %get SE zq/d50 for tau bins

%separate Qs into tau bins
Q_tau_bin_values = cell(N_Sites,1); %Qs into tau bins
Q_tau_bin_avg = cell(N_Sites,1); %get average Q for tau bins
Q_tau_bin_std = cell(N_Sites,1); %get std Q for tau bins
Q_tau_bin_SE = cell(N_Sites,1); %get SE Q for tau bins

%go through all sites
for i = 1:N_Sites
    
    %get tau bins
    bin_N_init = histcounts(tauRe_all{i},tau_bin_edges_init); %get number of tau in each initial bin
    tau_bin_min_Site = []; %initialize list of tau bin mins for site
    tau_bin_max_Site = []; %initialize list of tau bin maxes for site
    cum_N = 0; %initialize cumulative tau counts for site bin
    tau_bin_min_this = tau_bin_min_init(1); %initialize tau bin min for this bin
    for j = 1:N_tau_bins_init; %go through all initial tau bins to determine tau bins for site
        cum_N = cum_N+bin_N_init(j); %add to cumulative tau counts
        if cum_N>=bin_N_min; %if cumulative tau counts exceeds limit
            tau_bin_max_this = tau_bin_max_init(j); %set tau bin max for this bin
            tau_bin_min_Site = [tau_bin_min_Site, tau_bin_min_this]; %add to list of tau bin mins for Site
            tau_bin_max_Site = [tau_bin_max_Site, tau_bin_max_this]; %add to list of tau bin maxes for Site
            cum_N = 0; %reset cumulative tau counts for site bin
            tau_bin_min_this = tau_bin_max_init(j); %reset tau bin min for next bin
        end
    end
    N_tau_bins_Site = length(tau_bin_min_Site); %get number of tau bins for Site
   
    %initialize tau bins
    tau_bin_values{i} = cell(N_tau_bins_Site,1);
    N_tau_bin_values{i} = zeros(N_tau_bins_Site,1)*NaN;
    tau_bin_avg{i} = zeros(N_tau_bins_Site,1)*NaN;
    tau_bin_std{i} = zeros(N_tau_bins_Site,1)*NaN;
    tau_bin_SE{i} = zeros(N_tau_bins_Site,1)*NaN;

    %initialize fQ's into tau bins
    fQ_tau_bin_values{i} = cell(N_tau_bins_Site,1);
    fQ_tau_bin_avg{i} = zeros(N_tau_bins_Site,1)*NaN;
    fQ_tau_bin_std{i} = zeros(N_tau_bins_Site,1)*NaN;
    fQ_tau_bin_SE{i} = zeros(N_tau_bins_Site,1)*NaN;
    
    %initialize zq's into tau bins
    zq_tau_bin_values{i} = cell(N_tau_bins_Site,1);
    zq_tau_bin_avg{i} = zeros(N_tau_bins_Site,1)*NaN;
    zq_tau_bin_std{i} = zeros(N_tau_bins_Site,1)*NaN;
    zq_tau_bin_SE{i} = zeros(N_tau_bins_Site,1)*NaN;
    
    %initialize zqnorm's into tau bins
    zqnorm_tau_bin_values{i} = cell(N_tau_bins_Site,1);
    zqnorm_tau_bin_avg{i} = zeros(N_tau_bins_Site,1)*NaN;
    zqnorm_tau_bin_std{i} = zeros(N_tau_bins_Site,1)*NaN;
    zqnorm_tau_bin_SE{i} = zeros(N_tau_bins_Site,1)*NaN;
    
    %initialize Q's into tau bins
    Q_tau_bin_values{i} = cell(N_tau_bins_Site,1);
    Q_tau_bin_avg{i} = zeros(N_tau_bins_Site,1)*NaN;
    Q_tau_bin_std{i} = zeros(N_tau_bins_Site,1)*NaN;
    Q_tau_bin_SE{i} = zeros(N_tau_bins_Site,1)*NaN;
    
    %get values, avgs, and standard deviations for tau bins
    for j = 1:N_tau_bins_Site
        %get indices
        bin_ind = find(tauRe_all{i}>=tau_bin_min_Site(j)&tauRe_all{i}<=tau_bin_max_Site(j));

        %get tau in bin
        tau_bin_values{i}{j} = tauRe_all{i}(bin_ind);
        N_tau_bin_values{i}(j) = length(tau_bin_values{i}{j});
        tau_bin_avg{i}(j) = mean(tau_bin_values{i}{j});
        tau_bin_std{i}(j) = std(tau_bin_values{i}{j});
        tau_bin_SE{i}(j) = tau_bin_std{i}(j)./sqrt(N_tau_bin_values{i}(j));
        
        %get fQs for tau bins
        fQ_tau_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_tau_bin_avg{i}(j) = mean(fQ_tau_bin_values{i}{j});
        fQ_tau_bin_std{i}(j) = std(fQ_tau_bin_values{i}{j});
        fQ_tau_bin_SE{i}(j) = fQ_tau_bin_std{i}(j)./sqrt(N_tau_bin_values{i}(j));
        
        %get zqs for tau bins
        zq_tau_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_tau_bin_avg{i}(j) = mean(zq_tau_bin_values{i}{j});
        zq_tau_bin_std{i}(j) = std(zq_tau_bin_values{i}{j});
        zq_tau_bin_SE{i}(j) = zq_tau_bin_std{i}(j)./sqrt(N_tau_bin_values{i}(j));
        
        %get zqnorms for tau bins
        zqnorm_tau_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_tau_bin_avg{i}(j) = mean(zqnorm_tau_bin_values{i}{j});
        zqnorm_tau_bin_std{i}(j) = std(zqnorm_tau_bin_values{i}{j});
        zqnorm_tau_bin_SE{i}(j) = zqnorm_tau_bin_std{i}(j)./sqrt(N_tau_bin_values{i}(j));
        
        %get Qs for tau bins
        Q_tau_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_tau_bin_avg{i}(j) = mean(Q_tau_bin_values{i}{j});
        Q_tau_bin_std{i}(j) = std(Q_tau_bin_values{i}{j});
        Q_tau_bin_SE{i}(j) = Q_tau_bin_std{i}(j)./sqrt(N_tau_bin_values{i}(j));
    end
end


%% PERFORM FIT OF Q VERSUS TAU

%initialize tau bins with transport
tau_transport = cell(N_Sites,1);
sigma_tau_transport = cell(N_Sites,1);
Q_transport = cell(N_Sites,1);
sigma_Q_transport = cell(N_Sites,1);

%initialize fit values
C_Qtau_fit_all = zeros(N_Sites,1);
sigma_C_Qtau_fit_all = zeros(N_Sites,1);
tauit_fit_all = zeros(N_Sites,1);
sigma_tauit_fit_all = zeros(N_Sites,1);
ustit_fit_all = zeros(N_Sites,1);
sigma_ustit_fit_all = zeros(N_Sites,1);
chi2_Qtau_all = zeros(N_Sites,1);
p_Qtau_all = zeros(N_Sites,1);

for i = 1:N_Sites
    %get tau bins with transport
    tau_bin_ind = find(fQ_tau_bin_avg{i}>fQ_Q_tau);
    tau_transport{i} = tau_bin_avg{i}(tau_bin_ind);
    sigma_tau_transport{i} = tau_bin_SE{i}(tau_bin_ind);
    Q_transport{i} = Q_tau_bin_avg{i}(tau_bin_ind);
    sigma_Q_transport{i} = Q_tau_bin_SE{i}(tau_bin_ind);
    
    %get fit values
    [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(tau_transport{i}, Q_transport{i}, sigma_tau_transport{i}, sigma_Q_transport{i});
    C_Qtau_fit_all(i) = b;
    sigma_C_Qtau_fit_all(i) = sigma_b;
    tauit_fit_all(i) = -a/C_Qtau_fit_all(i);
    sigma_tauit_fit_all(i) = sigma_a/C_Qtau_fit_all(i);
    ustit_fit_all(i) = sqrt(tauit_fit_all(i)/rho_a);
    sigma_ustit_fit_all(i) = 0.5*sigma_tauit_fit_all(i)*sqrt(rho_a/tauit_fit_all(i));
    
    %compute chi2 and p-value
    sigma_Q_chi2 = sqrt(sigma_Q_transport{i}.^2 + b^2*sigma_tau_transport{i}.^2); %get sigma_Q values that combine uncertainties in Q and tau
    [chi2, p] = Chi2Calculation(Q_transport{i}, sigma_Q_chi2, Q_pred, length(Q_pred-2));
    chi2_Qtau_all(i) = chi2;
    p_Qtau_all(i) = p;
end


%% PLOT Q VERSUS TAU
figure; clf; hold on;
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    errorbar(tau_bin_avg{i},Q_tau_bin_avg{i},Q_tau_bin_SE{i},Markers_field{i},'MarkerSize',10,'Color',Colors_field{i});
    plot([tauit_fit_all(i) max(tau_transport{i})],[0 C_Qtau_fit_all(i)*(max(tau_transport{i})-tauit_fit_all(i))],'Color',Colors_field{i}); %plot fit for these values
    legend_items{2*i-1} = Sites{i}; %add to list of legend items
    legend_items{2*i} = 'fit'; %add to list of legend items
end
xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{$$Q$$ (g/m/s)}','Interpreter','Latex');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',PlotFont);
set(gcf, 'PaperPosition',[0 0 8 5]);
%print([folder_Plots,'Q_tau_',int2str(fQ_Q_tau*100),'.png'],'-dpng');
print([folder_Plots,'Q_tau.png'],'-dpng');


%% replot raw data
cla;
for i = 1:N_Sites
    ind_plot_all = find(fQ_all{i}>fQ_Q_tau); %get all possible points to plot based on fQ cutoff
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    errorbar(tauRe_all{i}(ind_plot),Q_all{i}(ind_plot),sigma_Q_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end
%replot first two sites for clariy
for i = 1:2
    ind_plot_all = find(fQ_all{i}>fQ_Q_tau); %get all possible points to plot based on fQ cutoff
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    errorbar(tauRe_all{i}(ind_plot),Q_all{i}(ind_plot),sigma_Q_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end
legend(Sites,'Location','NorthWest');

%print plot for draft
set(gcf, 'PaperPosition',[0 0 8 5]);
%print([folder_Plots,'Q_tau_raw_',int2str(fQ_Q_tau*100),'.png'],'-dpng');
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

%create initial tauex bins
tauex_bin_edges_init = tauex_bin_lower:tauex_bin_increment:tauex_bin_upper;
tauex_bin_min_init = tauex_bin_edges_init(1:end-1);
tauex_bin_max_init = tauex_bin_edges_init(2:end);
N_tauex_bins_init = length(tauex_bin_edges_init)-1;

%initialize lists of tauex bins by site
tauex_bin_values = cell(N_Sites,1); %tauex's in each bin
N_tauex_bin_values = cell(N_Sites,1); %N tauex's in each bin
tauex_bin_avg = cell(N_Sites,1); %avg tauex in bin
tauex_bin_std = cell(N_Sites,1); %std tauex in bin
tauex_bin_SE = cell(N_Sites,1); %SE tauex in bin
tauex_bin_min = cell(N_Sites,1); %minimum tauex in bin

%separate fQs into tauex bins
fQ_tauex_bin_values = cell(N_Sites,1); %fQs into tauex bins
fQ_tauex_bin_avg = cell(N_Sites,1); %get average fQ for tauex bins
fQ_tauex_bin_std = cell(N_Sites,1); %get std fQ for tauex bins
fQ_tauex_bin_SE = cell(N_Sites,1); %get SE fQ for tauex bins

%separate zq's into tauex bins
zq_tauex_bin_values = cell(N_Sites,1); %zq into tauex bins
zq_tauex_bin_avg = cell(N_Sites,1); %get average zq for tauex bins
zq_tauex_bin_std = cell(N_Sites,1); %get std zq for tauex bins
zq_tauex_bin_SE = cell(N_Sites,1); %get SE zq for tauex bins

%separate zqnorm's (zq/d50) into tauex bins
zqnorm_tauex_bin_values = cell(N_Sites,1); %zq/d50s into tauex bins
zqnorm_tauex_bin_avg = cell(N_Sites,1); %get average zq/d50 for tauex bins
zqnorm_tauex_bin_std = cell(N_Sites,1); %get std zq/d50 for tauex bins
zqnorm_tauex_bin_SE = cell(N_Sites,1); %get SE zq/d50 for tauex bins

%separate Qs into tauex bins
Q_tauex_bin_values = cell(N_Sites,1); %Qs into tauex bins
Q_tauex_bin_avg = cell(N_Sites,1); %get average Q for tauex bins
Q_tauex_bin_std = cell(N_Sites,1); %get std Q for tauex bins
Q_tauex_bin_SE = cell(N_Sites,1); %get SE Q for tauex bins

%separate Qhats into tauex bins
Qhat_tauex_bin_values = cell(N_Sites,1); %Qhats into tauex bins
Qhat_tauex_bin_avg = cell(N_Sites,1); %get average Qhat for tauex bins
Qhat_tauex_bin_std = cell(N_Sites,1); %get std Qhat for tauex bins
Qhat_tauex_bin_SE = cell(N_Sites,1); %get SE Qhat for tauex bins

%go through all sites
for i = 1:N_Sites
    
    %get tauex bins
    bin_N_init = histcounts(tauex_all{i},tauex_bin_edges_init); %get number of tauex in each initial bin
    tauex_bin_min_Site = []; %initialize list of tauex bin mins for site
    tauex_bin_max_Site = []; %initialize list of tauex bin maxes for site
    cum_N = 0; %initialize cumulative tauex counts for site bin
    tauex_bin_min_this = tauex_bin_min_init(1); %initialize tauex bin min for this bin
    for j = 1:N_tauex_bins_init; %go through all initial tauex bins to determine tauex bins for site
        cum_N = cum_N+bin_N_init(j); %add to cumulative tauex counts
        if cum_N>=bin_N_min; %if cumulative tauex counts exceeds limit
            tauex_bin_max_this = tauex_bin_max_init(j); %set tauex bin max for this bin
            tauex_bin_min_Site = [tauex_bin_min_Site, tauex_bin_min_this]; %add to list of tauex bin mins for Site
            tauex_bin_max_Site = [tauex_bin_max_Site, tauex_bin_max_this]; %add to list of tauex bin maxes for Site
            cum_N = 0; %reset cumulative tauex counts for site bin
            tauex_bin_min_this = tauex_bin_max_init(j); %reset tauex bin min for next bin
        end
    end
    N_tauex_bins_Site = length(tauex_bin_min_Site); %get number of tauex bins for Site
   
    %initialize tauex bins
    tauex_bin_values{i} = cell(N_tauex_bins_Site,1);
    N_tauex_bin_values{i} = zeros(N_tauex_bins_Site,1)*NaN;
    tauex_bin_avg{i} = zeros(N_tauex_bins_Site,1)*NaN;
    tauex_bin_std{i} = zeros(N_tauex_bins_Site,1)*NaN;
    tauex_bin_SE{i} = zeros(N_tauex_bins_Site,1)*NaN;
    tauex_bin_min{i} = zeros(N_tauex_bins_Site,1)*NaN;

    %initialize fQ's into tauex bins
    fQ_tauex_bin_values{i} = cell(N_tauex_bins_Site,1);
    fQ_tauex_bin_avg{i} = zeros(N_tauex_bins_Site,1)*NaN;
    fQ_tauex_bin_std{i} = zeros(N_tauex_bins_Site,1)*NaN;
    fQ_tauex_bin_SE{i} = zeros(N_tauex_bins_Site,1)*NaN;
    
    %initialize zq's into tauex bins
    zq_tauex_bin_values{i} = cell(N_tauex_bins_Site,1);
    zq_tauex_bin_avg{i} = zeros(N_tauex_bins_Site,1)*NaN;
    zq_tauex_bin_std{i} = zeros(N_tauex_bins_Site,1)*NaN;
    zq_tauex_bin_SE{i} = zeros(N_tauex_bins_Site,1)*NaN;
    
    %initialize zqnorm's into tauex bins
    zqnorm_tauex_bin_values{i} = cell(N_tauex_bins_Site,1);
    zqnorm_tauex_bin_avg{i} = zeros(N_tauex_bins_Site,1)*NaN;
    zqnorm_tauex_bin_std{i} = zeros(N_tauex_bins_Site,1)*NaN;
    zqnorm_tauex_bin_SE{i} = zeros(N_tauex_bins_Site,1)*NaN;
    
    %initialize Q's into tauex bins
    Q_tauex_bin_values{i} = cell(N_tauex_bins_Site,1);
    Q_tauex_bin_avg{i} = zeros(N_tauex_bins_Site,1)*NaN;
    Q_tauex_bin_std{i} = zeros(N_tauex_bins_Site,1)*NaN;
    Q_tauex_bin_SE{i} = zeros(N_tauex_bins_Site,1)*NaN;

    %initialize Qhat's into tauex bins
    Qhat_tauex_bin_values{i} = cell(N_tauex_bins_Site,1);
    Qhat_tauex_bin_avg{i} = zeros(N_tauex_bins_Site,1)*NaN;
    Qhat_tauex_bin_std{i} = zeros(N_tauex_bins_Site,1)*NaN;
    Qhat_tauex_bin_SE{i} = zeros(N_tauex_bins_Site,1)*NaN;
    
    %get values, avgs, and standard deviations for tauex bins
    for j = 1:N_tauex_bins_Site
        %get indices
        bin_ind = find(tauex_all{i}>=tauex_bin_min_Site(j)&tauex_all{i}<=tauex_bin_max_Site(j));

        %get tauex in bin
        tauex_bin_values{i}{j} = tauex_all{i}(bin_ind);
        N_tauex_bin_values{i}(j) = length(tauex_bin_values{i}{j});
        tauex_bin_avg{i}(j) = mean(tauex_bin_values{i}{j});
        tauex_bin_std{i}(j) = std(tauex_bin_values{i}{j});
        tauex_bin_SE{i}(j) = tauex_bin_std{i}(j)./sqrt(N_tauex_bin_values{i}(j));
        tauex_bin_min{i}(j) = min(tauex_bin_values{i}{j});
        
        %get fQ for tauex bins
        fQ_tauex_bin_values{i}{j} = fQ_all{i}(bin_ind);
        fQ_tauex_bin_avg{i}(j) = mean(fQ_tauex_bin_values{i}{j});
        fQ_tauex_bin_std{i}(j) = std(fQ_tauex_bin_values{i}{j});
        fQ_tauex_bin_SE{i}(j) = fQ_tauex_bin_std{i}(j)./sqrt(N_tauex_bin_values{i}(j));
        
        %get zqs for tauex bins
        zq_tauex_bin_values{i}{j} = zq_all{i}(bin_ind);
        zq_tauex_bin_avg{i}(j) = mean(zq_tauex_bin_values{i}{j});
        zq_tauex_bin_std{i}(j) = std(zq_tauex_bin_values{i}{j});
        zq_tauex_bin_SE{i}(j) = zq_tauex_bin_std{i}(j)./sqrt(N_tauex_bin_values{i}(j));
        
        %get zqnorms for tauex bins
        zqnorm_tauex_bin_values{i}{j} = zqnorm_all{i}(bin_ind);
        zqnorm_tauex_bin_avg{i}(j) = mean(zqnorm_tauex_bin_values{i}{j});
        zqnorm_tauex_bin_std{i}(j) = std(zqnorm_tauex_bin_values{i}{j});
        zqnorm_tauex_bin_SE{i}(j) = zqnorm_tauex_bin_std{i}(j)./sqrt(N_tauex_bin_values{i}(j));
        
        %get Qs for tauex bins
        Q_tauex_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_tauex_bin_avg{i}(j) = mean(Q_tauex_bin_values{i}{j});
        Q_tauex_bin_std{i}(j) = std(Q_tauex_bin_values{i}{j});
        Q_tauex_bin_SE{i}(j) = Q_tauex_bin_std{i}(j)./sqrt(N_tauex_bin_values{i}(j));
        
        %get Qhats for tauex bins
        Qhat_tauex_bin_values{i}{j} = Qhat_all{i}(bin_ind);
        Qhat_tauex_bin_avg{i}(j) = mean(Qhat_tauex_bin_values{i}{j});
        Qhat_tauex_bin_std{i}(j) = std(Qhat_tauex_bin_values{i}{j});
        Qhat_tauex_bin_SE{i}(j) = Qhat_tauex_bin_std{i}(j)./sqrt(N_tauex_bin_values{i}(j));
    end
end

%compute Qhat_bin from averaged values (instead of raw values)
for i=1:N_Sites
    Qhat_tauex_bin_avg{i} = Q_tauex_bin_avg{i}./((1/g)*ustit_fit_all(i).*1e3*tauex_bin_avg{i});
    Qhat_tauex_bin_SE{i} = Q_tauex_bin_SE{i}./((1/g)*ustit_fit_all(i).*1e3*tauex_bin_avg{i});
end


%% FITS TO Q AND QHAT VERSUS TAUEX

%scaling coefficient for Q versus tauex
C_Qtauex_fit_all = zeros(N_Sites,1);

%scaling coefficient for Qhat versus tauex
Qhat_bar_all = zeros(N_Sites,1); 
sigma_Qhat_bar_all = zeros(N_Sites,1);

%perform fits
for i = 1:N_Sites
    %set tauex bins 
    tauex_bin_ind = intersect(find(fQ_tauex_bin_avg{i}>fQ_Qhat_tauex),... %must have sufficient transport
        find(tauex_bin_min{i}>=N_sigma_tauit*sigma_tauit_fit_all(i))); %must exceed threshold
    
    %get values for fitting
    tauex_fit = tauex_bin_avg{i}(tauex_bin_ind);
    Q_fit = Q_tauex_bin_avg{i}(tauex_bin_ind);
    Qhat_fit = Qhat_tauex_bin_avg{i}(tauex_bin_ind);
    sigma_Qhat_fit = Qhat_tauex_bin_SE{i}(tauex_bin_ind);
    
    %compute fit values
    C_Qtauex_fit_all(i) = mean(Q_fit./tauex_fit); %for Q versus tauex
    %Qhat_bar_all(i) = mean(Qhat_fit); %mean of Qhat
    %sigma_Qhat_bar_all(i) = std(Qhat_fit)/sqrt(length(Qhat_fit));
    %sigma_Qhat_bar_all(i) = (Qhat_bar_all(i)/length(tauex_fit))*sqrt(sum(sigma_Qhat_fit.^2)); %uncertainty for mean of Qhat
    [Qhat_bar, sigma_Qhat_bar] = MeanUncertainty(Qhat_fit, sigma_Qhat_fit);
    Qhat_bar_all(i) = Qhat_bar;
    sigma_Qhat_bar_all(i) = sigma_Qhat_bar;
end


%% PLOT Q VERSUS TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    ind_plot_fQ = find(fQ_tauex_bin_avg{i}>fQ_Q_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_bin_min{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot = intersect(ind_plot_fQ,ind_plot_tauex);
    errorbar(tauex_bin_avg{i}(ind_plot),Q_tauex_bin_avg{i}(ind_plot),Q_tauex_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit, C_Qtauex_fit_all(i)*tauex_fit,'Color',Colors_field{i});
end

xlim([0 0.3]);
ylim([0 60]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%draft paper plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
%print([folder_Plots,'Q_tauex_',int2str(fQ_Q_tauex*100),'.png'],'-dpng');
print([folder_Plots,'Q_tauex.png'],'-dpng');


%% replot raw data
cla;
%for i = 1:2
for i = 1:N_Sites
    ind_plot_fQ = find(fQ_all{i}>fQ_Q_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot_all = intersect(ind_plot_fQ,ind_plot_tauex); %these are all possible points to plot
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    errorbar(tauex_all{i}(ind_plot),Q_all{i}(ind_plot),sigma_Q_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end

%replot first two sites for clarity
for i = 1:2
    ind_plot_fQ = find(fQ_all{i}>fQ_Q_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot_all = intersect(ind_plot_fQ,ind_plot_tauex); %these are all possible points to plot
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    errorbar(tauex_all{i}(ind_plot),Q_all{i}(ind_plot),sigma_Q_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end

legend(Sites,'Location','NorthWest');

%print plot for draft
set(gcf, 'PaperPosition',[0 0 6.5 5]);
%print([folder_Plots,'Q_tauex_raw_',int2str(fQ_Q_tauex*100),'.png'],'-dpng');
print([folder_Plots,'Q_tauex_raw.png'],'-dpng');

%% PLOT QHAT VS TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    ind_plot_fQ = find(fQ_tauex_bin_avg{i}>fQ_Q_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_bin_min{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot = intersect(ind_plot_fQ,ind_plot_tauex);
    errorbar(tauex_bin_avg{i}(ind_plot),Qhat_tauex_bin_avg{i}(ind_plot),Qhat_tauex_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%replot binned values for first two sites for clarity
for i = 1:2
    ind_plot_fQ = find(fQ_tauex_bin_avg{i}>fQ_Q_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_bin_min{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot = intersect(ind_plot_fQ,ind_plot_tauex);
    errorbar(tauex_bin_avg{i}(ind_plot),Qhat_tauex_bin_avg{i}(ind_plot),Qhat_tauex_bin_SE{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit,Qhat_bar_all(i)*ones(2,1),'Color',Colors_field{i});
end
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Dimensionless saltation flux, $$\hat{Q}$$}','Interpreter','Latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
legend(Sites,'Location','SouthEast');
xlim([0 0.3]);
ylim([0 9]);
set(gca,'FontSize',PlotFont);

%draft plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
%print([folder_Plots,'Qhat_tauex_',int2str(fQ_Qhat_tauex*100),'.png'],'-dpng');
print([folder_Plots,'Qhat_tauex.png'],'-dpng');

%% replot raw data
cla;
%for i = 1:2
for i = 1:N_Sites
    ind_plot_fQ = find(fQ_all{i}>fQ_Qhat_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot_all = intersect(ind_plot_fQ,ind_plot_tauex); %these are all possible points to plot
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    %now plot
    errorbar(tauex_all{i}(ind_plot),Qhat_all{i}(ind_plot),sigma_Qhat_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end

%replot first two sites for clarity
for i = 1:2
    ind_plot_fQ = find(fQ_all{i}>fQ_Qhat_tauex); %get tauex bins with sufficient transport
    ind_plot_tauex = find(tauex_all{i}>=N_sigma_tauit*sigma_tauit_fit_all(i)); %get tauex bins that exceed threshold
    ind_plot_all = intersect(ind_plot_fQ,ind_plot_tauex); %these are all possible points to plot
    
    %reduce it so that time intervals are non-overlapping
    ind_plot = ind_plot_all(1);
    for j = 2:length(ind_plot_all)
        if StartTime_all{i}(ind_plot_all(j))>=EndTime_all{i}(ind_plot(end))
            ind_plot = [ind_plot; ind_plot_all(j)];
        end
    end
    
    %now plot
    errorbar(tauex_all{i}(ind_plot),Qhat_all{i}(ind_plot),sigma_Qhat_all{i}(ind_plot),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
end

legend(Sites,'Location','NorthEast');
ylim([0 15]);

%print plot for draft
set(gcf, 'PaperPosition',[0 0 6.5 5]);
%print([folder_Plots,'Qhat_tauex_raw_',int2str(fQ_Qhat_tauex*100),'.png'],'-dpng');
print([folder_Plots,'Qhat_tauex_raw.png'],'-dpng');

%% FOR EACH SITE, PERFORM BINNING BY UST*TAUEX
%compute ust*tauex
usttauex_all = cell(N_Sites,1);
for i=1:N_Sites
    usttauex_all{i} = ustRe_all{i}.*tauex_all{i};
end

%create initial tauex bins
usttauex_bin_lower = 0;
usttauex_bin_upper = 0.145;
usttauex_bin_increment = 0.005;
usttauex_bin_edges_init = usttauex_bin_lower:usttauex_bin_increment:usttauex_bin_upper;
usttauex_bin_min_init = usttauex_bin_edges_init(1:end-1);
usttauex_bin_max_init = usttauex_bin_edges_init(2:end);
N_usttauex_bins_init = length(usttauex_bin_edges_init)-1;

%initialize lists of usttauex bins by site
usttauex_bin_values = cell(N_Sites,1); %usttauex's in each bin
N_usttauex_bin_values = cell(N_Sites,1); %N usttauex's in each bin
usttauex_bin_avg = cell(N_Sites,1); %avg usttauex in bin
usttauex_bin_std = cell(N_Sites,1); %std usttauex in bin
usttauex_bin_SE = cell(N_Sites,1); %SE usttauex in bin
usttauex_bin_min = cell(N_Sites,1); %minimum usttauex in bin

%separate Qs into usttauex bins
Q_usttauex_bin_values = cell(N_Sites,1); %Qs into usttauex bins
Q_usttauex_bin_avg = cell(N_Sites,1); %get average Q for usttauex bins
Q_usttauex_bin_std = cell(N_Sites,1); %get std Q for usttauex bins
Q_usttauex_bin_SE = cell(N_Sites,1); %get SE Q for usttauex bins

%go through all sites
for i = 1:N_Sites
    
    %get usttauex bins
    bin_N_init = histcounts(usttauex_all{i},usttauex_bin_edges_init); %get number of usttauex in each initial bin
    usttauex_bin_min_Site = []; %initialize list of usttauex bin mins for site
    usttauex_bin_max_Site = []; %initialize list of usttauex bin maxes for site
    cum_N = 0; %initialize cumulative usttauex counts for site bin
    usttauex_bin_min_this = usttauex_bin_min_init(1); %initialize usttauex bin min for this bin
    for j = 1:N_usttauex_bins_init; %go through all initial usttauex bins to determine usttauex bins for site
        cum_N = cum_N+bin_N_init(j); %add to cumulative usttauex counts
        if cum_N>=bin_N_min; %if cumulative usttauex counts exceeds limit
            usttauex_bin_max_this = usttauex_bin_max_init(j); %set usttauex bin max for this bin
            usttauex_bin_min_Site = [usttauex_bin_min_Site, usttauex_bin_min_this]; %add to list of usttauex bin mins for Site
            usttauex_bin_max_Site = [usttauex_bin_max_Site, usttauex_bin_max_this]; %add to list of usttauex bin maxes for Site
            cum_N = 0; %reset cumulative usttauex counts for site bin
            usttauex_bin_min_this = usttauex_bin_max_init(j); %reset usttauex bin min for next bin
        end
    end
    N_usttauex_bins_Site = length(usttauex_bin_min_Site); %get number of usttauex bins for Site
   
    %initialize usttauex bins
    usttauex_bin_values{i} = cell(N_usttauex_bins_Site,1);
    N_usttauex_bin_values{i} = zeros(N_usttauex_bins_Site,1)*NaN;
    usttauex_bin_avg{i} = zeros(N_usttauex_bins_Site,1)*NaN;
    usttauex_bin_std{i} = zeros(N_usttauex_bins_Site,1)*NaN;
    usttauex_bin_SE{i} = zeros(N_usttauex_bins_Site,1)*NaN;
    usttauex_bin_min{i} = zeros(N_usttauex_bins_Site,1)*NaN;

    %initialize Q's into usttauex bins
    Q_usttauex_bin_values{i} = cell(N_usttauex_bins_Site,1);
    Q_usttauex_bin_avg{i} = zeros(N_usttauex_bins_Site,1)*NaN;
    Q_usttauex_bin_std{i} = zeros(N_usttauex_bins_Site,1)*NaN;
    Q_usttauex_bin_SE{i} = zeros(N_usttauex_bins_Site,1)*NaN;

    %get values, avgs, and standard deviations for usttauex bins
    for j = 1:N_usttauex_bins_Site
        %get indices
        bin_ind = find(usttauex_all{i}>=usttauex_bin_min_Site(j)&usttauex_all{i}<=usttauex_bin_max_Site(j));

        %get usttauex in bin
        usttauex_bin_values{i}{j} = usttauex_all{i}(bin_ind);
        N_usttauex_bin_values{i}(j) = length(usttauex_bin_values{i}{j});
        usttauex_bin_avg{i}(j) = mean(usttauex_bin_values{i}{j});
        usttauex_bin_std{i}(j) = std(usttauex_bin_values{i}{j});
        usttauex_bin_SE{i}(j) = usttauex_bin_std{i}(j)./sqrt(N_usttauex_bin_values{i}(j));
        usttauex_bin_min{i}(j) = min(usttauex_bin_values{i}{j});
       
        %get Qs for usttauex bins
        Q_usttauex_bin_values{i}{j} = Q_all{i}(bin_ind);
        Q_usttauex_bin_avg{i}(j) = mean(Q_usttauex_bin_values{i}{j});
        Q_usttauex_bin_std{i}(j) = std(Q_usttauex_bin_values{i}{j});
        Q_usttauex_bin_SE{i}(j) = Q_usttauex_bin_std{i}(j)./sqrt(N_usttauex_bin_values{i}(j));
    end
end


%% FITS TO Q UST*TAUEX

%initialize scaling coefficient for Q versus ust*tauex
C_Qusttauex_fit = zeros(N_Sites,1);
Qpred_usttauex = cell(N_Sites,1);
chi2_Qusttauex = zeros(N_Sites,1);
P_Qusttauex = zeros(N_Sites,1);

%perform fits
for i = 1:N_Sites
   
    %compute scaling coefficient for Q versus ust*tauex
    C_Qusttauex_fit(i) = mean(Q_usttauex_bin_avg{i}./usttauex_bin_avg{i});
    
    %predict values of Q versus ust*tauex
    Qpred_usttauex{i} = usttauex_bin_avg{i}*C_Qusttauex_fit(i);
    
    %compute chi2 value based on prediction
    chi2_Qusttauex(i) = sum((1./Q_usttauex_bin_SE{i}.^2).*(Q_usttauex_bin_avg{i}-Qpred_usttauex{i}).^2);
    N_DF = length(usttauex_bin_avg{i})-2; %number of degrees of freedom
    
    P_Qusttauex(i) = chi2cdf(chi2_Qusttauex(i),N_DF);
end


%% PLOT Q VERSUS UST*TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    errorbar(usttauex_bin_avg{i},Q_usttauex_bin_avg{i},Q_usttauex_bin_SE{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
usttauex_fit = [0 0.14];
for i = 1:N_Sites
    plot(usttauex_bin_avg{i}, Qpred_usttauex{i},'Color',Colors_field{i});
end

xlim([0 0.14]);
ylim([0 45]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear velocity - excess stress product, $$u_{*} \tau_{ex}$$ (kg s$$^{-3}$$)}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
legend(Sites,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%draft paper plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'Q_usttauex.png'],'-dpng');



%% save analysis data
save([folder_AnalysisData,'StressFluxWindows_Analysis.mat'],'*all');


% %% DIFFERENT METRICS OF WIND
% i = 3;
% 
% %mean wind velocity
% figure; clf;
% subplot(3,1,1); hold on;
% plot(tauRe_all{i},Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/4);
% xlabel('Shear stress, $$\tau_{Re}$$ (Pa)','Interpreter','Latex');
% ylabel('Saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'FontSize',PlotFont);
% 
% subplot(3,1,2); hold on;
% plot(ubar_all{i},Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/4);
% xlabel('Mean wind velocity, $$\overline{u}$$ (m/s)','Interpreter','Latex');
% ylabel('Saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'FontSize',PlotFont);
% 
% subplot(3,1,3); hold on;
% plot(abs(uw_all{i}),Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/4);
% xlabel('Momentum flux, $$\overline{uw^{\prime}} \textrm{ (m}^2/\textrm{s}^2)$$ ','Interpreter','Latex');
% ylabel('Saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'FontSize',PlotFont);
% 
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
% print([folder_Plots,'Q_wind_raw.png'],'-dpng');
% 
%
% %% FOR EACH SITE, PERFORM BINNING BY fQ BINS
% 
% %separate tauth based on TFEM into frequency bins
% tauth_fQ_altbin_values = cell(N_Sites,1);
% tauth_fQ_altbin_avg = cell(N_Sites,1);
% tauth_fQ_altbin_std = cell(N_Sites,1);
% tauth_fQ_altbin_SE = cell(N_Sites,1);
% 
% for i = 1:N_Sites
%     tauth_fQ_altbin_values{i} = cell(N_fQ_bins,1);
%     tauth_fQ_altbin_avg{i} = zeros(N_fQ_bins,1)*NaN;
%     tauth_fQ_altbin_std{i} = zeros(N_fQ_bins,1)*NaN;
%     tauth_fQ_altbin_SE{i} = zeros(N_fQ_bins,1)*NaN;
%     
%     %get tauth's for fQ bins
%     for j = 1:N_fQ_bins
%         %get indices for bins
%         bin_ind = find(fQ1_all{i}>=fQ_bins_min(j)&fQ1_all{i}<=fQ_bins_max(j)); %use 1 second frequencies for binning
%                 
%         %tau_th
%         tauth_fQ_altbin_values{i}{j} = tauth_TFEM_all{i}(bin_ind);
%         tauth_fQ_altbin_values{i}{j} = tauth_fQ_altbin_values{i}{j}(~isnan(tauth_fQ_altbin_values{i}{j}));
%         tauth_fQ_altbin_avg{i}(j) = mean(tauth_fQ_altbin_values{i}{j});
%         tauth_fQ_altbin_std{i}(j) = std(tauth_fQ_altbin_values{i}{j});
%         tauth_fQ_altbin_SE{i}(j) = std(tauth_fQ_altbin_values{i}{j})/sqrt(length(tauth_fQ_altbin_values{i}{j}));
%     end 
% end
% 
%
% %% PERFORM FITS FOR FQ
% %set fQ bins
% fQ_bins_min = 0:0.05:0.95;
% fQ_bins_max = 0.05:0.05:1;
% fQ_bins_mid = mean([fQ_bins_min; fQ_bins_max]);
% N_fQ_bins = length(fQ_bins_mid);
% fQ_TFEM_fit_min = [0.65,0.7,0.1]; %mininum for fitting to obtain impact/fluid thresholds
% fQ_TFEM_fit_max = [0.9,0.85,0.9]; %maximum for fitting to obtain impact/fluid thresholds
% 
% %fit to TFEM thresholds
% fQ_TFEM_ind = cell(N_Sites,1);
% tauft_TFEM = zeros(N_Sites,1);
% tauit_TFEM = zeros(N_Sites,1);
% tauth_fit_TFEM = cell(N_Sites,1);
% sigma_tauth_fit_TFEM = cell(N_Sites,1);
% 
% for i = 1:N_Sites
%     fQ_TFEM_ind{i} = intersect(intersect(find(fQ_bins_mid>fQ_TFEM_fit_min(i)),find(fQ_bins_mid<fQ_TFEM_fit_max(i))),find(~isnan(tauth_fQ_altbin_avg{i})));
%     if i==3
%         [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_bins_mid(fQ_TFEM_ind{i})', tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}), fQ_bins_max(fQ_TFEM_ind{i})'-fQ_bins_min(fQ_TFEM_ind{i})', tauth_fQ_altbin_std{i}(fQ_TFEM_ind{i}));
%     else %don't use confidence intervals for Jeri and RG, because there is insufficient data for these
%         [a, b, ~, ~, tauth_fit, sigma_tauth_fit] = linearfit(fQ_bins_mid(fQ_TFEM_ind{i})', tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}));
%     end
%     tauft_TFEM(i) = a;
%     tauit_TFEM(i) = a+b;
%     tauth_fit_TFEM{i} = tauth_fit;
%     sigma_tauth_fit_TFEM{i} = sigma_tauth_fit;
% end 
% 
% %% PLOT TAUTH VERSUS FQ
% figure; clf; hold on; %initialize plot
% %legend_items = cell(N_Sites*2,1);
% for i = 1:N_Sites
%     errorbar(fQ_bins_mid(fQ_TFEM_ind{i}),tauth_fQ_altbin_avg{i}(fQ_TFEM_ind{i}),tauth_fQ_altbin_SE{i}(fQ_TFEM_ind{i}),Markers{i},'MarkerSize',10);
%     %legend_items{2*i-1} = Sites{i}; %add to list of legend items
%     %legend_items{2*i} = 'fit'; %add to list of legend items
% end
% for i = 1:N_Sites
%     plot([0 1],[tauft_TFEM(i) tauit_TFEM(i)],LineColors{i});
% end
% xlabel('frequency of transport');
% ylabel('TFEM inferred threshold stress, \tau_{th}');
% legend(Sites,'Location','NorthEast');
% set(gca,'FontSize',PlotFont);
% set(gcf, 'PaperPosition',[0 0 10 6]);
% print([folder_Plots,'tauth_fQ_all.png'],'-dpng');