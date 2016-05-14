%% ANALYZE SALTATION HEIGHT RELATIONSHIPS

%initialize
clearvars;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots

%load data
load(strcat(folder_ProcessedData,'StressFluxContinuousWindows')); %stress/flux continuous windows
load(strcat(folder_ProcessedData,'C_ustThr_sqrtd50')); %threshold scaling
load(strcat(folder_ProcessedData,'GreeleyNamikasData')); %literature data

%get info about sites, set markers for plotting
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%set parameter values
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

%% MAKE CALCULATIONS FOR LITERATURE VALUES
ustThr_Greeley96 = C_ustThr_sqrtd50*sqrt(d50_Greeley96); %get Greeley ustThr value
ustThr_Namikas03 = C_ustThr_sqrtd50*sqrt(d50_Namikas03); %get Namikas ustThr value
ust_ratio_Greeley96 = ust_Greeley96/ustThr_Greeley96; %get u* ratio values for Greeley
ust_ratio_Namikas03 = ust_Namikas03/ustThr_Namikas03; %get u* ratio values for Namikas
ust_ex_Greeley96 = ust_Greeley96-ustThr_Greeley96; %get u*ex values for Greeley
ust_ex_Namikas03 = ust_Namikas03-ustThr_Namikas03; %get u*ex values for Namikas
ze_mean_Greeley96 = mean(zbar_Greeley96); %mean saltation height for Greeley
ze_SE_Greeley96 = std(zbar_Greeley96)/sqrt(N_Greeley96); %standard error for mean saltation height for Greeley
ze_mean_Namikas03 = mean(zbar_Namikas03); %mean saltation height for Namikas
ze_SE_Namikas03 = std(zbar_Namikas03)/sqrt(N_Namikas03); %standard error for mean saltation height Namikas
ze_ratio_Greeley96 = 1e3*zbar_Greeley96/d50_Greeley96; %get ze/d50 values for Greeley
ze_ratio_Namikas03 = 1e3*zbar_Namikas03/d50_Namikas03; %get ze/d50 values for Greeley

%% SEPARATE SALTATION HEIGHT VALUES INTO U* BINS
%generate u* bins
ust_bins_min = 0.15:0.05:0.55;
ust_bins_max = 0.2:0.05:0.6;
ust_bins_mid = mean([ust_bins_min; ust_bins_max]);
N_bins = length(ust_bins_mid);

%initialize cell arrays of saltation heights in ust bins for all sites
ze_ust_bins_mean = cell(N_Sites,1);
ze_ust_bins_SE = cell(N_Sites,1);
for i=1:N_Sites
    
    %initialize lists of ust binned saltation heights for each site
    ze_ust_bins_mean{i} = zeros(N_bins,1);
    ze_ust_bins_SE{i} = zeros(N_bins,1);
    
    %go through each u* bin, get saltation heights for this bin
    for j = 1:N_bins
        ust_bin_ind = find(ustRe_cal_continuous{i}>=ust_bins_min(j)&ustRe_cal_continuous{i}<=ust_bins_max(j));
        ze_d50_bin_all = ze_continuous{i}(ust_bin_ind);
        not_outlier_ind = find(ze_d50_bin_all>0&ze_d50_bin_all<0.30);
        ze_ust_bins_mean{i}(j) = mean(ze_d50_bin_all(not_outlier_ind));
        ze_ust_bins_SE{i}(j) = std(ze_d50_bin_all(not_outlier_ind))./sqrt(length(not_outlier_ind));
    end
end

%% SEPARATE SALTATION HEIGHT RATIO VALUES INTO U* RATIO BINS
%generate ust_ratio bins
ust_ratio_bins_min = 0.5:0.25:1.75;
ust_ratio_bins_max = 0.75:0.25:2;
ust_ratio_bins_mid = mean([ust_ratio_bins_min; ust_ratio_bins_max]);
N_bins = length(ust_ratio_bins_mid);

%initialize cell arrays of ze/d50 ratios and ze/d90 versus ustratio bin all sites
ze_d50_ustratio_bins_mean = cell(N_Sites,1);
ze_d50_ustratio_bins_SE = cell(N_Sites,1);
ze_d90_ustratio_bins_mean = cell(N_Sites,1);
ze_d90_ustratio_bins_SE = cell(N_Sites,1);
for i = 1:N_Sites; %go through each site
    %initialize lists of saltation height ratios for each bin
    ze_d50_ustratio_bins_mean{i} = zeros(N_bins,1);
    ze_d50_ustratio_bins_SE{i} = zeros(N_bins,1);
    ze_d90_ustratio_bins_mean{i} = zeros(N_bins,1);
    ze_d90_ustratio_bins_SE{i} = zeros(N_bins,1);
    
    %compute ustThr, ust_ratio, ze/d50, and ze/d90 for all bins
    ustThr_Site = C_ustThr_sqrtd50*sqrt(d50_continuous{i});
    ust_ratio = ustRe_cal_continuous{i}./ustThr_Site;
    ze_d50_ratio = 1000*ze_continuous{i}./d50_continuous{i};
    ze_d90_ratio = 1000*ze_continuous{i}./d90_continuous{i};
    
    %separate data into individual bins
    for j = 1:N_bins
        ust_bin_ind = find(ust_ratio>=ust_ratio_bins_min(j)&ust_ratio<=ust_ratio_bins_max(j));
        ze_d50_bin_all = ze_d50_ratio(ust_bin_ind);
        ze_d90_bin_all = ze_d90_ratio(ust_bin_ind);
        not_outlier_ind = find(ze_d50_bin_all>0&ze_d50_bin_all<500);
        ze_d50_ustratio_bins_mean{i}(j) = mean(ze_d50_bin_all(not_outlier_ind));
        ze_d50_ustratio_bins_SE{i}(j) = std(ze_d50_bin_all(not_outlier_ind))/sqrt(length(not_outlier_ind));
        ze_d90_ustratio_bins_mean{i}(j) = mean(ze_d90_bin_all(not_outlier_ind));
        ze_d90_ustratio_bins_SE{i}(j) = std(ze_d90_bin_all(not_outlier_ind))/sqrt(length(not_outlier_ind));
    end
end

%% SEPARATE SALTATION HEIGHT VALUES INTO U*EX BINS
%generate ustex bins
ustex_bins_min = -0.05:0.025:0.225;
ustex_bins_max = -0.025:0.025:0.25;
ustex_bins_mid = mean([ustex_bins_min; ustex_bins_max]);
N_bins = length(ustex_bins_mid);

%initialize cell arrays of ze/d50 for u*ex for all sites
ze_d50_ustex_bins_mean = cell(N_Sites,1);
ze_d50_ustex_bins_SE = cell(N_Sites,1);
for i = 1:N_Sites;
    %initialize lists of saltation heights for u*ex for all bins
    ze_d50_ustex_bins_mean{i} = zeros(N_bins,1);
    ze_d50_ustex_bins_SE{i} = zeros(N_bins,1);
    ustThr_Site = C_ustThr_sqrtd50*sqrt(d50_continuous{i});
    ust_ustex = ustRe_cal_continuous{i}-ustThr_Site;
    ze_d50_ustex = 1000*ze_continuous{i}./d50_continuous{i};
    for j = 1:N_bins
        ust_bin_ind = find(ust_ustex>=ustex_bins_min(j)&ust_ustex<=ustex_bins_max(j));
        ze_d50_bin_all = ze_d50_ustex(ust_bin_ind);
        not_outlier_ind = find(ze_d50_bin_all>0&ze_d50_bin_all<500);
        ze_d50_ustex_bins_mean{i}(j) = mean(ze_d50_bin_all(not_outlier_ind));
        ze_d50_ustex_bins_SE{i}(j) = std(ze_d50_bin_all(not_outlier_ind))/sqrt(length(not_outlier_ind));
    end
end

%% SEPARATE SALTATION HEIGHT VALUES INTO D50 BINS
%initialize cell arrays of saltation heights and d50s for all sites
d50_bins = cell(N_Sites,1);
ze_d50_bins_mean = cell(N_Sites,1);
ze_d50_bins_SE = cell(N_Sites,1);
for i=1:N_Sites
    d50_bins{i} = unique(d50_continuous{i});
    N_bins = length(d50_bins{i});
    ze_d50_bins_mean{i} = zeros(N_bins,1);
    ze_d50_bins_SE{i} = zeros(N_bins,1);
    Q_positive_ind = find(Q_continuous{i}>0);
    for j = 1:N_bins
        d50_bin_ind = find(d50_continuous{i}==d50_bins{i}(j));
        ze_d50_bin_all = ze_continuous{i}(d50_bin_ind);
        not_outlier_ind = find(ze_d50_bin_all>0&ze_d50_bin_all<0.30);
        ze_d50_bins_mean{i}(j) = mean(ze_d50_bin_all(not_outlier_ind));
        ze_d50_bins_SE{i}(j) = std(ze_d50_bin_all(not_outlier_ind))./sqrt(length(not_outlier_ind));
    end
end

%% PLOT RESULTS
%plot z_salt data binned by u* range
figure(1); clf; hold on;
for i=1:N_Sites
    errorbar(ust_bins_mid,ze_ust_bins_mean{i},ze_ust_bins_SE{i},Markers{i});
end
plot(ust_Greeley96,zbar_Greeley96,'^k'); %add Greeley data
plot(ust_Namikas03,zbar_Namikas03,'dk'); %add Namikas data
xlabel('u_{*,Re} (m/s)','FontSize',16);
ylabel('z_{salt} (m)','FontSize',16);
xlim([0 0.65]);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthOutside');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightUst_Wenglor.png'],'-dpng');

%plot saltation height / d50_surface versus u*/u*th
figure(2); clf; hold on;
for i = 1:N_Sites;
    errorbar(ust_ratio_bins_mid,ze_d50_ustratio_bins_mean{i},ze_d50_ustratio_bins_SE{i},Markers{i});
end
plot(ust_ratio_Greeley96,ze_ratio_Greeley96,'^k'); %add Greeley data
plot(ust_ratio_Namikas03,ze_ratio_Namikas03,'dk'); %add Namikas data
xlabel('u_{*,Re}/u_{*,thr}','FontSize',16);
ylabel('z_{salt}/d_{50}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthOutside');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightzOverD50_Wenglor.png'],'-dpng');

%plot saltation height / d90_surface versus u*/u*th
figure(3); clf; hold on;
for i = 1:N_Sites;
    errorbar(ust_ratio_bins_mid,ze_d90_ustratio_bins_mean{i},ze_d90_ustratio_bins_SE{i},Markers{i});
end
xlabel('u_{*,Re}/u_{*,thr}','FontSize',16);
ylabel('z_{salt}/d_{90}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
h_legend = legend(legend_values,'Location','NorthOutside');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightzOverD90_Wenglor.png'],'-dpng');

%plot saltation height / d50_surface versus u*_ex
figure(4); clf; hold on;
for i = 1:N_Sites;
    errorbar(ustex_bins_mid,ze_d50_ustex_bins_mean{i},ze_d50_ustex_bins_SE{i},Markers{i},'MarkerSize',10);
end
plot(ust_ex_Greeley96,ze_ratio_Greeley96,'^k','MarkerSize',10); %plot Greeley values
plot(ust_ex_Namikas03,ze_ratio_Namikas03,'dk','MarkerSize',10); %plot Namikas values
xlabel('u_{*,Re}-u_{*,thr} (m/s)','FontSize',16);
ylabel('z_{salt}/d_{50}','FontSize',16);
%ylim([130 280]);
%xlim([-0.05 0.35]);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley et al. (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthOutside');
set(h_legend,'FontSize',16);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 5])
print([folder_Plots,'FluxHeightzOverD_ustEx_Wenglor.png'],'-dpng');

%plot z_salt data binned by d50 range
figure(5); clf; hold on;
for i=1:N_Sites
    errorbar(d50_bins{i},ze_d50_bins_mean{i},ze_d50_bins_SE{i},Markers{i});
end
errorbar(d50_Greeley96,ze_mean_Greeley96,ze_SE_Greeley96,'^k'); %add Greeley data
errorbar(d50_Namikas03,ze_mean_Namikas03,ze_SE_Namikas03,'dk'); %add Namikas data
xlabel('d_{50} (mm)','FontSize',16);
ylabel('z_{salt} (m)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['Greeley (1996)'];
legend_values{length(Sites)+2} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FluxHeightD50_Wenglor.png'],'-dpng');