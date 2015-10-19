%% LOOK AT VALUES OF WENGLOR CALIBRATION VERSUS INSTRUMENT HEIGHT

%initialize
clearvars;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/WenglorCalibration/'; %folder for plots

%load flux stress window data
load(strcat(folder_ProcessedData,'StressFluxWindows'));

%get info about number of sites, set markers
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%% GET RANGE OF CALIBRATION VALUES FOR Z BINS
%generate z bins
z_bins_min = 0:0.05:0.45;
z_bins_max = 0.05:0.05:0.5;
z_bins_mid = mean([z_bins_min; z_bins_max]);
N_bins = length(z_bins_mid);

%initialize lists of qcal binned for all sites
qcal_z_bins_mean = cell(N_Sites,1);
qcal_z_bins_std = cell(N_Sites,1);
qcal_z_bins_lower = cell(N_Sites,1);
qcal_z_bins_upper = cell(N_Sites,1);

%initialize lists of all calibration factors and heights by site
zq_all = cell(N_Sites,1);
qcal_all = cell(N_Sites,1);

for i=1:N_Sites
    %initialize lists of qcal binned within each site
    qcal_z_bins_mean{i} = zeros(N_bins,1)*NaN;
    qcal_z_bins_std{i} = zeros(N_bins,1)*NaN;
    qcal_z_bins_lower{i} = zeros(N_bins,1)*NaN;
    qcal_z_bins_upper{i} = zeros(N_bins,1)*NaN;
    
    %reshape z and qcal into lists for site
    zq_all{i} = cell2mat(zq_list{i}')';
    qcal_all{i} = cell2mat(qcal_list{i}')';
    
    %go through each z bin, get qcal for this bin
    for j = 1:N_bins
        z_bin_ind = find(zq_all{i}>=z_bins_min(j)&zq_all{i}<=z_bins_max(j));
        qcal_z_bin_all = qcal_all{i}(z_bin_ind);
        if ~isempty(qcal_z_bin_all)
            qcal_z_bins_mean{i}(j) = mean(qcal_z_bin_all);
            qcal_z_bins_std{i}(j) = std(qcal_z_bin_all);
            qcal_z_bins_lower{i}(j) = min(qcal_z_bin_all);
            qcal_z_bins_upper{i}(j) = max(qcal_z_bin_all);
        end
    end
end

%% plot all calibration factors versus height
figure(1); clf; hold on;
for i=1:N_Sites
    plot(zq_all{i},qcal_all{i},Markers{i});
end
set(gca,'yscale','log');
ylim([min(cell2mat(qcal_all)) max(cell2mat(qcal_all))])
xlabel('z_{Wenglor} (m)');
ylabel('C_{q/n} (g m^{-2} s^{-1} count^{-1})');
h_legend = legend(Sites,'Location','NorthOutside');
set(h_legend,'FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'WenglorCalibration_z.png'],'-dpng');

%% plot binned calibration factors versus height bin
figure(2); clf; hold on;
for i=1:N_Sites
    %errorbar(z_bins_mid,qcal_z_bins_mean{i},qcal_z_bins_std{i},Markers{i});
    errorbar(z_bins_mid,qcal_z_bins_mean{i},qcal_z_bins_lower{i},qcal_z_bins_upper{i},Markers{i});
end
set(gca,'yscale','log');
ylim([min(cell2mat(qcal_all)) max(cell2mat(qcal_all))])
xlabel('z_{Wenglor} (m)');
ylabel('C_{q/n} (g m^{-2} s^{-1} count^{-1})');
h_legend = legend(Sites,'Location','NorthOutside');
set(h_legend,'FontSize',16);
set(gca,'FontSize',16);
print([folder_Plots,'WenglorCalibration_zBins.png'],'-dpng');