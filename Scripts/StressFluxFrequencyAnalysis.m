%% ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/FluxFrequency/'; %folder for plots

%load flux stress window data
load(strcat(folder_ProcessedData,'StressFluxWindows'));

%get info about number of sites
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%PERFORM BINNING, COMPUTE STANDARD ERRORS
%create Q_bins
Q_bins_min = 0:2.5:55;
Q_bins_max = 2.5:2.5:57.5;
Q_bins_mid = mean([Q_bins_min;Q_bins_max]);
N_Q_bins = length(Q_bins_mid);
fQ_Q_bin_avg = cell(N_Sites,1);
fQ_Q_bin_SE = cell(N_Sites,1);

%create ust_bins
ust_bins_min = 0.025:.025:0.55;
ust_bins_max = 0.05:.025:0.575;
ust_bins_mid = mean([ust_bins_min;ust_bins_max]);
N_ust_bins = length(ust_bins_mid);
fQ_ust_bin_avg = cell(N_Sites,1);
fQ_ust_bin_SE = cell(N_Sites,1);

for i=1:N_Sites
    for j=1:N_Q_bins
        bin_ind = find(Q_list{i}>=Q_bins_min(j)&Q_list{i}<=Q_bins_max(j));
        fQ_Q_bin_avg{i}(j) = mean(fQ_list{i}(bin_ind));
        fQ_Q_bin_SE{i}(j) = std(fQ_list{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_ust_bins
        bin_ind = find(ustRe_cal_list{i}>=ust_bins_min(j)&ustRe_cal_list{i}<=ust_bins_max(j));
        fQ_ust_bin_avg{i}(j) = mean(fQ_list{i}(bin_ind));
        fQ_ust_bin_SE{i}(j) = std(fQ_list{i}(bin_ind))/sqrt(length(bin_ind));
    end
end

%CALCULATE EXPECTED fQ BASED ON RANGE OF U*INST
ust_thr = 0.32; %threshold, m/s
fQ_pred = zeros(N_ust_bins,1);
for j=1:N_ust_bins
    mu_ust_inst = ust_bins_mid(j);
    sigma_ust_inst = ust_bins_mid(j)/20;
    fQ_pred(j) = 1-normcdf(ust_thr,mu_ust_inst,sigma_ust_inst);
end

%PLOT THINGS
figure(1); clf; hold on;
for i=1:N_Sites
    errorbar(Q_bins_mid,fQ_Q_bin_avg{i},fQ_Q_bin_SE{i},Markers{i},'MarkerSize',5);
    %plot(Q_list{i},fQ_list{i},Markers{i});
end
xlim([0 max(cell2mat(Q_list))]);
xlabel('Q (g m^{-1} s^{-1})');
ylabel('f_Q');
set(gca,'FontSize',16);
legend_values = Sites;
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FrequencyFlux.png'],'-dpng');

figure(2); clf; hold on;
for i=1:N_Sites
    errorbar(ust_bins_mid,fQ_ust_bin_avg{i},fQ_ust_bin_SE{i},Markers{i},'MarkerSize',5);
    %plot(ustRe_cal_list{i},fQ_list{i},Markers{i});
end
plot(ust_bins_mid,fQ_pred,'k','LineWidth',2);
xlim([0 max(cell2mat(ustRe_cal_list))]);
xlabel('u_{*} (m/s)');
ylabel('f_Q');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = 'pred: u_{*,th}=0.32, \sigma_{u,*,inst}=u_{*}/20';
h_legend = legend(legend_values,'Location','NorthWest');
%set(h_legend,'FontSize',16);
print([folder_Plots,'FrequencyUst.png'],'-dpng');