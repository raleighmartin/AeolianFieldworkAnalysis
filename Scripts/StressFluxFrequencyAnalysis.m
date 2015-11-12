%% ANALYZE HOW SALTATION PROPERTIES VARY WITH TRANSPORT FREQUENCY

%initialize
clearvars;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_AnalysisData = '../AnalysisData/'; %folder for storing data output
folder_Plots = '../PlotOutput/FluxFrequency/'; %folder for plots

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows_all'));
load(strcat(folder_AnalysisData,'StressFluxWindows_continuous'));
load(strcat(folder_AnalysisData,'StressFluxWindows_intermittent'));

%get info about number of sites
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%PERFORM BINNING, COMPUTE AVERAGES AND STANDARD DEVIATIONS FOR BINS
%create Q_bins
Q_bins_min_all = [0, 0:2.5:55];
Q_bins_max_all = 0:2.5:57.5;
Q_bins_mid_all = mean([Q_bins_min_all;Q_bins_max_all]);
N_Q_bins_all = length(Q_bins_mid_all);
Q_bins_min = cell(N_Sites,1);
Q_bins_max = cell(N_Sites,1);
Q_bins_mid = cell(N_Sites,1);
N_Q_bins = zeros(N_Sites,1);
fQ_Q_bin_values = cell(N_Sites,1);
fQ_Q_bin_avg = cell(N_Sites,1);
fQ_Q_bin_std = cell(N_Sites,1);
fQ_Q_bin_25 = cell(N_Sites,1);
fQ_Q_bin_50 = cell(N_Sites,1);
fQ_Q_bin_75 = cell(N_Sites,1);

%create ust_bins
ust_bins_min_all = 0.025:.025:0.55;
ust_bins_max_all = 0.05:.025:0.575;
ust_bins_mid_all = mean([ust_bins_min_all;ust_bins_max_all]);
N_ust_bins_all = length(ust_bins_mid_all);
ust_bins_min = cell(N_Sites,1);
ust_bins_max = cell(N_Sites,1);
ust_bins_mid = cell(N_Sites,1);
N_ust_bins = zeros(N_Sites,1);
fQ_ust_bin_values = cell(N_Sites,1);
fQ_ust_bin_avg = cell(N_Sites,1);
fQ_ust_bin_std = cell(N_Sites,1);
fQ_ust_bin_25 = cell(N_Sites,1);
fQ_ust_bin_50 = cell(N_Sites,1);
fQ_ust_bin_75 = cell(N_Sites,1);

%apply values to bins
for i=1:N_Sites
    
    %flux bins
    fQ_Q_bin_values{i} = cell(N_Q_bins_all,1);
    fQ_Q_bin_avg{i} = zeros(N_Q_bins_all,1)*NaN;
    fQ_Q_bin_std{i} = zeros(N_Q_bins_all,1)*NaN;
    fQ_Q_bin_25{i} = zeros(N_Q_bins_all,1)*NaN;
    fQ_Q_bin_50{i} = zeros(N_Q_bins_all,1)*NaN;
    fQ_Q_bin_75{i} = zeros(N_Q_bins_all,1)*NaN;
    for j=1:N_Q_bins_all
        bin_ind = find(Q_all{i}>Q_bins_min_all(j)&Q_all{i}<=Q_bins_max_all(j));
        if ~isempty(bin_ind)
            fQ_Q_bin_values{i}{j} = fQ_norm_all{i}(bin_ind);
            %fQ_Q_bin_values{i}{j} = fQ_all{i}(bin_ind);
            fQ_Q_bin_avg{i}(j) = mean(fQ_Q_bin_values{i}{j});
            fQ_Q_bin_std{i}(j) = std(fQ_Q_bin_values{i}{j});
            fQ_Q_bin_25{i}(j) = prctile(fQ_Q_bin_values{i}{j},25);
            fQ_Q_bin_50{i}(j) = prctile(fQ_Q_bin_values{i}{j},50);
            fQ_Q_bin_75{i}(j) = prctile(fQ_Q_bin_values{i}{j},75);
        end
    end
    %keep only intervals that have something in them
    notempty_ind = find(~isnan(fQ_Q_bin_avg{i}));
    Q_bins_min{i} = Q_bins_min_all(notempty_ind);
    Q_bins_max{i} = Q_bins_max_all(notempty_ind);
    Q_bins_mid{i} = Q_bins_mid_all(notempty_ind);
    fQ_Q_bin_values{i} = fQ_Q_bin_values{i}(notempty_ind);
    fQ_Q_bin_avg{i} = fQ_Q_bin_avg{i}(notempty_ind);
    fQ_Q_bin_std{i} = fQ_Q_bin_std{i}(notempty_ind);
    fQ_Q_bin_25{i} = fQ_Q_bin_25{i}(notempty_ind);
    fQ_Q_bin_50{i} = fQ_Q_bin_50{i}(notempty_ind);
    fQ_Q_bin_75{i} = fQ_Q_bin_75{i}(notempty_ind);
    N_Q_bins(i) = length(notempty_ind);
    
    %shear velocity bins
    fQ_ust_bin_values{i} = cell(N_ust_bins_all,1);
    fQ_ust_bin_avg{i} = zeros(N_ust_bins_all,1)*NaN;
    fQ_ust_bin_std{i} = zeros(N_ust_bins_all,1)*NaN;
    fQ_ust_bin_25{i} = zeros(N_ust_bins_all,1)*NaN;
    fQ_ust_bin_50{i} = zeros(N_ust_bins_all,1)*NaN;
    fQ_ust_bin_75{i} = zeros(N_ust_bins_all,1)*NaN;
    for j=1:N_ust_bins_all
        bin_ind = find(ustRe_cal_all{i}>=ust_bins_min_all(j)&ustRe_cal_all{i}<=ust_bins_max_all(j));
        if ~isempty(bin_ind)
            fQ_ust_bin_values{i}{j} = fQ_norm_all{i}(bin_ind);
            %fQ_ust_bin_values{i}{j} = fQ_all{i}(bin_ind);
            fQ_ust_bin_avg{i}(j) = mean(fQ_ust_bin_values{i}{j});
            fQ_ust_bin_std{i}(j) = std(fQ_ust_bin_values{i}{j});
            fQ_ust_bin_25{i}(j) = prctile(fQ_ust_bin_values{i}{j},25);
            fQ_ust_bin_50{i}(j) = prctile(fQ_ust_bin_values{i}{j},50);
            fQ_ust_bin_75{i}(j) = prctile(fQ_ust_bin_values{i}{j},75);
        end
    end
    %keep only intervals that have something in them
    notempty_ind = find(~isnan(fQ_ust_bin_avg{i}));
    ust_bins_min{i} = ust_bins_min_all(notempty_ind);
    ust_bins_max{i} = ust_bins_max_all(notempty_ind);
    ust_bins_mid{i} = ust_bins_mid_all(notempty_ind);
    fQ_ust_bin_values{i} = fQ_ust_bin_values{i}(notempty_ind);
    fQ_ust_bin_avg{i} = fQ_ust_bin_avg{i}(notempty_ind);
    fQ_ust_bin_std{i} = fQ_ust_bin_std{i}(notempty_ind);
    fQ_ust_bin_25{i} = fQ_ust_bin_25{i}(notempty_ind);
    fQ_ust_bin_50{i} = fQ_ust_bin_50{i}(notempty_ind);
    fQ_ust_bin_75{i} = fQ_ust_bin_75{i}(notempty_ind);
    N_ust_bins(i) = length(notempty_ind);
end

%calculate expected fQ based on range of u*inst
ust_thr = 0.32; %threshold, m/s
fQ_pred = zeros(N_ust_bins_all,1);
for j=1:N_ust_bins_all
    mu_ust_inst = ust_bins_mid_all(j);
    sigma_ust_inst = ust_bins_mid_all(j)/20;
    fQ_pred(j) = 1-normcdf(ust_thr,mu_ust_inst,sigma_ust_inst);
end

%PLOT THINGS FOR ALL SITES
%frequency versus flux
figure(1); clf; hold on;
for i=1:N_Sites
    errorbar(Q_bins_mid{i},fQ_Q_bin_50{i},fQ_Q_bin_50{i}-fQ_Q_bin_25{i},fQ_Q_bin_75{i}-fQ_Q_bin_50{i},Markers{i},'MarkerSize',5);
    %errorbar(Q_bins_mid{i},fQ_Q_bin_avg{i},fQ_Q_bin_std{i},Markers{i},'MarkerSize',5);
    %plot(Q_all{i},fQ_all{i},Markers{i});
end
xlim([0 max(cell2mat(Q_all))]);
xlabel('Q (g m^{-1} s^{-1})');
ylabel('f_Q');
set(gca,'FontSize',16);
legend_values = Sites;
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'FrequencyFlux.png'],'-dpng');

%frequency versus u*
figure(2); clf; hold on;
for i=1:N_Sites
    errorbar(ust_bins_mid{i},fQ_ust_bin_50{i},fQ_ust_bin_50{i}-fQ_ust_bin_25{i},fQ_ust_bin_75{i}-fQ_ust_bin_50{i},Markers{i},'MarkerSize',5);
    %errorbar(ust_bins_mid{i},fQ_ust_bin_avg{i},fQ_ust_bin_std{i},Markers{i},'MarkerSize',5);
    %plot(ustRe_cal_all{i},fQ_all{i},Markers{i});
end
plot(ust_bins_mid_all,fQ_pred,'k','LineWidth',2);
xlim([0 max(cell2mat(ustRe_cal_all))]);
xlabel('u_{*} (m/s)');
ylabel('f_Q');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = 'pred: u_{*,th}=0.32, \sigma_{u,*,inst}=u_{*}/20';
h_legend = legend(legend_values,'Location','NorthWest');
%set(h_legend,'FontSize',16);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
print([folder_Plots,'FrequencyUst.png'],'-dpng');

%OCEANO ONLY ANALYSIS
%look at histograms of flux for certain u* ranges
i = 3; %examine only Oceano data
ust_bins_min_Q_hist = [0 0.25 0.3 0.35 0.4]; %lower u*s for histogram bins
ust_bins_max_Q_hist = [0.25 0.3 0.35 0.4 0.52]; %upper u*s for histogram bins
ust_bins_mid_Q_hist = mean([ust_bins_min_Q_hist;ust_bins_max_Q_hist]); %midpoint u*s for histogram bins
N_ust_bins_Q_hist = length(ust_bins_mid_Q_hist); %number of u* bins
Q_hist_ust_bin_continuous = zeros(N_ust_bins_Q_hist,N_Q_bins(i));
Q_hist_ust_bin_intermittent = zeros(N_ust_bins_Q_hist,N_Q_bins(i));

%go through each u* bin
for j = 1:N_ust_bins_Q_hist
    %for this u* bin, get all the fluxes associated with continuous transport
    bin_ind_continuous = ...
        find(ustRe_cal_continuous{i}>=ust_bins_min_Q_hist(j)&ustRe_cal_continuous{i}<=ust_bins_max_Q_hist(j));
    Q_bin_continuous = Q_continuous{i}(bin_ind_continuous);

    %for this u* bin, get all the fluxes associated with intermittent transport
    bin_ind_intermittent = ...
        find(ustRe_cal_intermittent{i}>=ust_bins_min_Q_hist(j)&ustRe_cal_intermittent{i}<=ust_bins_max_Q_hist(j));
    Q_bin_intermittent = Q_intermittent{i}(bin_ind_intermittent);

    %go through each Q bin
    for k = 1:N_Q_bins(i)
        %for this Q bin, get the number of continuous intervals
        Q_hist_ust_bin_continuous(j,k) = length(find(...
            Q_bin_continuous>=Q_bins_min{i}(k)&Q_bin_continuous<=Q_bins_max{i}(k)));
        
        %for this Q bin, get the number of intermittent intervals
        Q_hist_ust_bin_intermittent(j,k) = length(find(...
            Q_bin_intermittent>=Q_bins_min{i}(k)&Q_bin_intermittent<=Q_bins_max{i}(k)));
    end
end

%flux histogram for each u* bin
figure(3);
for j = 1:N_ust_bins_Q_hist
    subplot(N_ust_bins_Q_hist,1,j);
    Q_hist_barplot = [Q_hist_ust_bin_intermittent(j,:)' Q_hist_ust_bin_continuous(j,:)'];
    bar(Q_bins_mid{3},Q_hist_barplot,'stacked');
    xlim([0 35]);
    ylabel('N');
    title(['u_{*} = ',num2str(ust_bins_min_Q_hist(j)),'-',num2str(ust_bins_max_Q_hist(j)),' m/s']);
    set(gca,'FontSize',16);
end
xlabel('Q (g/m/s)','FontSize',16);
subplot(N_ust_bins_Q_hist,1,1);
h_legend = legend('Intermittent','Continuous');
set(h_legend,'FontSize',16);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 14])
print([folder_Plots,'FluxFrequencyHistograms.png'],'-dpng')