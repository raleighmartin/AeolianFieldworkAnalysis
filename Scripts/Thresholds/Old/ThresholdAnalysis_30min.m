%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% plotting info
PlotFont = 14;
PlotMarkers = {'s','d','o','<','>'};
PlotColors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};

%% site info
tauit_intercept = [0.137, 0.106, 0.086];
sigma_tauit_intercept = [0.015, 0.020, 0.008];

%% info about time scales for core analysis
delta_t_avg_analysis = duration(0,0,1); %delta_t for window-avg analysis
delta_t_lowpass_analysis = duration(0,0,25); %delta_t for lowpass analysis

%% diurnal info
diurnalrange_start_time = [duration(0,0,0); duration(14,0,0); duration(16,0,0)];
diurnalrange_end_time = [duration(14,0,0); duration(16,0,0); duration(24,0,0)];
diurnalrange_labels = {'midday','early afternoon','late afternoon'};
N_diurnalrange = length(diurnalrange_start_time);

%% date info
daterange_start = [datetime(2015,5,15),datetime(2015,5,23),datetime(2015,6,1)];
daterange_end = [datetime(2015,5,19),datetime(2015,5,28),datetime(2015,6,4)];
daterange_labels = {'May 15-19','May 23-28','June 1-4'};
N_daterange = length(daterange_start);

%% binning info
fQ_min = 0.05;
fQ_max = 0.95;
fQ_bin_minrange = 0.1;
fQ_bin_maxrange = 0.2;
fQ_bin_N_min = 3;

%% information about where to load/save data, plots, and functions
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for window data
LoadData_1_Path = strcat(folder_LoadData,'ActivityWindows_30min.mat'); %path for loading window data
LoadData_2_Path = strcat(folder_LoadData,'TimeWindows_30min.mat'); %path for loading time windows
SaveData_Path = strcat(folder_SaveData,'ThresholdWindows_30min.mat'); %path for saving output data

%% load window data
load(LoadData_1_Path);
load(LoadData_2_Path);

%% get info about delta t's
ind_delta_t_avg_analysis = find(delta_t_avg_window==delta_t_avg_analysis); %index of window-averaged delta t for core analysis
N_delta_t_avg = length(delta_t_avg_window);
delta_t_labels = cell(N_delta_t_avg,1); %initialize list of delta_t labels
for j=1:N_delta_t_avg
    delta_t_labels{j} = strcat('{\delta}t = ',num2str(seconds(delta_t_avg_window(j))),' s'); %add to list of delta_t labels
end
ind_delta_t_lowpass_analysis = find(delta_t_lowpass_window==delta_t_lowpass_analysis); %index of window-averaged delta t for core analysis

%% Separate data according to time of day
midpointtime_window = cell(N_Sites,1); %hour corresponding to midpoint time in window
diurnal_ind_window = cell(N_Sites,1); %index corresponding to diurnal range (1 for midday, 2 for early afternoon, 3 for late afternoon)
for i = 1:N_Sites
    midpointtime_window{i} = mean([StartTime_window{i},EndTime_window{i}],2)-Date_window{i}; %get midpoint times for windows
    diurnal_ind_window{i} = zeros(size(midpointtime_window{i})); %intialize list of indices for diurnal range
    for j = 1:N_diurnalrange
        ind_range = find(midpointtime_window{i}>=diurnalrange_start_time(j)&midpointtime_window{i}<diurnalrange_end_time(j)); %get indices of windows in diurnal range
        diurnal_ind_window{i}(ind_range)=j; %assign index for diurnal range
    end
end

%separate fQ and tauth into diurnal range, generate binned values
fQ_diurnal = cell(N_Sites,1);
fQ_diurnal_bin_avg = cell(N_Sites,1);
fQ_diurnal_bin_SE = cell(N_Sites,1);
tauth_diurnal = cell(N_Sites,1);
tauth_diurnal_bin_avg = cell(N_Sites,1);
tauth_diurnal_bin_SE = cell(N_Sites,1);
tauft_diurnal = cell(N_Sites,1);
sigma_tauft_diurnal = cell(N_Sites,1);
tauit_diurnal = cell(N_Sites,1);
sigma_tauit_diurnal = cell(N_Sites,1);
for i = 1:N_Sites
    fQ_diurnal{i} = cell(N_diurnalrange,1);
    fQ_diurnal_bin_avg{i} = cell(N_diurnalrange,1);
    fQ_diurnal_bin_SE{i} = cell(N_diurnalrange,1);
    tauth_diurnal{i} = cell(N_diurnalrange,1);
    tauth_diurnal_bin_avg{i} = cell(N_diurnalrange,1);
    tauth_diurnal_bin_SE{i} = cell(N_diurnalrange,1);
    tauit_diurnal{i} = zeros(N_diurnalrange,1);
    sigma_tauit_diurnal{i} = zeros(N_diurnalrange,1);
    tauft_diurnal{i} = zeros(N_diurnalrange,1);
    sigma_tauft_diurnal{i} = zeros(N_diurnalrange,1);
    for j=1:N_diurnalrange
        %get values
        fQ_diurnal{i}{j} = fQ_avg_window{i}(diurnal_ind_window{i}==j,ind_delta_t_avg_analysis);
        tauth_diurnal{i}{j} = tauth_fQ_avg_window{i}(diurnal_ind_window{i}==j,ind_delta_t_avg_analysis);

        fQ_diurnal_binning = fQ_diurnal{i}{j}(fQ_diurnal{i}{j}>=fQ_min&fQ_diurnal{i}{j}<=fQ_max);
        tauth_diurnal_binning = tauth_diurnal{i}{j}(fQ_diurnal{i}{j}>=fQ_min&fQ_diurnal{i}{j}<=fQ_max);
        
        %perform binning
        if ~isempty(fQ_diurnal_binning)
            [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_diurnal_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
            N_fQ_bins = length(fQ_bin_min);
            fQ_diurnal_bin_avg{i}{j} = fQ_bin_avg;
            fQ_diurnal_bin_SE{i}{j} = fQ_bin_SE;
            tauth_diurnal_bin_avg{i}{j} = zeros(N_fQ_bins,1);
            tauth_diurnal_bin_SE{i}{j} = zeros(N_fQ_bins,1);
            for k=1:N_fQ_bins   
                fQ_diurnal_bin_ind = find(fQ_diurnal_binning>=fQ_bin_min(k)&fQ_diurnal_binning<=fQ_bin_max(k));
                tauth_diurnal_bin_avg{i}{j}(k) = mean(tauth_diurnal_binning(fQ_diurnal_bin_ind));
                tauth_diurnal_bin_SE{i}{j}(k) = std(tauth_diurnal_binning(fQ_diurnal_bin_ind))/sqrt(length(fQ_diurnal_bin_ind));
            end
        end
        
        %get best fit values
        [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = linearfit(fQ_diurnal_bin_avg{i}{j}, tauth_diurnal_bin_avg{i}{j}, tauth_diurnal_bin_SE{i}{j});
        tauft_diurnal{i}(j) = a;
        sigma_tauft_diurnal{i}(j) = sigma_a;
        tauit_diurnal{i}(j) = a+b;
        sigma_tauit_diurnal{i}(j) = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);
    end
end

%plot binned values of fQ and tuath in separate time of day ranges
for i = 1:N_Sites;
    figure(1); clf; hold on; %initialize plot

    %plot average values
    for j = 1:N_diurnalrange
        plot(fQ_diurnal_bin_avg{i}{j},tauth_diurnal_bin_avg{i}{j},PlotMarkers{j},'Color',PlotColors{j});
    end

    %plot SE values
    for j = 1:N_diurnalrange
        N_fQ_bins = length(fQ_diurnal_bin_avg{i}{j});
        for k = 1:N_fQ_bins
            plot(fQ_diurnal_bin_avg{i}{j}(k)+fQ_diurnal_bin_SE{i}{j}(k)*[-1 1],tauth_diurnal_bin_avg{i}{j}(k)*[1 1],'Color',PlotColors{j});
            plot(fQ_diurnal_bin_avg{i}{j}(k)*[1 1],tauth_diurnal_bin_avg{i}{j}(k)+tauth_diurnal_bin_SE{i}{j}(k)*[-1 1],'Color',PlotColors{j});
        end
    end

    %plot fit values
    for j = 1:N_diurnalrange
        plot([0 1],[tauft_diurnal{i}(j) tauit_diurnal{i}(j)],'Color',PlotColors{j});
        plot([0 0]+j*0.005,tauft_diurnal{i}(j)+sigma_tauft_diurnal{i}(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
        plot([1 1]-j*0.005,tauit_diurnal{i}(j)+sigma_tauit_diurnal{i}(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
    end

    %annotate plot
    legend(diurnalrange_labels,'Location','NorthEast');
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('TFEM threshold, $$\tau_{th}$$','interpreter','latex');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
    print([folder_Plots,'tauth_fQ_diurnal_',Sites{i},'.png'],'-dpng');
end

%% Separate data according to dates

%analysis is only for Oceano
i = 3;
    
%get indices for date ranges
date_ind_window = zeros(size(Date_window{i})); %intialize list of indices for date range
for j = 1:N_daterange
    ind_range = find(Date_window{i}>=daterange_start(j)&Date_window{i}<daterange_end(j)); %get indices of windows in date range
    date_ind_window(ind_range)=j; %assign index for date range
end

%separate fQ and tauth into date range, generate binned values
fQ_date = cell(N_daterange,1);
fQ_date_bin_avg = cell(N_daterange,1);
fQ_date_bin_SE = cell(N_daterange,1);
tauth_date = cell(N_daterange,1);
tauth_date_bin_avg = cell(N_daterange,1);
tauth_date_bin_SE = cell(N_daterange,1);
tauit_date = zeros(N_daterange,1);
sigma_tauit_date = zeros(N_daterange,1);
tauft_date = zeros(N_daterange,1);
sigma_tauft_date = zeros(N_daterange,1);
for j=1:N_daterange
    %get values
    fQ_date{j} = fQ_avg_window{i}(date_ind_window==j,ind_delta_t_avg_analysis);
    tauth_date{j} = tauth_fQ_avg_window{i}(date_ind_window==j,ind_delta_t_avg_analysis);

    fQ_date_binning = fQ_date{j}(fQ_date{j}>=fQ_min&fQ_date{j}<=fQ_max);
    tauth_date_binning = tauth_date{j}(fQ_date{j}>=fQ_min&fQ_date{j}<=fQ_max);

    %perform binning
    if ~isempty(fQ_date_binning)
        [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_date_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
        N_fQ_bins = length(fQ_bin_min);
        fQ_date_bin_avg{j} = fQ_bin_avg;
        fQ_date_bin_SE{j} = fQ_bin_SE;
        tauth_date_bin_avg{j} = zeros(N_fQ_bins,1);
        tauth_date_bin_SE{j} = zeros(N_fQ_bins,1);
        for k=1:N_fQ_bins   
            fQ_date_bin_ind = find(fQ_date_binning>=fQ_bin_min(k)&fQ_date_binning<=fQ_bin_max(k));
            tauth_date_bin_avg{j}(k) = mean(tauth_date_binning(fQ_date_bin_ind));
            tauth_date_bin_SE{j}(k) = std(tauth_date_binning(fQ_date_bin_ind))/sqrt(length(fQ_date_bin_ind));
        end
    end

    %get best fit values
    [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = linearfit(fQ_date_bin_avg{j}, tauth_date_bin_avg{j}, tauth_date_bin_SE{j});
    tauft_date(j) = a;
    sigma_tauft_date(j) = sigma_a;
    tauit_date(j) = a+b;
    sigma_tauit_date(j) = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);
end

%plot binned values of fQ and tuath in separate time of day ranges
figure(2); clf; hold on; %initialize plot

%plot average values
for j = 1:N_daterange
    plot(fQ_date_bin_avg{j},tauth_date_bin_avg{j},PlotMarkers{j},'Color',PlotColors{j});
end

%plot SE values
for j = 1:N_daterange
    N_fQ_bins = length(fQ_date_bin_avg{j});
    for k = 1:N_fQ_bins
        plot(fQ_date_bin_avg{j}(k)+fQ_date_bin_SE{j}(k)*[-1 1],tauth_date_bin_avg{j}(k)*[1 1],'Color',PlotColors{j});
        plot(fQ_date_bin_avg{j}(k)*[1 1],tauth_date_bin_avg{j}(k)+tauth_date_bin_SE{j}(k)*[-1 1],'Color',PlotColors{j});
    end
end

%plot fit values
for j = 1:N_daterange
    plot([0 1],[tauft_date(j) tauit_date(j)],'Color',PlotColors{j});
    plot([0 0]+j*0.005,tauft_date(j)+sigma_tauft_date(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
    plot([1 1]-j*0.005,tauit_date(j)+sigma_tauit_date(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
end

%annotate plot
legend(daterange_labels,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('TFEM threshold, $$\tau_{th}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'tauth_fQ_date_',Sites{i},'.png'],'-dpng');

%% compare activity versus delta t - window-averaged values

%generate binned values
fQ_delta_t_avg_bin_avg = cell(N_Sites,1);
fQ_delta_t_avg_bin_SE = cell(N_Sites,1);
tauth_delta_t_avg_bin_avg = cell(N_Sites,1);
tauth_delta_t_avg_bin_SE = cell(N_Sites,1);
tauft_delta_t_avg = cell(N_Sites,1);
sigma_tauft_delta_t_avg = cell(N_Sites,1);
tauit_delta_t_avg = cell(N_Sites,1);
sigma_tauit_delta_t_avg = cell(N_Sites,1);
for i = 1:N_Sites
    fQ_delta_t_avg_bin_avg{i} = cell(N_delta_t_avg,1);
    fQ_delta_t_avg_bin_SE{i} = cell(N_delta_t_avg,1);
    tauth_delta_t_avg_bin_avg{i} = cell(N_delta_t_avg,1);
    tauth_delta_t_avg_bin_SE{i} = cell(N_delta_t_avg,1);
    tauit_delta_t_avg{i} = zeros(N_delta_t_avg,1);
    sigma_tauit_delta_t_avg{i} = zeros(N_delta_t_avg,1);
    tauft_delta_t_avg{i} = zeros(N_delta_t_avg,1);
    sigma_tauft_delta_t_avg{i} = zeros(N_delta_t_avg,1);
    for j=1:N_delta_t_avg
        %get values
        fQ_delta_t = fQ_avg_window{i}(:,j);
        tauth_delta_t = tauth_fQ_avg_window{i}(:,j);

        fQ_delta_t_binning = fQ_delta_t(fQ_delta_t>=fQ_min&fQ_delta_t<=fQ_max);
        tauth_delta_t_binning = tauth_delta_t(fQ_delta_t>=fQ_min&fQ_delta_t<=fQ_max);
        
        %perform binning
        if ~isempty(fQ_delta_t_binning)
            [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_delta_t_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
            N_fQ_bins = length(fQ_bin_min);
            fQ_delta_t_avg_bin_avg{i}{j} = fQ_bin_avg;
            fQ_delta_t_avg_bin_SE{i}{j} = fQ_bin_SE;
            tauth_delta_t_avg_bin_avg{i}{j} = zeros(N_fQ_bins,1);
            tauth_delta_t_avg_bin_SE{i}{j} = zeros(N_fQ_bins,1);
            for k=1:N_fQ_bins   
                fQ_delta_t_bin_ind = find(fQ_delta_t_binning>=fQ_bin_min(k)&fQ_delta_t_binning<=fQ_bin_max(k));
                tauth_delta_t_avg_bin_avg{i}{j}(k) = mean(tauth_delta_t_binning(fQ_delta_t_bin_ind));
                tauth_delta_t_avg_bin_SE{i}{j}(k) = std(tauth_delta_t_binning(fQ_delta_t_bin_ind))/sqrt(length(fQ_delta_t_bin_ind));
            end
        end
        
        %get best fit values
        [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = linearfit(fQ_delta_t_avg_bin_avg{i}{j}, tauth_delta_t_avg_bin_avg{i}{j}, tauth_delta_t_avg_bin_SE{i}{j});
        tauft_delta_t_avg{i}(j) = a;
        sigma_tauft_delta_t_avg{i}(j) = sigma_a;
        tauit_delta_t_avg{i}(j) = a+b;
        sigma_tauit_delta_t_avg{i}(j) = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);
    end
end

%plot activity versus delta t
for i = 1:N_Sites; %plot all sites
    figure(3); clf; hold on; %initialize plot

    %plot average values
    for j = 1:N_delta_t_avg
        plot(fQ_delta_t_avg_bin_avg{i}{j},tauth_delta_t_avg_bin_avg{i}{j},PlotMarkers{j},'Color',PlotColors{j});
    end

    %plot intercept values
    plot([1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'k','LineWidth',2);
    
    %plot SE values
    for j = 1:N_delta_t_avg
        N_fQ_bins = length(fQ_delta_t_avg_bin_avg{i}{j});
        for k = 1:N_fQ_bins
            plot(fQ_delta_t_avg_bin_avg{i}{j}(k)+fQ_delta_t_avg_bin_SE{i}{j}(k)*[-1 1],tauth_delta_t_avg_bin_avg{i}{j}(k)*[1 1],'Color',PlotColors{j});
            plot(fQ_delta_t_avg_bin_avg{i}{j}(k)*[1 1],tauth_delta_t_avg_bin_avg{i}{j}(k)+tauth_delta_t_avg_bin_SE{i}{j}(k)*[-1 1],'Color',PlotColors{j});
        end
    end

    %plot fit values
    for j = 1:N_delta_t_avg
        plot([0 1],[tauft_delta_t_avg{i}(j) tauit_delta_t_avg{i}(j)],'Color',PlotColors{j});
        plot([0 0]+j*0.005,tauft_delta_t_avg{i}(j)+sigma_tauft_delta_t_avg{i}(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
        plot([1 1]-j*0.005,tauit_delta_t_avg{i}(j)+sigma_tauit_delta_t_avg{i}(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
    end

    %annotate plot
    legend_values = delta_t_labels; legend_values{length(legend_values)+1} = 'intercept \tau_{it}';
    legend(legend_values,'Location','NorthEast');
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('TFEM threshold, $$\tau_{th}$$','interpreter','latex');
    title('using {\delta}t window-avg wind speed');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
    print([folder_Plots,'tauth_fQ_delta_t_avg_',Sites{i},'.png'],'-dpng');
end

%% compare activity versus delta t - low-pass values (for wind)

%generate binned values
fQ_delta_t_lowpass_bin_avg = cell(N_Sites,1);
fQ_delta_t_lowpass_bin_SE = cell(N_Sites,1);
tauth_delta_t_lowpass_bin_avg = cell(N_Sites,1);
tauth_delta_t_lowpass_bin_SE = cell(N_Sites,1);
tauft_delta_t_lowpass = cell(N_Sites,1);
sigma_tauft_delta_t_lowpass = cell(N_Sites,1);
tauit_delta_t_lowpass = cell(N_Sites,1);
sigma_tauit_delta_t_lowpass = cell(N_Sites,1);
for i = 1:N_Sites
    fQ_delta_t_lowpass_bin_avg{i} = cell(N_delta_t_avg,1);
    fQ_delta_t_lowpass_bin_SE{i} = cell(N_delta_t_avg,1);
    tauth_delta_t_lowpass_bin_avg{i} = cell(N_delta_t_avg,1);
    tauth_delta_t_lowpass_bin_SE{i} = cell(N_delta_t_avg,1);
    tauit_delta_t_lowpass{i} = zeros(N_delta_t_avg,1);
    sigma_tauit_delta_t_lowpass{i} = zeros(N_delta_t_avg,1);
    tauft_delta_t_lowpass{i} = zeros(N_delta_t_avg,1);
    sigma_tauft_delta_t_lowpass{i} = zeros(N_delta_t_avg,1);
    for j=1:N_delta_t_avg
        %get values
        fQ_delta_t = fQ_avg_window{i}(:,j);
        tauth_delta_t = tauth_fQ_lowpass_window{i}(:,j);

        fQ_delta_t_binning = fQ_delta_t(fQ_delta_t>=fQ_min&fQ_delta_t<=fQ_max);
        tauth_delta_t_binning = tauth_delta_t(fQ_delta_t>=fQ_min&fQ_delta_t<=fQ_max);
        
        %perform binning
        if ~isempty(fQ_delta_t_binning)
            [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_delta_t_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
            N_fQ_bins = length(fQ_bin_min);
            fQ_delta_t_lowpass_bin_avg{i}{j} = fQ_bin_avg;
            fQ_delta_t_lowpass_bin_SE{i}{j} = fQ_bin_SE;
            tauth_delta_t_lowpass_bin_avg{i}{j} = zeros(N_fQ_bins,1);
            tauth_delta_t_lowpass_bin_SE{i}{j} = zeros(N_fQ_bins,1);
            for k=1:N_fQ_bins   
                fQ_delta_t_bin_ind = find(fQ_delta_t_binning>=fQ_bin_min(k)&fQ_delta_t_binning<=fQ_bin_max(k));
                tauth_delta_t_lowpass_bin_avg{i}{j}(k) = mean(tauth_delta_t_binning(fQ_delta_t_bin_ind));
                tauth_delta_t_lowpass_bin_SE{i}{j}(k) = std(tauth_delta_t_binning(fQ_delta_t_bin_ind))/sqrt(length(fQ_delta_t_bin_ind));
            end
        end
        
        %get best fit values
        [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = linearfit(fQ_delta_t_lowpass_bin_avg{i}{j}, tauth_delta_t_lowpass_bin_avg{i}{j}, tauth_delta_t_lowpass_bin_SE{i}{j});
        tauft_delta_t_lowpass{i}(j) = a;
        sigma_tauft_delta_t_lowpass{i}(j) = sigma_a;
        tauit_delta_t_lowpass{i}(j) = a+b;
        sigma_tauit_delta_t_lowpass{i}(j) = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);
    end
end

%plot activity versus delta t - low pass values
for i = 1:N_Sites; %plot all sites
    figure(4); clf; hold on; %initialize plot

    %plot average values
    for j = 1:N_delta_t_avg
        plot(fQ_delta_t_lowpass_bin_avg{i}{j},tauth_delta_t_lowpass_bin_avg{i}{j},PlotMarkers{j},'Color',PlotColors{j});
    end

    %plot intercept values
    plot([1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'k','LineWidth',2);
    
    %plot SE values
    for j = 1:N_delta_t_avg
        N_fQ_bins = length(fQ_delta_t_lowpass_bin_avg{i}{j});
        for k = 1:N_fQ_bins
            plot(fQ_delta_t_lowpass_bin_avg{i}{j}(k)+fQ_delta_t_lowpass_bin_SE{i}{j}(k)*[-1 1],tauth_delta_t_lowpass_bin_avg{i}{j}(k)*[1 1],'Color',PlotColors{j});
            plot(fQ_delta_t_lowpass_bin_avg{i}{j}(k)*[1 1],tauth_delta_t_lowpass_bin_avg{i}{j}(k)+tauth_delta_t_lowpass_bin_SE{i}{j}(k)*[-1 1],'Color',PlotColors{j});
        end
    end

    %plot fit values
    for j = 1:N_delta_t_avg
        plot([0 1],[tauft_delta_t_lowpass{i}(j) tauit_delta_t_lowpass{i}(j)],'Color',PlotColors{j});
        plot([0 0]+j*0.005,tauft_delta_t_lowpass{i}(j)+sigma_tauft_delta_t_lowpass{i}(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
        plot([1 1]-j*0.005,tauit_delta_t_lowpass{i}(j)+sigma_tauit_delta_t_lowpass{i}(j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
    end

    %annotate plot
    legend_values = delta_t_labels; legend_values{length(legend_values)+1} = 'intercept \tau_{it}';
    legend(legend_values,'Location','NorthEast');
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('TFEM threshold, $$\tau_{th}$$','interpreter','latex');
    title('using 0.04 Hz low-pass wind speed');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
    print([folder_Plots,'tauth_fQ_delta_t_lowpass_',Sites{i},'.png'],'-dpng');
end

%plot tauth versus fQ for window-averaged and versus low pass values
for i = 1:N_Sites; %plot all sites
    figure(5); clf; hold on; %initialize plot

    %plot average values - window-average then low-pass - use values corresponding to core analysis
    plot(fQ_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis},tauth_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis},PlotMarkers{1},'Color',PlotColors{1});
    plot(fQ_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis},tauth_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis},PlotMarkers{2},'Color',PlotColors{2});
    
    %plot intercept values
    plot([1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'k','LineWidth',2);
    
    %plot SE values
    N_fQ_bins = length(fQ_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis});
    for k = 1:N_fQ_bins
        plot(fQ_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis}(k)+fQ_delta_t_avg_bin_SE{i}{ind_delta_t_avg_analysis}(k)*[-1 1],tauth_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis}(k)*[1 1],'Color',PlotColors{1});
        plot(fQ_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis}(k)*[1 1],tauth_delta_t_avg_bin_avg{i}{ind_delta_t_avg_analysis}(k)+tauth_delta_t_avg_bin_SE{i}{ind_delta_t_avg_analysis}(k)*[-1 1],'Color',PlotColors{1});
    end
    N_fQ_bins = length(fQ_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis});
    for k = 1:N_fQ_bins
        plot(fQ_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis}(k)+fQ_delta_t_lowpass_bin_SE{i}{ind_delta_t_lowpass_analysis}(k)*[-1 1],tauth_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis}(k)*[1 1],'Color',PlotColors{2});
        plot(fQ_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis}(k)*[1 1],tauth_delta_t_lowpass_bin_avg{i}{ind_delta_t_lowpass_analysis}(k)+tauth_delta_t_lowpass_bin_SE{i}{ind_delta_t_lowpass_analysis}(k)*[-1 1],'Color',PlotColors{2});
    end

    %plot fit values
    plot([0 1],[tauft_delta_t_avg{i}(ind_delta_t_avg_analysis) tauit_delta_t_avg{i}(ind_delta_t_avg_analysis)],'Color',PlotColors{1});
    plot([0 0]+j*0.005,tauft_delta_t_avg{i}(ind_delta_t_avg_analysis)+sigma_tauft_delta_t_avg{i}(ind_delta_t_avg_analysis)*[-1 1],'Color',PlotColors{1},'LineWidth',2);
    plot([1 1]-j*0.005,tauit_delta_t_avg{i}(ind_delta_t_avg_analysis)+sigma_tauit_delta_t_avg{i}(ind_delta_t_avg_analysis)*[-1 1],'Color',PlotColors{1},'LineWidth',2);
    plot([0 1],[tauft_delta_t_lowpass{i}(ind_delta_t_lowpass_analysis) tauit_delta_t_lowpass{i}(ind_delta_t_lowpass_analysis)],'Color',PlotColors{2});
    plot([0 0]+j*0.005,tauft_delta_t_lowpass{i}(ind_delta_t_lowpass_analysis)+sigma_tauft_delta_t_lowpass{i}(ind_delta_t_lowpass_analysis)*[-1 1],'Color',PlotColors{2},'LineWidth',2);
    plot([1 1]-j*0.005,tauit_delta_t_lowpass{i}(ind_delta_t_lowpass_analysis)+sigma_tauit_delta_t_lowpass{i}(ind_delta_t_lowpass_analysis)*[-1 1],'Color',PlotColors{2},'LineWidth',2);
 
    %annotate plot
    legend_values = {'1 s window-averaged u','0.04 Hz low-pass u','intercept \tau_{it}'};
    legend(legend_values,'Location','NorthEast');
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('TFEM threshold, $$\tau_{th}$$','interpreter','latex');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
    print([folder_Plots,'tauth_fQ_avg_vs_lowpass_',Sites{i},'.png'],'-dpng');
end

%% plot activity versus wind sector events for Oceano - using window averaged values

%analysis is only for Oceano
i = 3;

%get fQ's and tauth's only for Oceano and only delta_t of interest - using values from window averages
fQ_core_analysis_avg = fQ_avg_window{i}(:,ind_delta_t_avg_analysis);
tauth_core_analysis_avg = tauth_fQ_avg_window{i}(:,ind_delta_t_avg_analysis);
fplus_core_analysis_avg = fplus_avg_window{i}(:,ind_delta_t_avg_analysis);
fint_core_analysis_avg = fint_avg_window{i}(:,ind_delta_t_avg_analysis);
fint_down_core_analysis_avg = fint_down_avg_window{i}(:,ind_delta_t_avg_analysis);

%only include values in certain fQ range for binning - using values from window averages 
fQ_binning_avg = fQ_core_analysis_avg(fQ_core_analysis_avg>=fQ_min&fQ_core_analysis_avg<=fQ_max);
tauth_binning_avg = tauth_core_analysis_avg(fQ_core_analysis_avg>=fQ_min&fQ_core_analysis_avg<=fQ_max);
fplus_binning_avg = fplus_core_analysis_avg(fQ_core_analysis_avg>=fQ_min&fQ_core_analysis_avg<=fQ_max);
fint_binning_avg = fint_core_analysis_avg(fQ_core_analysis_avg>=fQ_min&fQ_core_analysis_avg<=fQ_max);
fint_down_binning_avg = fint_down_core_analysis_avg(fQ_core_analysis_avg>=fQ_min&fQ_core_analysis_avg<=fQ_max); 

%perform binning
[~, fQ_bin_N, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_binning_avg, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
N_fQ_bins = length(fQ_bin_min);
tauth_bin_avg = zeros(N_fQ_bins,1);
tauth_bin_SE = zeros(N_fQ_bins,1);
fQpred_ft_bin_avg = zeros(N_fQ_bins,1);
fQpred_ft_bin_SE = zeros(N_fQ_bins,1);
fQpred_it_bin_avg = zeros(N_fQ_bins,1);
fQpred_it_bin_SE = zeros(N_fQ_bins,1);
fQpred_hyst_bin_avg = zeros(N_fQ_bins,1);
fQpred_hyst_bin_SE = zeros(N_fQ_bins,1);
for k=1:N_fQ_bins   
    fQ_bin_ind = find(fQ_binning_avg>=fQ_bin_min(k)&fQ_binning_avg<=fQ_bin_max(k));
    tauth_bin_avg(k) = mean(tauth_binning_avg(fQ_bin_ind));
    tauth_bin_SE(k) = std(tauth_binning_avg(fQ_bin_ind))/sqrt(fQ_bin_N(k));
    fQpred_ft = fplus_binning_avg(fQ_bin_ind);
    fQpred_ft_bin_avg(k) = mean(fQpred_ft);
    fQpred_ft_bin_SE(k) = std(fQpred_ft)/sqrt(fQ_bin_N(k));
    fQpred_it = fplus_binning_avg(fQ_bin_ind)+fint_binning_avg(fQ_bin_ind);
    fQpred_it_bin_avg(k) = mean(fQpred_it);
    fQpred_it_bin_SE(k) = std(fQpred_it)/sqrt(fQ_bin_N(k));
    fQpred_hyst = fplus_binning_avg(fQ_bin_ind)+fint_down_binning_avg(fQ_bin_ind);
    fQpred_hyst_bin_avg(k) = mean(fQpred_hyst);
    fQpred_hyst_bin_SE(k) = std(fQpred_hyst)/sqrt(fQ_bin_N(k));
end

%plot binned values of fQ and wind ranges
figure(6); clf; hold on; %initialize plot

%plot mean predicted values
plot(fQ_bin_avg,fQpred_ft_bin_avg,PlotMarkers{1},'Color',PlotColors{1}); %fluid threshold prediction
plot(fQ_bin_avg,fQpred_it_bin_avg,PlotMarkers{2},'Color',PlotColors{2}); %impact threshold prediction
plot(fQ_bin_avg,fQpred_hyst_bin_avg,PlotMarkers{3},'Color',PlotColors{3}); %hysteresis prediction

%plot uncertainties
for k = 1:N_fQ_bins
    plot(fQ_bin_avg(k)*[1 1],fQpred_ft_bin_avg(k)+fQpred_ft_bin_SE(k)*[-1 1],'Color',PlotColors{1}); %fluid threshold prediction
    plot(fQ_bin_avg(k)+fQ_bin_SE(k)*[-1 1],fQpred_ft_bin_avg(k)*[1 1],'Color',PlotColors{1}); %fluid threshold prediction

    plot(fQ_bin_avg(k)*[1 1],fQpred_it_bin_avg(k)+fQpred_it_bin_SE(k)*[-1 1],'Color',PlotColors{2}); %impact threshold prediction
    plot(fQ_bin_avg(k)+fQ_bin_SE(k)*[-1 1],fQpred_it_bin_avg(k)*[1 1],'Color',PlotColors{2}); %impact threshold prediction

    plot(fQ_bin_avg(k)*[1 1],fQpred_hyst_bin_avg(k)+fQpred_hyst_bin_SE(k)*[-1 1],'Color',PlotColors{3}); %hysteresis prediction
    plot(fQ_bin_avg(k)+fQ_bin_SE(k)*[-1 1],fQpred_hyst_bin_avg(k)*[1 1],'Color',PlotColors{3}); %hysteresis prediction
end
    
plot([0,1],[0,1],'k--');
legend({'fluid threshold','impact threshold','hysteresis'},'Location','NorthWest');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('predicted activity, $$f_{Q,pred}$$','interpreter','latex');
title('using 1 s window-averaged wind speed');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'fQ_pred_delta_t_avg_',Sites{i},'.png'],'-dpng');

%% plot activity versus wind sector events for Oceano - using low-pass values

%analysis is only for Oceano
i = 3;

%get fQ's and tauth's only for Oceano and only delta_t of interest - using values from window averages
fQ_core_analysis_lowpass = fQ_avg_window{i}(:,ind_delta_t_lowpass_analysis);
tauth_core_analysis_lowpass = tauth_fQ_lowpass_window{i}(:,ind_delta_t_lowpass_analysis);
fplus_core_analysis_lowpass = fplus_lowpass_window{i}(:,ind_delta_t_lowpass_analysis);
fint_core_analysis_lowpass = fint_lowpass_window{i}(:,ind_delta_t_lowpass_analysis);
fint_down_core_analysis_lowpass = fint_down_lowpass_window{i}(:,ind_delta_t_lowpass_analysis);

%only include values in certain fQ range for binning - using values from window averages 
fQ_binning_lowpass = fQ_core_analysis_lowpass(fQ_core_analysis_lowpass>=fQ_min&fQ_core_analysis_lowpass<=fQ_max);
tauth_binning_lowpass = tauth_core_analysis_lowpass(fQ_core_analysis_lowpass>=fQ_min&fQ_core_analysis_lowpass<=fQ_max);
fplus_binning_lowpass = fplus_core_analysis_lowpass(fQ_core_analysis_lowpass>=fQ_min&fQ_core_analysis_lowpass<=fQ_max);
fint_binning_lowpass = fint_core_analysis_lowpass(fQ_core_analysis_lowpass>=fQ_min&fQ_core_analysis_lowpass<=fQ_max);
fint_down_binning_lowpass = fint_down_core_analysis_lowpass(fQ_core_analysis_lowpass>=fQ_min&fQ_core_analysis_lowpass<=fQ_max); 

%perform binning
[~, fQ_bin_N, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_binning_lowpass, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
N_fQ_bins = length(fQ_bin_min);
tauth_bin_avg = zeros(N_fQ_bins,1);
tauth_bin_SE = zeros(N_fQ_bins,1);
fQpred_ft_bin_avg = zeros(N_fQ_bins,1);
fQpred_ft_bin_SE = zeros(N_fQ_bins,1);
fQpred_it_bin_avg = zeros(N_fQ_bins,1);
fQpred_it_bin_SE = zeros(N_fQ_bins,1);
fQpred_hyst_bin_avg = zeros(N_fQ_bins,1);
fQpred_hyst_bin_SE = zeros(N_fQ_bins,1);
for k=1:N_fQ_bins   
    fQ_bin_ind = find(fQ_binning_lowpass>=fQ_bin_min(k)&fQ_binning_lowpass<=fQ_bin_max(k));
    tauth_bin_avg(k) = mean(tauth_binning_lowpass(fQ_bin_ind));
    tauth_bin_SE(k) = std(tauth_binning_lowpass(fQ_bin_ind))/sqrt(fQ_bin_N(k));
    fQpred_ft = fplus_binning_lowpass(fQ_bin_ind);
    fQpred_ft_bin_avg(k) = mean(fQpred_ft);
    fQpred_ft_bin_SE(k) = std(fQpred_ft)/sqrt(fQ_bin_N(k));
    fQpred_it = fplus_binning_lowpass(fQ_bin_ind)+fint_binning_lowpass(fQ_bin_ind);
    fQpred_it_bin_avg(k) = mean(fQpred_it);
    fQpred_it_bin_SE(k) = std(fQpred_it)/sqrt(fQ_bin_N(k));
    fQpred_hyst = fplus_binning_lowpass(fQ_bin_ind)+fint_down_binning_lowpass(fQ_bin_ind);
    fQpred_hyst_bin_avg(k) = mean(fQpred_hyst);
    fQpred_hyst_bin_SE(k) = std(fQpred_hyst)/sqrt(fQ_bin_N(k));
end

%plot binned values of fQ and wind ranges
figure(7); clf; hold on; %initialize plot

%plot mean predicted values
plot(fQ_bin_avg,fQpred_ft_bin_avg,PlotMarkers{1},'Color',PlotColors{1}); %fluid threshold prediction
plot(fQ_bin_avg,fQpred_it_bin_avg,PlotMarkers{2},'Color',PlotColors{2}); %impact threshold prediction
plot(fQ_bin_avg,fQpred_hyst_bin_avg,PlotMarkers{3},'Color',PlotColors{3}); %hysteresis prediction

%plot uncertainties
for k = 1:N_fQ_bins
    plot(fQ_bin_avg(k)*[1 1],fQpred_ft_bin_avg(k)+fQpred_ft_bin_SE(k)*[-1 1],'Color',PlotColors{1}); %fluid threshold prediction
    plot(fQ_bin_avg(k)+fQ_bin_SE(k)*[-1 1],fQpred_ft_bin_avg(k)*[1 1],'Color',PlotColors{1}); %fluid threshold prediction

    plot(fQ_bin_avg(k)*[1 1],fQpred_it_bin_avg(k)+fQpred_it_bin_SE(k)*[-1 1],'Color',PlotColors{2}); %impact threshold prediction
    plot(fQ_bin_avg(k)+fQ_bin_SE(k)*[-1 1],fQpred_it_bin_avg(k)*[1 1],'Color',PlotColors{2}); %impact threshold prediction

    plot(fQ_bin_avg(k)*[1 1],fQpred_hyst_bin_avg(k)+fQpred_hyst_bin_SE(k)*[-1 1],'Color',PlotColors{3}); %hysteresis prediction
    plot(fQ_bin_avg(k)+fQ_bin_SE(k)*[-1 1],fQpred_hyst_bin_avg(k)*[1 1],'Color',PlotColors{3}); %hysteresis prediction
end
    
plot([0,1],[0,1],'k--');
legend({'fluid threshold','impact threshold','hysteresis'},'Location','NorthWest');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('predicted activity, $$f_{Q,pred}$$','interpreter','latex');
title('using 0.04 Hz low-pass wind speed');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'fQ_pred_delta_t_lowpass_',Sites{i},'.png'],'-dpng');