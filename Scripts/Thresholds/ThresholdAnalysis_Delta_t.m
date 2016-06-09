%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% all possible subwindow times
delta_t_subwindow_all = [duration(0,0,25),duration(0,1,0),duration(0,3,0),duration(0,5,0),duration(0,10,0),duration(0,30,0)];
N_delta_t_subwindow = length(delta_t_subwindow_all);

%% plotting info
PlotFont = 12;
PlotMarkers = {'s','d','o','<','>'};
PlotColors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};

%% site info
tauit_intercept = [0.137, 0.106, 0.086];
sigma_tauit_intercept = [0.015, 0.020, 0.008];

%% info about time scales for core analysis
delta_t_avg_analysis = duration(0,0,1); %delta_t for window-avg analysis
delta_t_lowpass_analysis = duration(0,0,25); %delta_t for lowpass analysis

%% diurnal info
diurnalrange_start_hour = [0; 14; 16];
diurnalrange_end_hour = [14; 16; 24];
diurnalrange_labels = {'midday','early afternoon','late afternoon'};
N_diurnalrange = length(diurnalrange_start_hour);

%% binning info
fQ_min = 0.05;
fQ_max = 0.95;
fQ_bin_minrange = 0.1;
fQ_bin_maxrange = 0.2;
fQ_bin_N_min = 3;

%% go through all possible durations of subwindows for plotting
for m = 1:N_delta_t_subwindow
    
    %% get specific subwindow time
    delta_t_subwindow = delta_t_subwindow_all(m);

    %% information about where to load/save data, plots, and functions
    folder_LoadData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
    folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
    folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots
    folder_Functions = '../Functions/'; %folder with functions
    addpath(folder_Functions); %point MATLAB to location of functions

    %% Specific information for window data
    LoadData_1_Path = strcat(folder_LoadData,'ActivityWindows_30min_',int2str(seconds(delta_t_subwindow)),'s.mat'); %path for loading window data
    LoadData_2_Path = strcat(folder_LoadData,'TimeWindows_30min.mat'); %path for loading time windows
    SaveData_Path = strcat(folder_SaveData,'ThresholdWindows_30min_',int2str(seconds(delta_t_subwindow)),'s.mat'); %path for saving output data

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

    %% bin activity subwindow values - using window averages

    %generate binned values - threshold
    fQ_avg_bin_avg = cell(N_Sites,1); %transport activity - average
    fQ_avg_bin_SE = cell(N_Sites,1); %transport activity - standard error
    uth_avg_bin_avg = cell(N_Sites,1); %TFEM threshold wind - average
    uth_avg_bin_SE = cell(N_Sites,1); %TFEM threshold wind - standard error
    ustth_avg_bin_avg = cell(N_Sites,1); %TFEM threshold shear velocity - average
    ustth_avg_bin_SE = cell(N_Sites,1); %TFEM threshold shear velocity - standard error
    tauth_avg_bin_avg = cell(N_Sites,1); %TFEM threshold shear stress - average
    tauth_avg_bin_SE = cell(N_Sites,1); %TFEM threshold shear stress - standard error

    %generate binned values - hysteresis
    fplus_avg_bin_avg = cell(N_Sites,1); %fraction of time above fluid threshold - average
    fplus_avg_bin_SE = cell(N_Sites,1); %fraction of time above fluid threshold - standard error
    fint_avg_bin_avg = cell(N_Sites,1); %fraction of time between impact and fluid threshold - average
    fint_avg_bin_SE = cell(N_Sites,1); %fraction of time between impact and fluid threshold - standard error
    fint_down_avg_bin_avg = cell(N_Sites,1); %fraction of time between impact and fluid threshold from downward crossing - average
    fint_down_avg_bin_SE = cell(N_Sites,1); %fraction of time between impact and fluid threshold from downward crossing - standard error
    fQpred_ft_avg_bin_avg = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold - average
    fQpred_ft_avg_bin_SE = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold - standard error
    fQpred_it_avg_bin_avg = cell(N_Sites,1); %predicted transport fraction based on time above impact threshold - average
    fQpred_it_avg_bin_SE = cell(N_Sites,1); %predicted transport fraction based on time above impact threshold - standard error
    fQpred_hyst_avg_bin_avg = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold and hysteresis - average
    fQpred_hyst_avg_bin_SE = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold and hysteresis - standard error

    %generate best fit values
    tauft_avg = zeros(N_Sites,1); %inferred fluid threshold stress - best fit
    sigma_tauft_avg = zeros(N_Sites,1); %inferred fluid threshold stress - uncertainty
    tauit_avg = zeros(N_Sites,1); %inferred impact threshold stress - best fit
    sigma_tauit_avg = zeros(N_Sites,1); %inferred impact threshold stress - uncertainty
    ustitftratio_avg = zeros(N_Sites,1); %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_avg = zeros(N_Sites,1); %ratio of impact and fluid threshold shear velocities - uncertainty

    %generate binned values tagged to diurnal cycle
    fQ_avg_diurnal_bin_avg = cell(N_Sites,1); %transport activity - average
    fQ_avg_diurnal_bin_SE = cell(N_Sites,1); %transport activity - standard error
    uth_avg_diurnal_bin_avg = cell(N_Sites,1); %TFEM threshold wind - average
    uth_avg_diurnal_bin_SE = cell(N_Sites,1); %TFEM threshold wind - standard error
    ustth_avg_diurnal_bin_avg = cell(N_Sites,1); %TFEM threshold shear velocity - average
    ustth_avg_diurnal_bin_SE = cell(N_Sites,1); %TFEM threshold shear velocity - standard error
    tauth_avg_diurnal_bin_avg = cell(N_Sites,1); %TFEM threshold shear stress - average
    tauth_avg_diurnal_bin_SE = cell(N_Sites,1); %TFEM threshold shear stress - standard error

    %generate best fit values tagged to diurnal cycle
    tauft_avg_diurnal = zeros(N_Sites,N_diurnalrange); %inferred fluid threshold stress - best fit
    sigma_tauft_avg_diurnal = zeros(N_Sites,N_diurnalrange); %inferred fluid threshold stress - uncertainty
    tauit_avg_diurnal = zeros(N_Sites,N_diurnalrange); %inferred impact threshold stress - best fit
    sigma_tauit_avg_diurnal = zeros(N_Sites,N_diurnalrange); %inferred impact threshold stress - uncertainty
    ustitftratio_avg_diurnal = zeros(N_Sites,N_diurnalrange); %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_avg_diurnal = zeros(N_Sites,N_diurnalrange); %ratio of impact and fluid threshold shear velocities - uncertainty

    %go through sites
    for i = 1:N_Sites

        %get values
        timeofday = reshape(timeofday_subwindow{i},numel(timeofday_subwindow{i}),1);
        fQ_avg = reshape(fQ_avg_subwindow{i},numel(fQ_avg_subwindow{i}),1);
        uth_avg = reshape(uth_fQ_avg_subwindow{i},numel(uth_fQ_avg_subwindow{i}),1);
        fplus_avg = reshape(fplus_avg_subwindow{i},numel(fplus_avg_subwindow{i}),1);
        fint_avg = reshape(fint_avg_subwindow{i},numel(fint_avg_subwindow{i}),1);
        fint_down_avg = reshape(fint_down_avg_subwindow{i},numel(fint_down_avg_subwindow{i}),1);

        %keep only values in fQ range
        timeofday_binning = timeofday(fQ_avg>=fQ_min&fQ_avg<=fQ_max);
        fQ_avg_binning = fQ_avg(fQ_avg>=fQ_min&fQ_avg<=fQ_max);
        uth_avg_binning = uth_avg(fQ_avg>=fQ_min&fQ_avg<=fQ_max);
        fplus_avg_binning = fplus_avg(fQ_avg>=fQ_min&fQ_avg<=fQ_max);
        fint_avg_binning = fint_avg(fQ_avg>=fQ_min&fQ_avg<=fQ_max);
        fint_down_avg_binning = fint_down_avg(fQ_avg>=fQ_min&fQ_avg<=fQ_max);

        %get predicted values for fQ
        fQpred_ft_avg_binning = fplus_avg_binning;
        fQpred_it_avg_binning = fint_avg_binning + fplus_avg_binning;
        fQpred_hyst_avg_binning = fint_down_avg_binning + fplus_avg_binning;

        %% PERFORM BINNING - TOTAL
        if ~isempty(fQ_avg_binning)
            %create activity bins
            [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_avg_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
            N_fQ_bins = length(fQ_bin_min);
            fQ_avg_bin_avg{i} = fQ_bin_avg;
            fQ_avg_bin_SE{i} = fQ_bin_SE;

            %initialize binned threshold values
            uth_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            uth_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            ustth_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            ustth_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            tauth_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            tauth_avg_bin_SE{i} = zeros(N_fQ_bins,1);

            %initialize binned hysteresis values
            fplus_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            fplus_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            fint_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            fint_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            fint_down_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            fint_down_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            fQpred_ft_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            fQpred_ft_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            fQpred_it_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            fQpred_it_avg_bin_SE{i} = zeros(N_fQ_bins,1);
            fQpred_hyst_avg_bin_avg{i} = zeros(N_fQ_bins,1);
            fQpred_hyst_avg_bin_SE{i} = zeros(N_fQ_bins,1);

            for k=1:N_fQ_bins
                %generate binned threshold values
                fQ_bin_ind = find(fQ_avg_binning>=fQ_bin_min(k)&fQ_avg_binning<=fQ_bin_max(k));
                uth_avg_bin_avg{i}(k) = mean(uth_avg_binning(fQ_bin_ind));
                uth_avg_bin_SE{i}(k) = std(uth_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));

                %convert from uth to ustth and tauth
                ustth_avg_bin_avg{i}(k) = (kappa*uth_avg_bin_avg{i}(k))/log(zU_window{i}(j)/z0); %threshold shear velocity - avg
                ustth_avg_bin_SE{i}(k) = (kappa*uth_avg_bin_SE{i}(k))/log(zU_window{i}(j)/z0); %threshold shear velocity - uncertainty
                tauth_avg_bin_avg{i}(k) = rho_a*ustth_avg_bin_avg{i}(k)^2; %threshold shear stress - avg
                tauth_avg_bin_SE{i}(k) = 2*rho_a*ustth_avg_bin_avg{i}(k)*ustth_avg_bin_SE{i}(k); %threshold shear stress - uncertainty

                %generate binned hysteresis values
                fplus_avg_bin_avg{i}(k) = mean(fplus_avg_binning(fQ_bin_ind));
                fplus_avg_bin_SE{i}(k) = std(fplus_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
                fint_avg_bin_avg{i}(k) = mean(fint_avg_binning(fQ_bin_ind));
                fint_avg_bin_SE{i}(k) = std(fint_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
                fint_down_avg_bin_avg{i}(k) = mean(fint_down_avg_binning(fQ_bin_ind));
                fint_down_avg_bin_SE{i}(k) = std(fint_down_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
                fQpred_ft_avg_bin_avg{i}(k) = mean(fQpred_ft_avg_binning(fQ_bin_ind));
                fQpred_ft_avg_bin_SE{i}(k) = std(fQpred_ft_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
                fQpred_it_avg_bin_avg{i}(k) = mean(fQpred_it_avg_binning(fQ_bin_ind));
                fQpred_it_avg_bin_SE{i}(k) = std(fQpred_it_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
                fQpred_hyst_avg_bin_avg{i}(k) = mean(fQpred_hyst_avg_binning(fQ_bin_ind));
                fQpred_hyst_avg_bin_SE{i}(k) = std(fQpred_hyst_avg_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
            end
        end

        %get best fit values
        [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = linearfit(fQ_avg_bin_avg{i}, tauth_avg_bin_avg{i}, tauth_avg_bin_SE{i});
        tauft_avg(i) = a;
        sigma_tauft_avg(i) = sigma_a;
        tauit_avg(i) = a+b;
        sigma_tauit_avg(i) = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);

        %get ratio of ustit / ustft
        ustitftratio_avg(i) = sqrt(tauit_avg(i)/tauft_avg(i));
        sigma_ustitftratio_avg(i) = (1/2)*sqrt((sigma_tauft_avg(i)^2/(tauft_avg(i)*tauit_avg(i)))+...
            (sigma_tauit_avg(i)^2*tauft_avg(i)/tauit_avg(i)^3));

        %% PERFORM BINNING - DIURNAL

        %intialize values for diurnal range
        fQ_avg_diurnal_bin_avg{i} = cell(N_diurnalrange,1);
        fQ_avg_diurnal_bin_SE{i} = cell(N_diurnalrange,1);
        uth_avg_diurnal_bin_avg{i} = cell(N_diurnalrange,1);
        uth_avg_diurnal_bin_SE{i} = cell(N_diurnalrange,1);
        ustth_avg_diurnal_bin_avg{i} = cell(N_diurnalrange,1);
        ustth_avg_diurnal_bin_SE{i} = cell(N_diurnalrange,1);
        tauth_avg_diurnal_bin_avg{i} = cell(N_diurnalrange,1);
        tauth_avg_diurnal_bin_SE{i} = cell(N_diurnalrange,1);

        for j=1:N_diurnalrange
            ind_diurnalrange = find(timeofday_binning>=diurnalrange_start_hour(j) & ...
                timeofday_binning<diurnalrange_end_hour(j)); %get indices of windows in diurnal range
            %get values for binning
            fQ_avg_binning_diurnal = fQ_avg_binning(ind_diurnalrange);
            uth_avg_binning_diurnal = uth_avg_binning(ind_diurnalrange);
            if ~isempty(fQ_avg_binning_diurnal)
                [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_avg_binning_diurnal, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
                N_fQ_bins = length(fQ_bin_min);
                fQ_avg_diurnal_bin_avg{i}{j} = fQ_bin_avg;
                fQ_avg_diurnal_bin_SE{i}{j} = fQ_bin_SE;
                uth_avg_diurnal_bin_avg{i}{j} = zeros(N_fQ_bins,1);
                uth_avg_diurnal_bin_SE{i}{j} = zeros(N_fQ_bins,1);
                ustth_avg_diurnal_bin_avg{i}{j} = zeros(N_fQ_bins,1);
                ustth_avg_diurnal_bin_SE{i}{j} = zeros(N_fQ_bins,1);
                tauth_avg_diurnal_bin_avg{i}{j} = zeros(N_fQ_bins,1);
                tauth_avg_diurnal_bin_SE{i}{j} = zeros(N_fQ_bins,1);

                for k=1:N_fQ_bins   
                    fQ_bin_ind = find(fQ_avg_binning_diurnal>=fQ_bin_min(k)&fQ_avg_binning_diurnal<=fQ_bin_max(k));
                    uth_avg_diurnal_bin_avg{i}{j}(k) = mean(uth_avg_binning_diurnal(fQ_bin_ind));
                    uth_avg_diurnal_bin_SE{i}{j}(k) = std(uth_avg_binning_diurnal(fQ_bin_ind))/sqrt(length(fQ_bin_ind));

                    %convert from uth to ustth and tauth
                    ustth_avg_diurnal_bin_avg{i}{j}(k) = (kappa*uth_avg_diurnal_bin_avg{i}{j}(k))/log(zU_window{i}(j)/z0); %threshold shear velocity - avg
                    ustth_avg_diurnal_bin_SE{i}{j}(k) = (kappa*uth_avg_diurnal_bin_SE{i}{j}(k))/log(zU_window{i}(j)/z0); %threshold shear velocity - uncertainty
                    tauth_avg_diurnal_bin_avg{i}{j}(k) = rho_a*ustth_avg_diurnal_bin_avg{i}{j}(k)^2; %threshold shear stress - avg
                    tauth_avg_diurnal_bin_SE{i}{j}(k) = 2*rho_a*ustth_avg_diurnal_bin_avg{i}{j}(k)*ustth_avg_diurnal_bin_SE{i}{j}(k); %threshold shear stress - uncertainty
                end
            end

            %get best fit values
            [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = linearfit(fQ_avg_diurnal_bin_avg{i}{j}, tauth_avg_diurnal_bin_avg{i}{j}, tauth_avg_diurnal_bin_SE{i}{j});
            tauft_avg_diurnal(i,j) = a;
            sigma_tauft_avg_diurnal(i,j) = sigma_a;
            tauit_avg_diurnal(i,j) = a+b;
            sigma_tauit_avg_diurnal(i,j) = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);

            %get ratio of ustit / ustft
            ustitftratio_avg(i,j) = sqrt(tauit_avg_diurnal(i,j)/tauft_avg_diurnal(i,j));
            sigma_ustitftratio_avg(i,j) = (1/2)*sqrt((sigma_tauft_avg_diurnal(i,j)^2/(tauft_avg_diurnal(i,j)*tauit_avg_diurnal(i,j)))+...
                (sigma_tauit_avg_diurnal(i,j)^2*tauft_avg_diurnal(i,j)/tauit_avg_diurnal(i,j)^3));
        end
    end

    %% plot tau_th versus activity for sites
    figure(1); clf; hold on; %initialize plot

    %plot average values
    for i = 1:N_Sites
        plot(fQ_avg_bin_avg{i},tauth_avg_bin_avg{i},PlotMarkers{i},'Color',PlotColors{i});
    end

    %plot SE values
    for i = 1:N_Sites
        N_fQ_bins = length(fQ_avg_bin_avg{i});
        for k = 1:N_fQ_bins
            plot(fQ_avg_bin_avg{i}(k)+fQ_avg_bin_SE{i}(k)*[-1 1],tauth_avg_bin_avg{i}(k)*[1 1],'Color',PlotColors{i});
            plot(fQ_avg_bin_avg{i}(k)*[1 1],tauth_avg_bin_avg{i}(k)+tauth_avg_bin_SE{i}(k)*[-1 1],'Color',PlotColors{i});
        end
    end

    %plot fit lines
    for i = 1:N_Sites
        plot([0 1],[tauft_avg(i) tauit_avg(i)],'Color',PlotColors{i});
    end

    %plot fit intercept values
    for i = 1:N_Sites
        plot([0 0]+0.01,tauft_avg(i)+sigma_tauft_avg(i)*[-1 1],'Color',PlotColors{i},'LineWidth',2);
        plot([1 1]-0.01,tauit_avg(i)+sigma_tauit_avg(i)*[-1 1],'Color',PlotColors{i},'LineWidth',2);
    end

    %plot independent intercept thresholds
    for i = 1:N_Sites
        plot([1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'--','Color',PlotColors{i},'LineWidth',2);
    end

    %annotate plot
    legend(SiteNames,'Location','NorthEast');
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('TFEM threshold, $$\tau_{th}$$ (Pa)','interpreter','latex');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
    print([folder_Plots,'tauth_fQ_',int2str(seconds(delta_t_subwindow)),'s.png'],'-dpng');

    %% plot tau_th versus activity for sites - AGAIN FOR ARO
    figure(1); clf; hold on; %initialize plot

    %plot average values
    for i = 1:N_Sites
        plot(fQ_avg_bin_avg{i},tauth_avg_bin_avg{i},PlotMarkers{i},'Color',PlotColors{i});
    end

    %plot SE values
    for i = 1:N_Sites
        N_fQ_bins = length(fQ_avg_bin_avg{i});
        for k = 1:N_fQ_bins
            plot(fQ_avg_bin_avg{i}(k)+fQ_avg_bin_SE{i}(k)*[-1 1],tauth_avg_bin_avg{i}(k)*[1 1],'Color',PlotColors{i});
            plot(fQ_avg_bin_avg{i}(k)*[1 1],tauth_avg_bin_avg{i}(k)+tauth_avg_bin_SE{i}(k)*[-1 1],'Color',PlotColors{i});
        end
    end

    %plot fit lines
    for i = 1:N_Sites
        plot([0 1],[tauft_avg(i) tauit_avg(i)],'Color',PlotColors{i});
    end

    %plot fit intercept values
    for i = 1:N_Sites
        plot([0 0]+0.01,tauft_avg(i)+sigma_tauft_avg(i)*[-1 1],'Color',PlotColors{i},'LineWidth',2);
        plot([1 1]-0.01,tauit_avg(i)+sigma_tauit_avg(i)*[-1 1],'Color',PlotColors{i},'LineWidth',2);
    end

    %annotate plot
    legend(SiteNames,'Location','NorthEast');
    xlabel('Fraction of time that saltation is occurring');
    ylabel('Wind stress at observed threshold (Pa)');
    annotation('textarrow',[0.11 0.11],[0.3 0.43]);
    text(0.01,0.105,'fluid threshold for saltation initiation');
    annotation('textarrow',[0.98 0.98],[0.67 0.49]);
    text(0.99,0.161,'impact threshold for saltation cessation','HorizontalAlignment','Right');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
    print([folder_Plots,'tauth_fQ_',int2str(seconds(delta_t_subwindow)),'s_ARO.png'],'-dpng');
    
    %% plot binned values of fQ and tuath in separate diurnal ranges
    figure(2); clf; 
    for i = 1:N_Sites
        subplot(1,N_Sites,i); hold on; %initialize subplot

        %plot average values
        for j = 1:N_diurnalrange
            plot(fQ_avg_diurnal_bin_avg{i}{j},tauth_avg_diurnal_bin_avg{i}{j},PlotMarkers{j},'Color',PlotColors{j});
        end

        %plot SE values
        for j = 1:N_diurnalrange
            N_fQ_bins = length(fQ_avg_diurnal_bin_avg{i}{j});
            for k = 1:N_fQ_bins
                plot(fQ_avg_diurnal_bin_avg{i}{j}(k)+fQ_avg_diurnal_bin_SE{i}{j}(k)*[-1 1],tauth_avg_diurnal_bin_avg{i}{j}(k)*[1 1],'Color',PlotColors{j});
                plot(fQ_avg_diurnal_bin_avg{i}{j}(k)*[1 1],tauth_avg_diurnal_bin_avg{i}{j}(k)+tauth_avg_diurnal_bin_SE{i}{j}(k)*[-1 1],'Color',PlotColors{j});
            end
        end

        %plot fit values - diurnal
        for j = 1:N_diurnalrange
            plot([0 1],[tauft_avg_diurnal(i,j) tauit_avg_diurnal(i,j)],'Color',PlotColors{j});
            plot([0 0]+0.01,tauft_avg_diurnal(i,j)+sigma_tauft_avg_diurnal(i,j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
            plot([1 1]-0.01,tauit_avg_diurnal(i,j)+sigma_tauit_avg_diurnal(i,j)*[-1 1],'Color',PlotColors{j},'LineWidth',2);
        end

        %plot fit values - average
        plot([0 1],[tauft_avg(i) tauit_avg(i)],'Color','k','LineWidth',2);
        plot([0 0],tauft_avg(i)+sigma_tauft_avg(i)*[-1 1],'Color','k','LineWidth',2);
        plot([1 1],tauit_avg(i)+sigma_tauit_avg(i)*[-1 1],'Color','k','LineWidth',2);

        %annotate plot
        legend(diurnalrange_labels,'Location','NorthEast');
        xlabel('transport activity, $$f_Q$$','interpreter','latex');
        ylabel('TFEM threshold, $$\tau_{th}$$','interpreter','latex');
        set(gca,'FontSize',PlotFont);
        set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    end

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 5]);
    print([folder_Plots,'tauth_fQ_diurnal_',int2str(seconds(delta_t_subwindow)),'s.png'],'-dpng');

    %% plot fQ_pred versus fQ for sites

    %initialize plot
    figure(3); clf; 
    for i = 1:N_Sites
        subplot(1,N_Sites,i); clf;
    end

    %plot average values
    for i = 1:N_Sites
        subplot(1,N_Sites,i); hold on;
        plot(fQ_avg_bin_avg{i},fQpred_ft_avg_bin_avg{i},PlotMarkers{1},'Color',PlotColors{1});
        plot(fQ_avg_bin_avg{i},fQpred_it_avg_bin_avg{i},PlotMarkers{2},'Color',PlotColors{2});
        plot(fQ_avg_bin_avg{i},fQpred_hyst_avg_bin_avg{i},PlotMarkers{3},'Color',PlotColors{3});
    end

    %plot SE values
    for i = 1:N_Sites
        subplot(1,N_Sites,i); hold on;
        N_fQ_bins = length(fQ_avg_bin_avg{i});
        for k = 1:N_fQ_bins
            plot(fQ_avg_bin_avg{i}(k)+fQ_avg_bin_SE{i}(k)*[-1 1],fQpred_ft_avg_bin_avg{i}(k)*[1 1],'Color',PlotColors{1});
            plot(fQ_avg_bin_avg{i}(k)*[1 1],fQpred_ft_avg_bin_avg{i}(k)+fQpred_ft_avg_bin_SE{i}(k)*[-1 1],'Color',PlotColors{1});
            plot(fQ_avg_bin_avg{i}(k)+fQ_avg_bin_SE{i}(k)*[-1 1],fQpred_it_avg_bin_avg{i}(k)*[1 1],'Color',PlotColors{2});
            plot(fQ_avg_bin_avg{i}(k)*[1 1],fQpred_it_avg_bin_avg{i}(k)+fQpred_it_avg_bin_SE{i}(k)*[-1 1],'Color',PlotColors{2});
            plot(fQ_avg_bin_avg{i}(k)+fQ_avg_bin_SE{i}(k)*[-1 1],fQpred_hyst_avg_bin_avg{i}(k)*[1 1],'Color',PlotColors{3});
            plot(fQ_avg_bin_avg{i}(k)*[1 1],fQpred_hyst_avg_bin_avg{i}(k)+fQpred_hyst_avg_bin_SE{i}(k)*[-1 1],'Color',PlotColors{3});
        end
    end

    %plot 1-1 lines
    for i = 1:N_Sites
        subplot(1,N_Sites,i); hold on;
        plot([0,1],[0,1],'k--');
    end

    %annotate plot
    for i = 1:N_Sites
        subplot(1,N_Sites,i);
        if i==N_Sites
            legend({'fluid threshold','impact threshold','hysteresis'},'Location','SouthEast');
        end
        xlabel('transport activity, $$f_Q$$','interpreter','latex');
        ylabel('predicted activity, $$f_{Q,pred}$$','interpreter','latex');
        title(SiteNames{i});
        set(gca,'FontSize',PlotFont);
        set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    end

    %print plot
    set(gca, 'LooseInset', get(gca,'TightInset'));
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 6]);
    print([folder_Plots,'fQpred_fQ_',int2str(seconds(delta_t_subwindow)),'s.png'],'-dpng');
end