% %% LOOK AT CALIBRATION VALUES FOR WENGLORS TO DETERMINE SYSTEMATIC TRENDS
% 
% %% clear existing data and load processed data and metadata
% clearvars;
% 
% %% indicate where to find processed data
% folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/';
% 
% %% indicate where to save analysis data
% folder_SaveData = '../DataOutput/';
% 
% %% indicate where to save plots
% folder_Plots = '../PlotOutput/WenglorCalibration/';
% 
% %% information about sites for analysis
% Site_Names = {'Jericoacoara','RanchoGuadalupe','Oceano'};
% Folder_Names = {'Jericoacoara2014','RanchoGuadalupe2015','Oceano2015'};
% N_Sites = length(Site_Names);
% 
% %% load processed data and metadata for each site, keep only Wenglor and BSNE data
% Data_Wenglor = cell(N_Sites,1); %initialize cell array containing Wenglor data
% Data_BSNE = cell(N_Sites,1); %initialize cell array containing BSNE data
% Metadata = cell(N_Sites,1); %initialize cell array containing metadata
% for i = 1:N_Sites
%     ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Site_Names{i});
%     Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Site_Names{i});
%     load(ProcessedData_Path); %load processed data
%     Data_Wenglor{i} = ProcessedData.Wenglor; %get Wenglor processed data
%     Data_BSNE{i} = ProcessedData.BSNE; %get BSNE processed data
%     Metadata{i} = load(Metadata_Path); %load metadata
%     clear ProcessedData;
% end
% 
% %% get Wenglor calibration factors by site, Wenglor name, and interval number
% 
% %initialize cell arrays of values for each site
% SiteNames = {'Jeri','Rancho','Oceano'};
% z_Site = cell(N_Sites,1); %cell array with Wenglor heights for each site
% C_Site = cell(N_Sites,1); %cell array with calibration factors for each site
% counts_Site = cell(N_Sites,1); %cell array with particle counts for each site
% flux_Wenglor_Site = cell(N_Sites,1); %cell array with expected fluxes for each site
% flux_BSNE_Site = cell(N_Sites,1); %cell array of BSNE fluxes for each site
% WID_Site = cell(N_Sites,1); %cell array of Wenglor identifier
% 
% for i = 1:N_Sites
%     
%     %get start and end times of BSNE intervals
%     StartTimes_BSNE = [Data_BSNE{i}(:).StartTime];
%     EndTimes_BSNE = [Data_BSNE{i}.EndTime];
%     flux_BSNE = [Data_BSNE{i}.Q];
%     N_Intervals = length(StartTimes_BSNE);
%     
%     %get list of Wenglors
%     Wenglors = fieldnames(Data_Wenglor{i});
%     N_Wenglors = length(Wenglors);
%     
%     %initialize cell arrays of heights, calibration factors, mean particle counts, expected fluxes, and IDs for each Wenglor
%     z_Wenglor = cell(N_Wenglors,1);
%     C_Wenglor = cell(N_Wenglors,1);
%     counts_Wenglor = cell(N_Wenglors,1);
%     flux_Wenglor = cell(N_Wenglors,1);
%     WID_Wenglor = cell(N_Wenglors,1);
%     
%     %get calibration factors for Wenglors in each interval
%     for j=1:N_Wenglors
%         
%         %initialize lists of heights, calibration factors, counts, expected fluxes, and IDs for all intervals
%         z_Intervals = zeros(N_Intervals,1);
%         C_Intervals = zeros(N_Intervals,1);
%         counts_Intervals = zeros(N_Intervals,1);
%         flux_Intervals = zeros(N_Intervals,1);
%         WID_Intervals = cell(N_Intervals,1);
%         
%         for k=1:N_Intervals
%             [counts, ~, IntervalNumber, IntervalInd] = ExtractVariableTimeInterval(Data_Wenglor{i}.(Wenglors{j}),StartTimes_BSNE(k),EndTimes_BSNE(k),'n','int','int');
%                         
%             %if interval contains data
%             if ~isempty(counts)
%                 
%                 %get Wenglor height
%                 z_Intervals(k) = Data_Wenglor{i}.(Wenglors{j})(IntervalNumber(1)).InstrumentHeight.z;
%                 
%                 %get Wenglor calibration factor
%                 C = mode(Data_Wenglor{i}.(Wenglors{j})(IntervalNumber(1)).flux.qzPerCount);
%                 if ~isinf(C)
%                     C_Intervals(k) = C;
%                 else
%                     C_Intervals(k) = NaN;
%                 end
%                 
%                 %get Wenglor mean counts per second
%                 counts_Intervals(k) = mean(counts)*25;
%                 
%                 %get Wenglor expected fluxes
%                 flux_Intervals(k) = mean(Data_Wenglor{i}.(Wenglors{j})(IntervalNumber(1)).flux.qz(IntervalInd{1}));
%                 
%                 %get Wenglor ID
%                 WID_Intervals{k} = Metadata{i}.InstrumentMetadata.InstrumentID(intersect(...
%                     find(strcmp(Metadata{i}.InstrumentMetadata.Instrument,Wenglors{j})),... %match in metadata list for specific Wenglor
%                     find(Metadata{i}.InstrumentMetadata.StartTime==Data_Wenglor{i}.(Wenglors{j})(IntervalNumber(1)).StartTime))); %match in metadata list for specific start time
%             
%             %if interval does not contain data
%             else
%                 z_Intervals(k) = NaN;
%                 C_Intervals(k) = NaN;
%                 counts_Intervals(k) = NaN;
%                 flux_Intervals(k) = NaN;
%             end
%         end
%         
%         %add list of calibration factors to cell array for Wenglor
%         z_Wenglor{j} = z_Intervals;
%         C_Wenglor{j} = C_Intervals;
%         counts_Wenglor{j} = counts_Intervals;
%         flux_Wenglor{j} = flux_Intervals;
%         WID_Wenglor{j} = WID_Intervals;
%     end
% 
%     %add cell array of calibration factors for each site
%     z_Site{i} = z_Wenglor;
%     C_Site{i} = C_Wenglor;
%     counts_Site{i} = counts_Wenglor;
%     flux_Wenglor_Site{i} = flux_Wenglor;
%     flux_BSNE_Site{i} = flux_BSNE;
%     WID_Site{i} = WID_Wenglor;
% end
% 
% %% Separate calibration factors and heights by each unique Wenglor
% WIDs = cell(N_Sites,1);
% WID_z = cell(N_Sites,1);
% WID_C = cell(N_Sites,1);
% WID_counts = cell(N_Sites,1);
% WID_Wenglor_flux = cell(N_Sites,1);
% WID_BSNE_flux = cell(N_Sites,1);
% 
% for i = 1:N_Sites;
%     
%     %get all possible Wenglor ID's
%     WIDs{i} = unique(Metadata{i}.InstrumentMetadata.InstrumentID(...
%         strcmp(Metadata{i}.InstrumentMetadata.InstrumentType,'Wenglor')));
%     N_WIDs = length(WIDs{i});
%     
%     WID_z{i} = cell(N_WIDs,1);
%     WID_C{i} = cell(N_WIDs,1);
%     WID_counts{i} = cell(N_WIDs,1);
%     WID_Wenglor_flux{i} = cell(N_WIDs,1);
%     WID_BSNE_flux{i} = cell(N_WIDs,1);
%     
%     for j = 1:length(WID_Site{i})
%         for k = 1:length(WID_Site{i}{j})
%             if ~isnan(C_Site{i}{j}(k))
%                 ind_WID = find(strcmp(WIDs{i},WID_Site{i}{j}{k})); %get index of WID within list
%                 WID_z{i}{ind_WID}(length(WID_z{i}{ind_WID})+1)=z_Site{i}{j}(k); %add Wenglor height to list for this WID
%                 WID_C{i}{ind_WID}(length(WID_C{i}{ind_WID})+1)=C_Site{i}{j}(k); %add Wenglor calibration factor to list for this WID
%                 WID_counts{i}{ind_WID}(length(WID_counts{i}{ind_WID})+1)=counts_Site{i}{j}(k); %add Wenglor mean counts/second to list for this WID
%                 WID_Wenglor_flux{i}{ind_WID}(length(WID_Wenglor_flux{i}{ind_WID})+1)=flux_Wenglor_Site{i}{j}(k); %add Wenglor expected flux to list for this WID
%                 WID_BSNE_flux{i}{ind_WID}(length(WID_BSNE_flux{i}{ind_WID})+1)=flux_BSNE_Site{i}(k); %add BSNE flux to list for this WID
%             end
%         end
%     end
% end

%% plot calibration factor versus Wenglor height
figure(1); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    
    for j = 1:length(WIDs{i})
        [C_sort,sort_ind] = sort(WID_C{i}{j});
        z_sort = WID_z{i}{j}(sort_ind);
        [z_sort,sort_ind] = sort(z_sort);
        C_sort = WID_C{i}{j}(sort_ind);
        plot(z_sort,C_sort,'-o','MarkerSize',2);
    end
    %legend(WIDs{i});
    
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    
    xlim([min(horzcat(WID_z{i}{:})) max(horzcat(WID_z{i}{:}))])
    ylim([1e1 1e4]);
    xlabel('z (cm)');
    ylabel('C (q/n)');
    title(SiteNames{i});
    set(gca,'FontSize',16);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print([folder_Plots,'CalibrationFactor_Height.png'],'-dpng');
end

%% plot calibration factor versus mean counts during interval
figure(2); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    
    for j = 1:length(WIDs{i})
        [counts_sort,sort_ind] = sort(WID_counts{i}{j});
        C_sort = WID_C{i}{j}(sort_ind);
        plot(counts_sort,C_sort,'-o','MarkerSize',2);
    end
    %legend(WIDs{i});
    
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    
    xlim([min(horzcat(WID_counts{i}{:})) max(horzcat(WID_counts{i}{:}))])
    ylim([1e1 1e4]);
    xlabel('counts (n/s)');
    ylabel('C (q/n)');
    title(SiteNames{i});
    set(gca,'FontSize',16);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print([folder_Plots,'CalibrationFactor_Counts.png'],'-dpng');
end

%% plot calibration factor versus expected flux during interval
figure(3); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    
    for j = 1:length(WIDs{i})
        [flux_sort, sort_ind] = sort(WID_Wenglor_flux{i}{j});
        C_sort = WID_C{i}{j}(sort_ind);
        plot(flux_sort,C_sort,'-o','MarkerSize',2);
    end
    %legend(WIDs{i});
  
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    
    xlim([min(horzcat(WID_Wenglor_flux{i}{:})) max(horzcat(WID_Wenglor_flux{i}{:}))])
    ylim([1e1 1e4]);
    xlabel('expected flux (g/m^2/s)');
    ylabel('C (q/n)');
    title(SiteNames{i});
    set(gca,'FontSize',16);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print([folder_Plots,'CalibrationFactor_ExpectedFlux.png'],'-dpng');
end


%% plot calibration factor versus BSNE flux during interval
figure(4); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    
    for j = 1:length(WIDs{i})
        [flux_sort, sort_ind] = sort(WID_BSNE_flux{i}{j});
        C_sort = WID_C{i}{j}(sort_ind);
        plot(flux_sort,C_sort,'-o','MarkerSize',2);
    end
    %legend(WIDs{i});
  
    set(gca,'yscale','log');
    %set(gca,'xscale','log');
    
    xlim([min(horzcat(WID_BSNE_flux{i}{:})) max(horzcat(WID_BSNE_flux{i}{:}))])
    ylim([1e1 1e4]);
    xlabel('BSNE flux (g/m/s)');
    ylabel('C (q/n)');
    title(SiteNames{i});
    set(gca,'FontSize',16);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 5])
    print([folder_Plots,'CalibrationFactor_BSNE_Flux.png'],'-dpng');
end


%% Separate calibration factors into height bins
N_zbins = 10;
zbins_edges = logspace(log10(0.02),log10(0.48),N_zbins+1);
zbins_min = zbins_edges(1:N_zbins);
zbins_max = zbins_edges(2:N_zbins+1);
zbins_mid = sqrt(zbins_min.*zbins_max);

zbins_C = cell(N_Sites,1);
zbins_counts = cell(N_Sites,1);
zbins_BSNE_flux = cell(N_Sites,1);

for i = 1:N_Sites;
    zbins_C{i} = cell(N_zbins,1);
    zbins_counts{i} = cell(N_zbins,1);
    zbins_BSNE_flux{i} = cell(N_zbins,1);
    for j = 1:length(z_Site{i})
        for k = 1:length(z_Site{i}{j})
            if ~isnan(C_Site{i}{j}(k))
                ind_z = intersect(find(zbins_min<=z_Site{i}{j}(k)),find(zbins_max>=z_Site{i}{j}(k))); %get index of z within list
                if ~isempty(ind_z)
                    zbins_C{i}{ind_z}(length(zbins_C{i}{ind_z})+1)=C_Site{i}{j}(k); %add Wenglor calibration factor to list for this z bin
                    zbins_counts{i}{ind_z}(length(zbins_counts{i}{ind_z})+1)=counts_Site{i}{j}(k); %add mean counts to list for this z bin
                    zbins_BSNE_flux{i}{ind_z}(length(zbins_BSNE_flux{i}{ind_z})+1)=flux_BSNE_Site{i}(k); %add overall BSNE flux to list for this z bin
                end
            end
        end
    end
end

%% plot calibration factor versus mean counts during interval for each z bin
figure(5); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    
    zbins_legend = {};
    for j = 1:length(zbins_mid)
        if ~isempty(zbins_counts{i}{j})
            [counts_sort,sort_ind] = sort(zbins_counts{i}{j});
            C_sort = zbins_C{i}{j}(sort_ind);
            plot(counts_sort,C_sort,'-o','MarkerSize',2);
            zbins_legend{length(zbins_legend)+1} = ['z = ',num2str(round(zbins_mid(j)*1000)/10),' cm'];
        end
    end
    h_legend = legend(zbins_legend);
    if i == 2;
        set(h_legend,'Location','SouthWest');
    end
    
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    
    xlim([min(horzcat(zbins_counts{i}{:})) max(horzcat(zbins_counts{i}{:}))])
    ylim([1e1 1e4]);
    xlabel('counts (n/s)');
    ylabel('C (q/n)');
    title(SiteNames{i});
    set(gca,'FontSize',12);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 5])
    print([folder_Plots,'CalibrationFactor_Height_Counts.png'],'-dpng');
end

%% plot calibration factor versus BSNE flux during interval for each z bin
figure(6); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    
    zbins_legend = {};
    for j = 1:length(zbins_mid)
        if ~isempty(zbins_counts{i}{j})
            [BSNE_flux_sort,sort_ind] = sort(zbins_BSNE_flux{i}{j});
            C_sort = zbins_C{i}{j}(sort_ind);
            plot(BSNE_flux_sort,C_sort,'-o','MarkerSize',2);
            zbins_legend{length(zbins_legend)+1} = ['z = ',num2str(round(zbins_mid(j)*1000)/10),' cm'];
        end
    end
    h_legend = legend(zbins_legend);
    if i == 2;
        set(h_legend,'Location','SouthWest');
    end
    
    set(gca,'yscale','log');
%    set(gca,'xscale','log');
    
    xlim([min(horzcat(zbins_BSNE_flux{i}{:})) max(horzcat(zbins_BSNE_flux{i}{:}))])
    ylim([1e1 1e4]);
    xlabel('BSNE flux (g/m/s)');
    ylabel('C (q/n)');
    title(SiteNames{i});
    set(gca,'FontSize',12);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 5])
    print([folder_Plots,'CalibrationFactor_Height_BSNE_flux.png'],'-dpng');
end

%% plot calibration factor versus mean counts during interval for each z bin, only for Oceano
figure(7); clf; hold on;
    
i = 3;
zbins_legend = {};
for j = 1:length(zbins_mid)
    if ~isempty(zbins_counts{i}{j})
        [counts_sort,sort_ind] = sort(zbins_counts{i}{j});
        C_sort = zbins_C{i}{j}(sort_ind);
        plot(counts_sort,C_sort,'-o','MarkerSize',2);
        zbins_legend{length(zbins_legend)+1} = ['z = ',num2str(round(zbins_mid(j)*1000)/10),' cm'];
    end
end
h_legend = legend(zbins_legend,'Location','SouthWest');

set(gca,'yscale','log');
set(gca,'xscale','log');

xlim([min(horzcat(zbins_counts{i}{:})) max(horzcat(zbins_counts{i}{:}))])
ylim([min(horzcat(zbins_C{i}{:})) max(horzcat(zbins_C{i}{:}))]);
xlabel('counts (n/s)');
ylabel('C (q/n)');
title(SiteNames{i});
set(gca,'FontSize',16);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 6])
print([folder_Plots,'CalibrationFactor_Height_Counts_Oceano.png'],'-dpng');

%% pick a site and day, plot mean calibration
i = 3; %Oceano
k_June3 = [19,18,18,17,12,12,12,13,13]; %second to last value in each list - June 3
%k_day2 = [20,19,19,18,13,13,13,14,14]; %last value in each list - June 4
Wenglors = fieldnames(Data_Wenglor{i});
N_Wenglors = length(Wenglors);

figure(8); clf; hold on;
subplot(2,1,1); hold on;
subplot(2,1,2); hold on;
for j = 1:(N_Wenglors-2)
    n = Data_Wenglor{i}.(Wenglors{j})(k_June3(j)).n.int;
    n_norm = n/mean(n);
    t = Data_Wenglor{i}.(Wenglors{j})(k_June3(j)).t.int;
    [n_wa, t_wa] = window_average(n, t, duration(0,10,0));
    [n_norm_wa, t_wa] = window_average(n_norm, t, duration(0,10,0));
    subplot(2,1,1);
    plot(t_wa,n_wa*25,'LineWidth',2);
    subplot(2,1,2);
    plot(t_wa,n_norm_wa,'LineWidth',2);
end

subplot(2,1,1);
xlim([datenum(2015,6,3,12,30,0), datenum(2015,6,3,18,30,0)])
ylabel('n (counts/s)')
title('Oceano, June 3');
legend(Wenglors,'Location','NorthWest')
set(gca,'FontSize',16)

subplot(2,1,2);
xlim([datenum(2015,6,3,12,30,0), datenum(2015,6,3,18,30,0)])
ylabel('n / \mu_{n}')
title('Oceano, June 3');
legend(Wenglors,'Location','NorthWest')
set(gca,'FontSize',16)

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 9])
print([folder_Plots,'Counts_Oceano_June3.png'],'-dpng');