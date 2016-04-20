%% initialize
clearvars;

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/GrainSize/'; %folder for outputs
SaveData_Path = strcat(folder_AnalysisData,'MeanGrainSize'); %path for saving output data
folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% Information about clustered dates for Oceano
cluster_StartDate = [datetime(2015,5,15); datetime(2015,5,23); datetime(2015,6,2)];
cluster_EndDate = [datetime(2015,5,19); datetime(2015,5,28); datetime(2015,6,4)];
N_clusters = length(cluster_StartDate);
cluster_names = cell(N_clusters,1);
for i=1:N_clusters
    cluster_names{i} = [datestr(cluster_StartDate(i),'mmm dd'),' - ',datestr(cluster_EndDate(i),'mmm dd')];
end

%% plotting information
PlotFont = 12;
LineWidth_Surface = 1;

%% load grain size / processed data
GrainSizeData_all = cell(N_Sites,1);
FluxBSNE_all = cell(N_Sites,1);
for i = 1:N_Sites
    GrainSizeData_Path = strcat(folder_ProcessedData,'GrainSize_',Sites{i});
    GrainSizeData_all{i} = load(GrainSizeData_Path);

    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path);
    FluxBSNE_all{i} = ProcessedData.FluxBSNE;
    clear ProcessedData;
end


%% initialize combined surface samples
d_surface_site = cell(N_Sites,1);
dlogd_surface_site = cell(N_Sites,1);
dV_surface_site = cell(N_Sites,1);
dVdlogd_surface_site = cell(N_Sites,1);
d10_surface_site = zeros(N_Sites,1);
d50_surface_site = zeros(N_Sites,1);
d90_surface_site = zeros(N_Sites,1);
dates_surface = cell(N_Sites,1);
d10_surface_date = cell(N_Sites,1);
d50_surface_date = cell(N_Sites,1);
d90_surface_date = cell(N_Sites,1);

%% initialize combined airborne samples
d_airborne_site = cell(N_Sites,1);
dlogd_airborne_site = cell(N_Sites,1);
dV_airborne_site = cell(N_Sites,1);
dVdlogd_airborne_site = cell(N_Sites,1);
d10_airborne_site = zeros(N_Sites,1);
d50_airborne_site = zeros(N_Sites,1);
d90_airborne_site = zeros(N_Sites,1);
dates_airborne = cell(N_Sites,1);
d10_airborne_date = cell(N_Sites,1);
d50_airborne_date = cell(N_Sites,1);
d90_airborne_date = cell(N_Sites,1);

%% go through each site
for i = 1:N_Sites
    
    %% GET SURFACE SAMPLES
    GrainSize_Surface = GrainSizeData_all{i}.GrainSize_Surface;
    
    %get size bins from first sample
    d_surface = [GrainSize_Surface(1).gsd(2:end-1).Sizeclass_mid_mm];
    d_surface_lower = [GrainSize_Surface(1).gsd(2:end-1).Sizeclass_lower_mm];
    d_surface_upper = [GrainSize_Surface(1).gsd(2:end-1).Sizeclass_upper_mm];
    dlogd_surface = log(d_surface_upper) - log(d_surface_lower);
    
    %initialize matrix of surface samples
    N_surface = length(GrainSize_Surface);
    N_d = length(d_surface);
    dV_surface = zeros(N_surface,N_d);
    dVdlogd_surface = zeros(N_surface,N_d);
    
    %get each surface size distribution
    for j = 1:N_surface
        dV_surface(j,:) = [GrainSize_Surface(j).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
        dVdlogd_surface(j,:) = dV_surface(j,:)./dlogd_surface; %by normalized volume
    end
    
    %% SURFACE - GET SITE MEAN
    dV_surface_bar = mean(dV_surface);
    dVdlogd_surface_bar = mean(dVdlogd_surface);

    %add to cell array
    d_surface_site{i} = d_surface;
    dlogd_surface_site{i} = dlogd_surface;
    dV_surface_site{i} = dV_surface_bar;
    dVdlogd_surface_site{i} = dVdlogd_surface_bar;
    
    %compute benchmark grain sizes
    CV = cumsum(dV_surface_bar); %cumulative distribution by volume
    
    %get d10
    ind_d10 = find(CV>=0.1, 1);
    d_below = d_surface_lower(ind_d10);
    d_above = d_surface_upper(ind_d10);
    wt_below = (CV(ind_d10)-0.1)/dV_surface_bar(ind_d10);
    wt_above = 1-wt_below;
    d10_surface_site(i) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    
    %get d50
    ind_d50 = find(CV>=0.5, 1);
    d_below = d_surface_lower(ind_d50);
    d_above = d_surface_upper(ind_d50);
    wt_below = (CV(ind_d50)-0.5)/dV_surface_bar(ind_d50);
    wt_above = 1-wt_below;
    d50_surface_site(i) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    
    %get d90
    ind_d90 = find(CV>=0.9, 1);
    d_below = d_surface_lower(ind_d90);
    d_above = d_surface_upper(ind_d90);
    wt_below = (CV(ind_d90)-0.9)/dV_surface_bar(ind_d90);
    wt_above = 1-wt_below;
    d90_surface_site(i) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    
    
    %% SURFACE - GET DAILY MEAN
    dates_samples = [GrainSize_Surface.Date];
    dates_unique = unique(dates_samples);
    N_dates = length(dates_unique);
    dV_surface_date = zeros(N_dates,N_d);
    dVdlogd_surface_date = zeros(N_dates,N_d);
    for j = 1:N_dates
        ind_date = find(dates_samples==dates_unique(j));
        dV_surface_date(j,:) = mean(dV_surface(ind_date,:),1);
        dVdlogd_surface_date(j,:) = mean(dVdlogd_surface(ind_date,:),1);
    end

    %get daily d10, d50, d90
    CV = cumsum(dV_surface_date')'; %cumulative distribution by volume
    dates_surface{i} = dates_unique';
    d10_surface_date{i} = zeros(N_dates,1);
    d50_surface_date{i} = zeros(N_dates,1);
    d90_surface_date{i} = zeros(N_dates,1);
    
    for j = 1:N_dates
        %get d10
        ind_d10 = find(CV(j,:)>=0.1, 1);
        d_below = d_surface_lower(ind_d10);
        d_above = d_surface_upper(ind_d10);
        wt_below = (CV(j,ind_d10)-0.1)/dV_surface_date(j,ind_d10);
        wt_above = 1-wt_below;
        d10_surface_date{i}(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

        %get d50
        ind_d50 = find(CV(j,:)>=0.5, 1);
        d_below = d_surface_lower(ind_d50);
        d_above = d_surface_upper(ind_d50);
        wt_below = (CV(j,ind_d50)-0.5)/dV_surface_date(j,ind_d50);
        wt_above = 1-wt_below;
        d50_surface_date{i}(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

        %get d90
        ind_d90 = find(CV(j,:)>=0.9, 1);
        d_below = d_surface_lower(ind_d90);
        d_above = d_surface_upper(ind_d90);
        wt_below = (CV(j,ind_d90)-0.9)/dV_surface_date(j,ind_d90);
        wt_above = 1-wt_below;
        d90_surface_date{i}(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    end
    
    
    %% SURFACE - GET CLUSTER MEAN FOR OCEANO
    if i == 3;
        d_surface_cluster = d_surface;
        dlogd_surface_cluster = dlogd_surface;
        dV_surface_cluster = zeros(N_clusters,N_d);
        dVdlogd_surface_cluster = zeros(N_clusters,N_d);
        for j = 1:N_clusters
            ind_cluster = intersect(find(dates_samples>=cluster_StartDate(j)),find(dates_samples<=cluster_EndDate(j)));
            dV_surface_cluster(j,:) = mean(dV_surface(ind_cluster,:),1);
            dVdlogd_surface_cluster(j,:) = mean(dVdlogd_surface(ind_cluster,:),1);
        end

        %get daily d10, d50, d90
        CV = cumsum(dV_surface_cluster')'; %cumulative distribution by volume
        d10_surface_cluster = zeros(N_clusters,1);
        d50_surface_cluster = zeros(N_clusters,1);
        d90_surface_cluster = zeros(N_clusters,1);

        for j = 1:N_clusters
            %get d10
            ind_d10 = find(CV(j,:)>=0.1, 1);
            d_below = d_surface_lower(ind_d10);
            d_above = d_surface_upper(ind_d10);
            wt_below = (CV(j,ind_d10)-0.1)/dV_surface_cluster(j,ind_d10);
            wt_above = 1-wt_below;
            d10_surface_cluster(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

            %get d50
            ind_d50 = find(CV(j,:)>=0.5, 1);
            d_below = d_surface_lower(ind_d50);
            d_above = d_surface_upper(ind_d50);
            wt_below = (CV(j,ind_d50)-0.5)/dV_surface_cluster(j,ind_d50);
            wt_above = 1-wt_below;
            d50_surface_cluster(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

            %get d90
            ind_d90 = find(CV(j,:)>=0.9, 1);
            d_below = d_surface_lower(ind_d90);
            d_above = d_surface_upper(ind_d90);
            wt_below = (CV(j,ind_d90)-0.9)/dV_surface_cluster(j,ind_d90);
            wt_above = 1-wt_below;
            d90_surface_cluster(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));
        end
    end
    
    
    %% BSNE SAMPLES
    GrainSize_BSNE = GrainSizeData_all{i}.GrainSize_BSNE;
    FluxBSNE = FluxBSNE_all{i};
    
    %get indices of airborne samples
    ind_airborne = []; %initialize list of airborne sample indices
    N_FluxBSNE = length(FluxBSNE);
    for j = 1:N_FluxBSNE
        dz_BSNE = abs(FluxBSNE(j).z.z-FluxBSNE(j).z.zq); %get difference in BSNE heights from zq for time interval
        Name_dzmin_BSNE = FluxBSNE(j).name(dz_BSNE==min(dz_BSNE),1); %get name of BSNE with minimal distance difference (take only first if there is more than one)
        StartTime_BSNE = FluxBSNE(j).StartTime; %get start time for BSNE interval
        ind_GrainSize_BSNE = find([GrainSize_BSNE.StartTime]<=StartTime_BSNE & [GrainSize_BSNE.EndTime]>StartTime_BSNE); %get indices within GrainSize_BSNE corresponding to time interval
        ind_airborne_BSNE = ind_GrainSize_BSNE(strcmp({GrainSize_BSNE(ind_GrainSize_BSNE).NameBSNE},Name_dzmin_BSNE)==1); %get index within GrainSize_BSNE for airborne sample
        ind_airborne = [ind_airborne; ind_airborne_BSNE]; %add to list of airborne indices
    end
    ind_airborne = unique(ind_airborne); %remove repeated values
    
    %get size bins from first sample
    d_airborne = [GrainSize_BSNE(ind_airborne(1)).gsd(2:end-1).Sizeclass_mid_mm]; 
    d_airborne_lower = [GrainSize_Surface(1).gsd(2:end-1).Sizeclass_lower_mm];
    d_airborne_upper = [GrainSize_Surface(1).gsd(2:end-1).Sizeclass_upper_mm];
    dlogd_airborne = log(d_airborne_upper) - log(d_airborne_lower);
    
    %initialize matrix of airborne samples
    N_airborne = length(ind_airborne);
    N_d = length(d_airborne);
    dV_airborne = zeros(N_airborne,N_d);
    dVdlogd_airborne = zeros(N_airborne,N_d);
    
    %get each airborne size distribution
    for j = 1:N_airborne
        dV_airborne(j,:) = [GrainSize_BSNE(ind_airborne(j)).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
        dVdlogd_airborne(j,:) = dV_airborne(j,:)./dlogd_airborne; %by normalized volume
    end
    
    %get mean of airborne sizes
    dV_airborne_bar = mean(dV_airborne);
    dVdlogd_airborne_bar = mean(dVdlogd_airborne);
    
    %add to cell array
    d_airborne_site{i} = d_airborne;
    dlogd_airborne_site{i} = dlogd_airborne;
    dV_airborne_site{i} = dV_airborne_bar;
    dVdlogd_airborne_site{i} = dVdlogd_airborne_bar;
    
    %compute benchmark grain sizes
    CV = cumsum(dV_airborne_bar); %cumulative distribution by volume
    
    %get d10
    ind_d10 = find(CV>=0.1, 1);
    d_below = d_airborne_lower(ind_d10);
    d_above = d_airborne_upper(ind_d10);
    wt_below = (CV(ind_d10)-0.1)/dV_airborne_bar(ind_d10);
    wt_above = 1-wt_below;
    d10_airborne_site(i) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    
    %get d50
    ind_d50 = find(CV>=0.5, 1);
    d_below = d_airborne_lower(ind_d50);
    d_above = d_airborne_upper(ind_d50);
    wt_below = (CV(ind_d50)-0.5)/dV_airborne_bar(ind_d50);
    wt_above = 1-wt_below;
    d50_airborne_site(i) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    
    %get d90
    ind_d90 = find(CV>=0.9, 1);
    d_below = d_airborne_lower(ind_d90);
    d_above = d_airborne_upper(ind_d90);
    wt_below = (CV(ind_d90)-0.9)/dV_airborne_bar(ind_d90);
    wt_above = 1-wt_below;
    d90_airborne_site(i) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    
    
    %% BSNE - GET DAILY MEAN
    dates_samples = [GrainSize_BSNE(ind_airborne).Date];
    dates_unique = unique(dates_samples);
    N_dates = length(dates_unique);
    dV_airborne_date = zeros(N_dates,N_d);
    dVdlogd_airborne_date = zeros(N_dates,N_d);
    for j = 1:N_dates
        ind_date = find(dates_samples==dates_unique(j));
        dV_airborne_date(j,:) = mean(dV_airborne(ind_date,:),1);
        dVdlogd_airborne_date(j,:) = mean(dVdlogd_airborne(ind_date,:),1);
    end

    %get daily d10, d50, d90
    CV = cumsum(dV_airborne_date')'; %cumulative distribution by volume
    dates_airborne{i} = dates_unique';
    d10_airborne_date{i} = zeros(N_dates,1);
    d50_airborne_date{i} = zeros(N_dates,1);
    d90_airborne_date{i} = zeros(N_dates,1);
    
    for j = 1:N_dates
        %get d10
        ind_d10 = find(CV(j,:)>=0.1, 1);
        d_below = d_airborne_lower(ind_d10);
        d_above = d_airborne_upper(ind_d10);
        wt_below = (CV(j,ind_d10)-0.1)/dV_airborne_date(j,ind_d10);
        wt_above = 1-wt_below;
        d10_airborne_date{i}(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

        %get d50
        ind_d50 = find(CV(j,:)>=0.5, 1);
        d_below = d_airborne_lower(ind_d50);
        d_above = d_airborne_upper(ind_d50);
        wt_below = (CV(j,ind_d50)-0.5)/dV_airborne_date(j,ind_d50);
        wt_above = 1-wt_below;
        d50_airborne_date{i}(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

        %get d90
        ind_d90 = find(CV(j,:)>=0.9, 1);
        d_below = d_airborne_lower(ind_d90);
        d_above = d_airborne_upper(ind_d90);
        wt_below = (CV(j,ind_d90)-0.9)/dV_airborne_date(j,ind_d90);
        wt_above = 1-wt_below;
        d90_airborne_date{i}(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));
    end
    
    
    %% BSNE - GET CLUSTER MEAN FOR OCEANO
    if i == 3;
        d_airborne_cluster = d_airborne;
        dlogd_airborne_cluster = dlogd_airborne;
        dV_airborne_cluster = zeros(N_clusters,N_d);
        dVdlogd_airborne_cluster = zeros(N_clusters,N_d);
        for j = 1:N_clusters
            ind_cluster = intersect(find(dates_samples>=cluster_StartDate(j)),find(dates_samples<=cluster_EndDate(j)));
            dV_airborne_cluster(j,:) = mean(dV_airborne(ind_cluster,:),1);
            dVdlogd_airborne_cluster(j,:) = mean(dVdlogd_airborne(ind_cluster,:),1);
        end

        %get daily d10, d50, d90
        CV = cumsum(dV_airborne_cluster')'; %cumulative distribution by volume
        d10_airborne_cluster = zeros(N_clusters,1);
        d50_airborne_cluster = zeros(N_clusters,1);
        d90_airborne_cluster = zeros(N_clusters,1);

        for j = 1:N_clusters
            %get d10
            ind_d10 = find(CV(j,:)>=0.1, 1);
            d_below = d_airborne_lower(ind_d10);
            d_above = d_airborne_upper(ind_d10);
            wt_below = (CV(j,ind_d10)-0.1)/dV_airborne_cluster(j,ind_d10);
            wt_above = 1-wt_below;
            d10_airborne_cluster(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

            %get d50
            ind_d50 = find(CV(j,:)>=0.5, 1);
            d_below = d_airborne_lower(ind_d50);
            d_above = d_airborne_upper(ind_d50);
            wt_below = (CV(j,ind_d50)-0.5)/dV_airborne_cluster(j,ind_d50);
            wt_above = 1-wt_below;
            d50_airborne_cluster(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));

            %get d90
            ind_d90 = find(CV(j,:)>=0.9, 1);
            d_below = d_airborne_lower(ind_d90);
            d_above = d_airborne_upper(ind_d90);
            wt_below = (CV(j,ind_d90)-0.9)/dV_airborne_cluster(j,ind_d90);
            wt_above = 1-wt_below;
            d90_airborne_cluster(j) = exp(wt_below*log(d_below)+wt_above*log(d_above));
        end
    end
    
end

%plot just surface GSD for each site
figure(1); clf; hold on;
for i = 1:N_Sites
    plot(d_surface_site{i},dVdlogd_surface_site{i},'LineWidth',LineWidth_Surface);
end
set(gca,'xscale','log','yscale','log');
legend(Sites,'Location','NorthWest');
ylim([1e-2 2]);
xlim([6e-2 2]);
ax = gca;
ax.XTick = [0.06:0.01:0.1, 0.2:0.1:1, 2];
ax.XTickLabel = {'0.06','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
ax.YTick = [0.01:0.01:0.1, 0.2:0.1:1, 2];
ax.YTickLabel = {'0.01','0.02','0.03','','0.05','','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
xlabel('Particle diameter, d (mm)');
ylabel('Norm. surface size distr., dV/dln(d)');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'MeanGSD_Surface.png'],'-dpng');


%plot mean GSD for each site
figure(2); clf;
subplot(2,1,1); hold on;
for i = 1:N_Sites
    plot(d_surface_site{i},dVdlogd_surface_site{i},'LineWidth',LineWidth_Surface);
end
set(gca,'xscale','log','yscale','log');
legend(Sites,'Location','NorthWest');
ylim([1e-2 2]);
xlim([6e-2 2]);
ax = gca;
ax.XTick = [0.06:0.01:0.1, 0.2:0.1:1, 2];
ax.XTickLabel = {'0.06','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
ax.YTick = [0.01:0.01:0.1, 0.2:0.1:1, 2];
ax.YTickLabel = {'0.01','0.02','0.03','','0.05','','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
title('Surface');
xlabel('Particle size, d (mm)');
ylabel('Norm. size distr., dV/dln(d)');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

subplot(2,1,2); hold on;
for i = 1:N_Sites
    plot(d_airborne_site{i},dVdlogd_airborne_site{i},'LineWidth',LineWidth_Surface);
end
set(gca,'xscale','log','yscale','log');
ylim([1e-2 2]);
xlim([6e-2 2]);
ax = gca;
ax.XTick = [0.06:0.01:0.1, 0.2:0.1:1, 2];
ax.XTickLabel = {'0.06','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
ax.YTick = [0.01:0.01:0.1, 0.2:0.1:1, 2];
ax.YTickLabel = {'0.01','0.02','0.03','','0.05','','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
title('Airborne');
xlabel('Particle size, d (mm)');
ylabel('Norm. size distr., dV/dln(d)');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

%print plot for draft
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
print([folder_Plots,'MeanGSD.png'],'-dpng');


%plot cluster gsd for Oceano
figure(3); clf;
subplot(2,1,1); hold on;
for i = 1:N_clusters
    plot(d_surface_cluster,dVdlogd_surface_cluster(i,:),'LineWidth',LineWidth_Surface);
end
set(gca,'xscale','log','yscale','log');
legend(cluster_names,'Location','NorthWest');
ylim([1e-2 2]);
xlim([6e-2 2]);
ax = gca;
ax.XTick = [0.06:0.01:0.1, 0.2:0.1:1, 2];
ax.XTickLabel = {'0.06','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
ax.YTick = [0.01:0.01:0.1, 0.2:0.1:1, 2];
ax.YTickLabel = {'0.01','0.02','0.03','','0.05','','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
title('Surface');
xlabel('Particle size, d (mm)');
ylabel('Norm. size distr., dV/dln(d)');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

subplot(2,1,2); hold on;
for i = 1:N_clusters
    plot(d_airborne_cluster,dVdlogd_airborne_cluster(i,:),'LineWidth',LineWidth_Surface);
end
set(gca,'xscale','log','yscale','log');
legend(cluster_names,'Location','NorthWest');
ylim([1e-2 2]);
xlim([6e-2 2]);
ax = gca;
ax.XTick = [0.06:0.01:0.1, 0.2:0.1:1, 2];
ax.XTickLabel = {'0.06','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
ax.YTick = [0.01:0.01:0.1, 0.2:0.1:1, 2];
ax.YTickLabel = {'0.01','0.02','0.03','','0.05','','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
title('Airborne');
xlabel('Particle size, d (mm)');
ylabel('Norm. size distr., dV/dln(d)');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

%print plot for draft
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
print([folder_Plots,'DateClusterGSD.png'],'-dpng');


%% SAVE DATA
save(SaveData_Path,'Sites','*site','*date','*cluster','dates*','cluster*','d10*','d50*','d90*'); 