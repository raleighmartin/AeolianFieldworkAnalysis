%% initialize
clearvars;

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/GrainSize/'; %folder for outputs
MeanGrainSizeData_Path = strcat(folder_AnalysisData,'MeanGrainSize'); %path for loading mean grain size data
folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% plotting information
PlotFont = 14;
LineWidth_Surface = 1;

%% load grain size / BSNE data
GrainSizeBSNE_all = cell(N_Sites,1);
FluxBSNE_all = cell(N_Sites,1);
for i = 1:N_Sites
    %load grain size data
    GrainSizeData_Path = strcat(folder_ProcessedData,'GrainSize_',Sites{i});
    load(GrainSizeData_Path);
    GrainSizeBSNE_all{i} = GrainSize_BSNE;

    %load BSNE flux data
    BSNEData_Path = strcat(folder_ProcessedData,'FluxBSNE_',Sites{i});
    load(BSNEData_Path);
    FluxBSNE_all{i} = FluxBSNE;
end

%load mean grain size data
load(MeanGrainSizeData_Path);

%% initialize information for BSNE intervals
d_airborne_interval = cell(N_Sites,1);
dlogd_airborne_interval = cell(N_Sites);
Q_interval = cell(N_Sites,1);
dV_airborne_interval_bar = cell(N_Sites,1);
dVdlogd_airborne_interval_bar = cell(N_Sites,1);
d10_airborne_interval_bar = cell(N_Sites,1);
d50_airborne_interval_bar = cell(N_Sites,1);
d90_airborne_interval_bar = cell(N_Sites,1);

%% go through each site
for i = 1:N_Sites
    
    %get data for site
    GrainSize_BSNE = GrainSizeBSNE_all{i};
    Flux_BSNE = FluxBSNE_all{i};
    N_FluxBSNE = length(Flux_BSNE);
    
    %get size bins from first sample
    d_airborne_interval{i} = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_mid_mm];
    d_airborne_lower = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_lower_mm];
    d_airborne_upper = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_upper_mm];
    dlogd_airborne_interval{i} = log(d_airborne_upper) - log(d_airborne_lower);
    N_d = length(d_airborne_interval{i});
    
    %initialize matrices
    Q_interval{i} = zeros(N_FluxBSNE,1); %total fluxes for grain size samples
    dV_airborne_interval_bar{i} = zeros(N_FluxBSNE,N_d); %size distributions for grain size samples
    dVdlogd_airborne_interval_bar{i} = zeros(N_FluxBSNE,N_d); %size distributions for grain size samples
    d10_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d10 for grain size samples
    d50_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d50 for grain size samples
    d90_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d90 for grain size samples
    
    %get indices of airborne samples
    for j = 1:N_FluxBSNE
        Q_interval{i}(j) = Flux_BSNE(j).Q.Q;
        qz_BSNE = Flux_BSNE(j).qz.qz;
        Name_BSNE = Flux_BSNE(j).name;
        StartTime_BSNE = Flux_BSNE(j).StartTime; %get start time for BSNE interval
        ind_GrainSize = find([GrainSize_BSNE.StartTime]<=StartTime_BSNE & [GrainSize_BSNE.EndTime]>StartTime_BSNE); %get indices within GrainSize_BSNE corresponding to time interval
    
        %initialize matrix of airborne samples
        N_airborne_interval = length(ind_GrainSize);
        dV_airborne_interval = zeros(N_airborne_interval,N_d);
        dVdlogd_airborne_interval = zeros(N_airborne_interval,N_d);
        qz_interval = zeros(N_airborne_interval,1); %partial fluxes for grain size samples
        
        %get each airborne size distribution
        for k = 1:N_airborne_interval
            dV_airborne_interval(k,:) = [GrainSize_BSNE(ind_GrainSize(k)).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
            dVdlogd_airborne_interval(k,:) = dV_airborne_interval(k,:)./dlogd_airborne_interval{i}; %by normalized volume
            try
                qz_interval(k) = qz_BSNE(strcmp(Name_BSNE,GrainSize_BSNE(ind_GrainSize(k)).NameBSNE)); %get flux associated with BSNE
            catch
                qz_interval(k) = 0;
            end
        end
        
        %get mean of airborne sizes
        dV_airborne_interval_bar{i}(j,:) = sum(dV_airborne_interval.*qz_interval)/sum(qz_interval);
        dVdlogd_airborne_interval_bar{i}(j,:) = sum(dVdlogd_airborne_interval.*qz_interval)/sum(qz_interval);
                
        %compute benchmark grain sizes
        if(~isnan(dV_airborne_interval_bar{i}(j,1)))
            [d10, d50, d90] = ReferenceGrainSizes(dV_airborne_interval_bar{i}(j,:), d_airborne_lower, d_airborne_upper);
            d10_airborne_interval_bar{i}(j) = d10;
            d50_airborne_interval_bar{i}(j) = d50;
            d90_airborne_interval_bar{i}(j) = d90;
        end
    end
end

figure(1); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    h90 = plot(Q_interval{i},d90_airborne_interval_bar{i},'^');
    h50 = plot(Q_interval{i},d50_airborne_interval_bar{i},'o');
    h10 = plot(Q_interval{i},d10_airborne_interval_bar{i},'v');
    c90 = get(h90,'Color');
    c50 = get(h50,'Color');
    c10 = get(h10,'Color');
    plot([0 max(Q_interval{i})],d90_surface_site(i)*[1 1],'Color',c90);
    plot([0 max(Q_interval{i})],d50_surface_site(i)*[1 1],'Color',c50);
    plot([0 max(Q_interval{i})],d10_surface_site(i)*[1 1],'Color',c10);
    ylim([0 0.9]);
    xlabel('saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$','Interpreter','Latex')
    ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
    title(Sites{i});
end
legend('d_{90}','d_{50}','d_{10}');

%print plot for draft
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'FluxGrainSize.png'],'-dpng');