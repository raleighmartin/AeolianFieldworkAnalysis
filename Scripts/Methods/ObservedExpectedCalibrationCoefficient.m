%% initialize
clearvars;

%% physical parameters
rho_s = 2650*1e3; %sediment density, g/m^3
A_w = 30*0.6*(1e-3)^2; %area of Wenglor, m^2
d_min_detection = 0; %minimum detected grain size, mm
max_dz = 0.1; %maximum distance between observed and predicted Cqn heights, m

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/BSNE/'; %folder for mean grain size data
folder_SaltationData = '../../AnalysisData/Windowing/'; %folder for saltation flux data
MeanGrainSizeData_Path = strcat(folder_AnalysisData,'MeanGrainSize'); %path for loading mean grain size data
SaltationFluxData_Path = strcat(folder_SaltationData,'DataWindowCalcs_30min_Restricted'); %path for loading saltation data
folder_Plots = '../../PlotOutput/Methods/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% plotting information
PlotFont = 10;
LineWidth_Plot = 1;
Marker_Site = {'s','d','o'};
Color_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

%% load grain size / BSNE data
GrainSizeBSNE_all = cell(N_Sites,1); %cell array containing grain size data for all sites
FluxBSNE_all = cell(N_Sites,1); %cell array containing BSNE flux data for all sites
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

%load Wenglor saltation flux data
load(SaltationFluxData_Path);

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize information for BSNE intervals
z_Cqn_pred = cell(N_Sites,1); %heights for Cqn prediction
z_Cqn_obs = cell(N_Sites,1); %heights for observed Cqn
Cqn_pred = cell(N_Sites,1); %expected calibration coefficient for BSNE interval
Cqn_obs = cell(N_Sites,1); %observed grain sizes for BSNE interval

%% go through each site
for i = 1:N_Sites
    
    %get data for site
    GrainSize_BSNE = GrainSizeBSNE_all{i};
    Flux_BSNE = FluxBSNE_all{i};
    N_FluxBSNE = length(Flux_BSNE);
    
    %get size bins from first sample
    d_grainsize = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_mid_mm];
    N_d = length(d_grainsize); %number of grain size bins
    
    %get indices of grain size distribution within Wenglor sensitivity range
    ind_sensitivity = find(d_grainsize>=d_min_detection);
    
    %initialize matrices for calibration coefficients
    z_Cqn_pred{i} = [];
    z_Cqn_obs{i} = [];
    Cqn_pred{i} = []; %expected calibration coefficient
    Cqn_obs{i} = []; %observed calibration coefficient
    
    %go through each BSNE interval to get calibration coefficients
    for j = 1:N_FluxBSNE
        
        %extract information for BSNEs
        z_BSNE = Flux_BSNE(j).z.z; %heights of BSNEs
        Name_BSNE = Flux_BSNE(j).name;
        StartTime_BSNE = Flux_BSNE(j).StartTime; %get start time for BSNE interval
        EndTime_BSNE = Flux_BSNE(j).EndTime; %get start time for BSNE interval
        
        %get indices of GrainSize_BSNE time intervals containing BSNE time interval
        ind_GrainSize = find([GrainSize_BSNE.StartTime]<=StartTime_BSNE & [GrainSize_BSNE.EndTime]>StartTime_BSNE);

        %get index of first Wenglor interval within BSNE time interval, then get Wenglor heights and calibration coefficients
        ind_Wenglor = find(StartTimes_all{i}>=StartTime_BSNE & StartTimes_all{i}<=EndTime_BSNE,1);
        z_Wenglor = zW_all{i}{ind_Wenglor};
        Cqn_Wenglor = Cqnbar_all{i}{ind_Wenglor};
        
        %initialize matrix of grain size samples
        N_grainsize_interval = length(ind_GrainSize);
        dV_airborne_interval = zeros(N_grainsize_interval,N_d);
        
        %get each airborne size distribution
        for k = 1:N_grainsize_interval
            dV_airborne_interval(k,:) = [GrainSize_BSNE(ind_GrainSize(k)).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
            
            %get index for BSNE
            ind_BSNE = find(strcmp(Name_BSNE,GrainSize_BSNE(ind_GrainSize(k)).NameBSNE));
                       
            if ~isempty(ind_BSNE)
                %calculate predicted Cqn, add to list
                z_pred = z_BSNE(ind_BSNE);
                z_Cqn_pred{i} = [z_Cqn_pred{i}; z_pred]; %add height of Cqn prediction to list
                d_bar_interval = sum(dV_airborne_interval(k,ind_sensitivity).*d_grainsize(ind_sensitivity))/sum(dV_airborne_interval(k,ind_sensitivity)); %get mean grain size for BSNE interval (mm)
                Cqn_pred_interval = (pi/6)*(1e-3*d_bar_interval).^3*(rho_s)/A_w; %get expected calibration coefficient
                Cqn_pred{i} = [Cqn_pred{i}; Cqn_pred_interval]; %add expected calibration coefficient to list

                %find closest Wenglor z, get observed Cqn, add to list
                z_diff = abs(z_Wenglor - z_pred);
                ind_Wenglor_closest = find(z_diff == min(z_diff));
                z_Cqn_obs{i} = [z_Cqn_obs{i}; z_Wenglor(1)]; %add Wenglor height to list for observed Cqn
                Cqn_obs{i} = [Cqn_obs{i}; mean(Cqn_Wenglor(ind_Wenglor_closest))]; %add Wenglor calibration coefficient to list of observed Cqn
                %take mean value here if there are more than one
            end           
        end
    end
end

%% perform fit to calibration factor versus counts rate and wind direction
%initialize values
nbar_combined = cell(N_Sites,1);
Cqnbar_combined = cell(N_Sites,1);
theta_combined = cell(N_Sites,1);
std_theta_combined = cell(N_Sites,1);
Cqnbar_pred_nbar = cell(N_Sites,1);
Cqnbar_pred_theta = cell(N_Sites,1);
Cqnbar_pred_std_theta = cell(N_Sites,1);
sigma_Cqnbar_combined = cell(N_Sites,1);
slope_Cqnbar_nbar = zeros(N_Sites,1);
slope_Cqnbar_theta = zeros(N_Sites,1);
slope_Cqnbar_std_theta = zeros(N_Sites,1);
sigma_slope_Cqnbar_nbar = zeros(N_Sites,1);
sigma_slope_Cqnbar_theta = zeros(N_Sites,1);
sigma_slope_Cqnbar_std_theta = zeros(N_Sites,1);
for i = 1:N_Sites
    %combine values
    ind_flux = find(Q_all{i}>0);
    for j = 1:length(ind_flux)
        nbar_combined{i} = [nbar_combined{i} nbar_all{i}{ind_flux(j)}];
        Cqnbar_combined{i} = [Cqnbar_combined{i} Cqnbar_all{i}{ind_flux(j)}];
        sigma_Cqnbar_combined{i} = [sigma_Cqnbar_combined{i} sigma_Cqnbar_all{i}{ind_flux(j)}];
        theta_combined{i} = [theta_combined{i} abs(theta_all{i}(ind_flux(j)))*ones(size(nbar_all{i}{ind_flux(j)}))];
        %std_theta_combined{i} = [std_theta_combined{i} std_theta_all{i}(ind_flux(j))*ones(size(nbar_all{i}{ind_flux(j)}))];
        std_theta_combined{i} = [std_theta_combined{i} sigma_theta_all{i}(ind_flux(j))*ones(size(nbar_all{i}{ind_flux(j)}))];
    end
    %remove zero values
    ind_nonzero = find(nbar_combined{i}~=0);
    nbar_combined{i} = nbar_combined{i}(ind_nonzero);
    Cqnbar_combined{i} = Cqnbar_combined{i}(ind_nonzero);
    sigma_Cqnbar_combined{i} = sigma_Cqnbar_combined{i}(ind_nonzero);  
    theta_combined{i} = theta_combined{i}(ind_nonzero);
    std_theta_combined{i} = std_theta_combined{i}(ind_nonzero);
    
    %perform linear fit to log of nbar versus log of Cqnbar
    [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = ...
        linearfit(log10(nbar_combined{i}), log10(Cqnbar_combined{i}), log10(sigma_Cqnbar_combined{i}));
    Cqnbar_pred_nbar{i} = 10.^yfit;
    slope_Cqnbar_nbar(i) = b;
    sigma_slope_Cqnbar_nbar(i) = sigma_b;
    
    %perform linear fit to log of Cqnbar versus theta
    [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = ...
        linearfit(theta_combined{i}, log10(Cqnbar_combined{i}), log10(sigma_Cqnbar_combined{i}));
    Cqnbar_pred_theta{i} = 10.^yfit;
    slope_Cqnbar_theta(i) = b;
    sigma_slope_Cqnbar_theta(i) = sigma_b;
    
    %perform linear fit to log of Cqnbar versus std_theta
    [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab, da_dy, db_dy] = ...
        linearfit(std_theta_combined{i}, log10(Cqnbar_combined{i}), log10(sigma_Cqnbar_combined{i}));
    Cqnbar_pred_std_theta{i} = 10.^yfit;
    slope_Cqnbar_std_theta(i) = b;
    sigma_slope_Cqnbar_std_theta(i) = sigma_b;    
end

%%%%%%%%%
%% PLOT %
%%%%%%%%%

figure(1); clf; 

%% observed calibration factor
subplot('Position',[0.07 0.13 0.235 0.84]); hold on;
for i = 1:N_Sites
    for j = 1:length(zW_all{i})
        if ~isnan(zW_all{i}{j})
            plot(zW_all{i}{j},Cqnbar_all{i}{j},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
        end
    end
end
xlabel('HF sensor height, $$z_{HF,i}$$ (m)','Interpreter','Latex')
ylabel('Observed calibration factor, $$C_{qn,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
text(0.012,80,'(a)','FontSize',PlotFont);
xlim([1e-2, 1]); 
ylim([1e-1 1e2]);
set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont);

%% expected calibration factor
subplot('Position',[0.39 0.13 0.235 0.84]); hold on;
for i = 1:N_Sites
    plot(z_Cqn_pred{i},Cqn_pred{i},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
end
xlabel('LF trap height, $$z_{LF,i}$$ (m)','Interpreter','Latex')
ylabel('Predicted calibration factor, $$C_{qn,pred,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
text(0.012,80,'(b)','FontSize',PlotFont);
xlim([1e-2, 1]); 
ylim([1e-1 1e2]);
set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont)

%% comparison of calibration factors
subplot('Position',[0.715 0.13 0.235 0.84]); hold on;

%plot 1:1 line
plot([1e-1 1e2],[1e-1 1e2],'--')

for i = 1:N_Sites
    ind_dz = find(abs(z_Cqn_pred{i} - z_Cqn_obs{i})<= max_dz); %determine which values have dz < max_dz
    plot(Cqn_obs{i}(ind_dz),Cqn_pred{i}(ind_dz),Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
end

%format plot
xlim([1e-1 1e2])
ylim([1e-1 1e2]);
text(1.2e-1,80,'(c)','FontSize',PlotFont);
xlabel('Observed calibration factor, $$C_{qn,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
ylabel('Predicted calibration factor, $$C_{qn,pred,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont);

%create legend
legend_items = {'1-1 line',SiteNames{1},SiteNames{2},SiteNames{3}};
legend(legend_items,'Location','NorthEast');

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8.5 3.5],'PaperPosition',[0 0 8.5 3.5],'PaperPositionMode','Manual');
print([folder_Plots,'CalibrationObservedPredicted.png'],'-dpng');
print([folder_Plots,'CalibrationObservedPredicted.tif'],'-dtiff');

%% plot observed calibration factor versus number counts rate
figure(2); clf; hold on;

%plot values
for i = 1:N_Sites
    plot(nbar_combined{i},Cqnbar_combined{i},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
end

%plot fit lines
for i = 1:N_Sites
    plot(nbar_combined{i},Cqnbar_pred_nbar{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot*2);
end

xlabel('HF sensor counts rate over calibration interval, $$n_{cal,HF,i}$$ (s $$^{-1}$$)','Interpreter','Latex')
ylabel('Observed calibration factor, $$C_{qn,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
legend(SiteNames,'Location','SouthWest');
set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont);
set(gca,'LooseInset',get(gca,'TightInset'));

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 4.5],'PaperPosition',[0 0 6 4.5],'PaperPositionMode','Manual');
print([folder_Plots,'Calibration_CountsRate.png'],'-dpng');
print([folder_Plots,'Calibration_CountsRate.tif'],'-dtiff');

%% plot observed calibration factor versus wind direction

figure(3); clf; hold on;

%% subplot versus mean wind direction
subplot('Position',[0.08 0.12 0.42 0.86]); hold on;

%plot values
for i = 1:N_Sites
    plot(theta_combined{i},Cqnbar_combined{i},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
end

%plot fit lines
for i = 1:N_Sites
    plot(theta_combined{i},Cqnbar_pred_theta{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot*2);
end

xlabel('Mean absolute wind direction, $$|\theta|$$ ($$^{\circ}$$)','Interpreter','Latex');
ylabel('Observed calibration factor, $$C_{qn,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
set(gca,'xscale','linear','yscale','log','box','on','FontSize',PlotFont);
set(gca,'LooseInset',get(gca,'TightInset'));
ylims = ylim;
text(1.7,ylims(1)+0.85*range(ylims),'(a)','FontSize',PlotFont);

%% subplot versus std of wind variation
subplot('Position',[0.56 0.12 0.42 0.86]); hold on;

%plot values
for i = 1:N_Sites
    plot(std_theta_combined{i},Cqnbar_combined{i},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
end

%plot fit lines
for i = 1:N_Sites
    plot(std_theta_combined{i},Cqnbar_pred_std_theta{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot*2);
end

xlabel('Std. dev. of 2-second wind dir., $$\sigma_{\theta,2s}$$ ($$^{\circ}$$)','Interpreter','Latex');
legend(SiteNames,'Location','NorthEast');
set(gca,'xscale','linear','yscale','log','box','on','FontSize',PlotFont);
set(gca,'LooseInset',get(gca,'TightInset'));
ylims = ylim;
text(4.5,ylims(1)+0.85*range(ylims),'(b)','FontSize',PlotFont);


%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 4],'PaperPosition',[0 0 8 4],'PaperPositionMode','Manual');
print([folder_Plots,'Calibration_WindDirection.png'],'-dpng');
print([folder_Plots,'Calibration_WindDirection.tif'],'-dtiff');

%% plot observed calibration factor versus saltation flux
figure(5); clf; hold on;
for i = 1:N_Sites
    plot(Q_all{i}(1)*ones(size(Cqnbar_all{i}{1})),Cqnbar_all{i}{1},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
end

for i = 1:N_Sites
    for j = 1:length(Cqnbar_all{i})
        if ~isnan(Cqnbar_all{i}{j})
            plot(Q_all{i}(j)*ones(size(Cqnbar_all{i}{j})),Cqnbar_all{i}{j},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
        end
    end
end
xlabel('Observed saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','Latex');
ylabel('Observed calibration factor, $$C_{qn,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');
legend(SiteNames,'Location','SouthWest');
set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont);
set(gca,'LooseInset',get(gca,'TightInset'));

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 4.5],'PaperPosition',[0 0 6 4.5],'PaperPositionMode','Manual');
print([folder_Plots,'Calibration_SaltationFlux.png'],'-dpng');
print([folder_Plots,'Calibration_SaltationFlux.tif'],'-dtiff');