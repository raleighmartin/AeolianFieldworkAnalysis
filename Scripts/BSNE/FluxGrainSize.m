%% initialize
clearvars;

%% physical parameters
rho_s = 2650*1e3; %sediment density, g/m^3
A_w = 30*0.6*(1e-3)^2; %area of Wenglor, m^2

%% partial flux grain-size bin values
%d_bin_lower = [0.1, 0.18, 0.32, 0.56];
%d_bin_upper = [0.18, 0.32, 0.56, 1];
d_bin_lower = [0.1, 0.18, 0.32];
d_bin_upper = [0.18, 0.32, 1];
d_bin_mid = geomean([d_bin_lower; d_bin_upper]);
N_d_bin = length(d_bin_lower);

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/GrainSize/'; %folder for mean grain size data
folder_SaltationData = '../../AnalysisData/Windowing/'; %folder for saltation flux data
MeanGrainSizeData_Path = strcat(folder_AnalysisData,'MeanGrainSize'); %path for loading mean grain size data
SaltationFluxData_Path = strcat(folder_SaltationData,'DataWindowCalcs_30min_Restricted'); %path for loading saltation data
folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots
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

%load saltation flux data
load(SaltationFluxData_Path); %path for loading saltation data

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% compute fraction of surface grains in each bin
f_d_surface = cell(N_Sites,1);
for i = 1:N_Sites
    f_d_surface{i} = zeros(1,N_d_bin);
    for j = 1:N_d_bin
        ind_full_bin = find(d_surface_lower_site{i}>=d_bin_lower(j) & d_surface_upper_site{i}<=d_bin_upper(j)); %grain size bins completely in bin
        ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
        wt_bottom_bin = (log(d_surface_upper_site{i}(ind_bottom_bin))-log(d_bin_lower(j)))/...
            (log(d_surface_upper_site{i}(ind_bottom_bin))-log(d_surface_lower_site{i}(ind_bottom_bin))); %fraction of partial bin for this bin
        ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
        wt_top_bin = (log(d_bin_upper(j))-log(d_surface_lower_site{i}(ind_top_bin)))/...
            (log(d_surface_upper_site{i}(ind_top_bin))-log(d_surface_lower_site{i}(ind_top_bin))); %fraction of partial bin for this bin
        f_d_surface{i}(j) = sum(dV_surface_site{i}(ind_full_bin))+...
            wt_bottom_bin*dV_surface_site{i}(ind_bottom_bin)+...
            wt_top_bin*dV_surface_site{i}(ind_top_bin); %compute fraction of grains by volume in this bin
    end
end

%% initialize information for BSNE intervals
d_airborne_interval = cell(N_Sites,1);
dlogd_airborne_interval = cell(N_Sites);
Q_interval = cell(N_Sites,1);
ust_interval = cell(N_Sites,1); %shear velocity for interval
tau_interval = cell(N_Sites,1); %shear stress for interval
dV_airborne_interval_bar = cell(N_Sites,1);
dVdlogd_airborne_interval_bar = cell(N_Sites,1);
d10_airborne_interval_bar = cell(N_Sites,1);
d50_airborne_interval_bar = cell(N_Sites,1);
d90_airborne_interval_bar = cell(N_Sites,1);      
f_d_airborne = cell(N_Sites,1); %compute fraction of airborne sample in each bin
f_d_airborne_surface_ratio = cell(N_Sites,1); %fraction of airborne sample in size bin versus surface sample

%information about grain size variation with height
d10_airborne_z = cell(N_Sites,1);
d50_airborne_z = cell(N_Sites,1);
d90_airborne_z = cell(N_Sites,1);
dV_airborne = cell(N_Sites,1);
z_airborne = cell(N_Sites,1);
z_airborne_ind = cell(N_Sites,1); %indices of intervals associated with z's
Cqn_airborne = cell(N_Sites,1); %expected calibration coefficient

% information about size-conditioned BSNE fluxes
q_profile_bin = cell(N_Sites,1); %size conditioned flux profile
Q_bin = cell(N_Sites,1); %size conditioned total flux

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
    d_airborne_mid = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_mid_mm];
    dlogd_airborne_interval{i} = log(d_airborne_upper) - log(d_airborne_lower);
    N_d = length(d_airborne_interval{i});
    
    %initialize matrices
    Q_interval{i} = zeros(N_FluxBSNE,1); %total fluxes for grain size samples
    ust_interval{i} = zeros(N_FluxBSNE,1); %shear velocities for grain size samples
    tau_interval{i} = zeros(N_FluxBSNE,1); %shear stress for interval
    dV_airborne_interval_bar{i} = zeros(N_FluxBSNE,N_d); %size distributions for grain size samples
    dVdlogd_airborne_interval_bar{i} = zeros(N_FluxBSNE,N_d); %size distributions for grain size samples
    d10_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d10 for grain size samples
    d50_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d50 for grain size samples
    d90_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d90 for grain size samples
    f_d_airborne{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %compute fraction of airborne sample in each bin
    f_d_airborne_surface_ratio{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %compute fraction of airborne sample in each bin versus surface fraction
    
    %information about grain size variation with height
    d10_airborne_z{i} = [];
    d50_airborne_z{i} = [];
    d90_airborne_z{i} = [];
    dV_airborne{i} = [];
    z_airborne{i} = [];
    z_airborne_ind{i} = []; %indices of intervals associated with z's
    Cqn_airborne{i} = []; %expected calibration coefficient
    
    %get indices of airborne samples
    for j = 1:N_FluxBSNE
        Q_interval{i}(j) = Flux_BSNE(j).Q.Q;
        z_BSNE = Flux_BSNE(j).z.z; %heights of BSNEs
        qz_BSNE = Flux_BSNE(j).qz.qz;
        Name_BSNE = Flux_BSNE(j).name;
        StartTime_BSNE = Flux_BSNE(j).StartTime; %get start time for BSNE interval
        EndTime_BSNE = Flux_BSNE(j).EndTime; %get start time for BSNE interval
        ind_GrainSize = find([GrainSize_BSNE.StartTime]<=StartTime_BSNE & [GrainSize_BSNE.EndTime]>StartTime_BSNE); %get indices within GrainSize_BSNE corresponding to time interval
    
        %get ust for interval
        ind_interval = find(StartTimes_all{i}>=StartTime_BSNE & EndTimes_all{i}<=EndTime_BSNE);
        ust_interval{i}(j)=mean(ustRe_all{i}(ind_interval)); %shear velocity for interval
        tau_interval{i}(j)=mean(tauRe_all{i}(ind_interval)); %shear stress for interval
        
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
                ind_BSNE = find(strcmp(Name_BSNE,GrainSize_BSNE(ind_GrainSize(k)).NameBSNE)); %get index for BSNE
                qz_interval(k) = qz_BSNE(ind_BSNE); %get flux associated with BSNE
            catch
                ind_BSNE = [];
                qz_interval(k) = 0;
            end
            
            %get reference sizes, add to list
            if ~isempty(ind_BSNE)
                [d10, d50, d90] = ReferenceGrainSizes(dV_airborne_interval(k,:), d_airborne_lower, d_airborne_upper);
                d10_airborne_z{i} = [d10_airborne_z{i}; d10];
                d50_airborne_z{i} = [d50_airborne_z{i}; d50];
                d90_airborne_z{i} = [d90_airborne_z{i}; d90];
                dV_airborne{i} = [dV_airborne{i}; dV_airborne_interval(k,:)];
                z_interval = z_BSNE(ind_BSNE); %get height associated with BSNE
                z_airborne{i} = [z_airborne{i}; z_interval];
                z_airborne_ind{i} = [z_airborne_ind{i}; j];
                d_bar = sum(dV_airborne_interval(k,:).*d_airborne_mid)/sum(dV_airborne_interval(k,:)); %get mean grain size, mm
                Cqn = (pi/6)*(1e-3*d_bar).^3*(rho_s)/A_w; %get expected calibration coefficient
                Cqn_airborne{i} = [Cqn_airborne{i}; Cqn]; %expected calibration coefficient
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
        
        %compute fraction of airborne sample in each bin
        for k = 1:N_d_bin
            ind_full_bin = find(d_airborne_lower>=d_bin_lower(k) & d_airborne_upper<=d_bin_upper(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_airborne_upper(ind_bottom_bin))-log(d_bin_lower(k)))/...
                (log(d_airborne_upper(ind_bottom_bin))-log(d_airborne_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper(k))-log(d_airborne_lower(ind_top_bin)))/...
                (log(d_airborne_upper(ind_top_bin))-log(d_airborne_lower(ind_top_bin))); %fraction of partial bin for this bin
            f_d_airborne{i}(j,k) = sum(dV_airborne_interval_bar{i}(j,ind_full_bin))+...
                wt_bottom_bin*dV_airborne_interval_bar{i}(j,ind_bottom_bin)+...
                wt_top_bin*dV_airborne_interval_bar{i}(j,ind_top_bin); %compute fraction of grains by volume in this bin
        end
        f_d_airborne_surface_ratio{i}(j,:) = f_d_airborne{i}(j,:)./f_d_surface{i};
    end

    %% SIZE-CONDITIONED BSNE FLUXES
    %information about size-conditioned BSNE fluxes
    q_profile_bin{i} = cell(N_FluxBSNE,N_d_bin); %size conditioned flux profile
    Q_bin{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %size conditioned total flux
    
    %go through each time interval
    for j = 1:N_FluxBSNE
        
        %grain size information for this interval
        ind_interval = find(z_airborne_ind{i} == j);
        z_interval = z_airborne{i}(ind_interval);
        dV_interval = dV_airborne{i}(ind_interval,:);
        
        %BSNE flux profile for this interval
        z_BSNE = Flux_BSNE(j).z.z;
        sigma_z_BSNE = Flux_BSNE(j).z.sigma_z;
        qz_BSNE = Flux_BSNE(j).qz.qz;
        sigma_qz_BSNE = Flux_BSNE(j).qz.sigma;
        
        %compute size conditioned flux
        for k = 1:N_d_bin
            ind_full_bin = find(d_airborne_lower>=d_bin_lower(k) & d_airborne_upper<=d_bin_upper(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_airborne_upper(ind_bottom_bin))-log(d_bin_lower(k)))/...
                (log(d_airborne_upper(ind_bottom_bin))-log(d_airborne_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper(k))-log(d_airborne_lower(ind_top_bin)))/...
                (log(d_airborne_upper(ind_top_bin))-log(d_airborne_lower(ind_top_bin))); %fraction of partial bin for this bin
            
            %initialize q_profile_bin for specific BSNE interval
            q_profile_bin{i}{j,k} = zeros(size(qz_BSNE));
            sigma_q_profile_bin = zeros(size(qz_BSNE));
            
            for l = 1:length(z_BSNE)
                dz = abs(z_BSNE(l)-z_interval);
                ind_closest = find(dz==min(dz));
                if ~isempty(ind_closest)
                    f_d = sum(dV_interval(ind_closest,ind_full_bin))+...
                        wt_bottom_bin*dV_interval(ind_closest,ind_bottom_bin)+...
                        wt_top_bin*dV_interval(ind_closest,ind_top_bin); %compute fraction of grains by volume in this bin
                    %q_profile_bin{i}{j,k}(l) = qz_BSNE(l)*f_d/f_d_surface{i}(k);
                    q_profile_bin{i}{j,k}(l) = qz_BSNE(l)*f_d;
                    sigma_q_profile_bin(l) = sigma_qz_BSNE(l)*f_d/f_d_surface{i}(k);
                end
            end
            
            %fit profile
            [~,~,Q,~,~,sigma_Q] = qz_profilefit(q_profile_bin{i}{j,k}, z_BSNE, sigma_q_profile_bin, sigma_z_BSNE);
            Q_bin{i}(j,k) = Q; %size conditioned total flux
        end
    end
end

% %%plot variation in reference grain sizes with saltation flux
% figure(1); clf;
% for i = 1:N_Sites
%     subplot(1,N_Sites,i); hold on;
%     h90 = plot(Q_interval{i},d90_airborne_interval_bar{i},'^');
%     h50 = plot(Q_interval{i},d50_airborne_interval_bar{i},'o');
%     h10 = plot(Q_interval{i},d10_airborne_interval_bar{i},'v');
%     c90 = get(h90,'Color');
%     c50 = get(h50,'Color');
%     c10 = get(h10,'Color');
%     plot([0 max(Q_interval{i})],d90_surface_site(i)*[1 1],'Color',c90);
%     plot([0 max(Q_interval{i})],d50_surface_site(i)*[1 1],'Color',c50);
%     plot([0 max(Q_interval{i})],d10_surface_site(i)*[1 1],'Color',c10);
%     ylim([0 0.9]);
%     xlabel('saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','Latex')
%     ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     title(Sites{i});
% end
% legend('d_{90}','d_{50}','d_{10}');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
% print([folder_Plots,'FluxGrainSize.png'],'-dpng');
% 
% %%plot variation in reference grain sizes with shear velocity
% figure(2); clf;
% for i = 1:N_Sites
%     subplot(1,N_Sites,i); hold on;
%     h90 = plot(ust_interval{i},d90_airborne_interval_bar{i},'^');
%     h50 = plot(ust_interval{i},d50_airborne_interval_bar{i},'o');
%     h10 = plot(ust_interval{i},d10_airborne_interval_bar{i},'v');
%     c90 = get(h90,'Color');
%     c50 = get(h50,'Color');
%     c10 = get(h10,'Color');
%     plot([0 max(ust_interval{i})],d90_surface_site(i)*[1 1],'Color',c90);
%     plot([0 max(ust_interval{i})],d50_surface_site(i)*[1 1],'Color',c50);
%     plot([0 max(ust_interval{i})],d10_surface_site(i)*[1 1],'Color',c10);
%     ylim([0 0.9]);
%     xlabel('shear velocity, $$u_{*}$$ (m s$$^{-1}$$)','Interpreter','Latex')
%     ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     title(Sites{i});
% end
% legend('d_{90}','d_{50}','d_{10}');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
% print([folder_Plots,'ShearVelocityGrainSize.png'],'-dpng');
% 
% %%plot variation in partial fluxes with total flux
% figure(3); clf;
% for k = 1:N_d_bin
%     subplot(1,N_d_bin,k); hold on;
%     for i = 1:N_Sites
%         plot(Q_interval{i},f_d_airborne_surface_ratio{i}(:,k),Marker_Site{i},'Color',Color_Site{i})
%         title(['d = ',num2str(d_bin_lower(k)),' - ',num2str(d_bin_upper(k)),' mm']);
%     end
%     xlabel('saltation flux, $$Q$$','interpreter','latex');
%     ylabel('fraction of size in air versus surface, $$f_{d,air}/f_{d,surf}$$','interpreter','latex');
%     set(gca,'YScale','Log','YTickLabel',{'0.01','0.1','1','10'},'Box','On');
%     ylim([1e-2 1e1]);
% end
% legend(Sites,'Location','NorthEast');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 4]);
% print([folder_Plots,'Flux_GrainSizeRelativeFlux.png'],'-dpng');
% 
% %%plot variation in partial fluxes with shear velocity
% figure(4); clf;
% for k = 1:N_d_bin
%     subplot(1,N_d_bin,k); hold on;
%     for i = 1:N_Sites
%         plot(ust_interval{i},f_d_airborne_surface_ratio{i}(:,k),Marker_Site{i},'Color',Color_Site{i})
%         title(['d = ',num2str(d_bin_lower(k)),' - ',num2str(d_bin_upper(k)),' mm']);
%     end
%     xlabel('shear velocity, $$u_{*}$$ (m s$$^{-1}$$)','Interpreter','Latex')
%     ylabel('fraction of size in air versus surface, $$f_{d,air}/f_{d,surf}$$','interpreter','latex');
%     set(gca,'YScale','Log','YTickLabel',{'0.01','0.1','1','10'},'Box','On');
%     ylim([1e-2 1e1]);
% end
% legend(Sites,'Location','NorthEast');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 4]);
% print([folder_Plots,'ShearVelocity_GrainSizeRelativeFlux.png'],'-dpng');
% 
% %%plot variation in airborne grain sizes with height
% figure(5); clf;
% for i = 1:N_Sites
%     %subplot(1,N_Sites,i); hold on;
%     subplot('Position',[-0.26+0.325*i 0.13 0.27 0.82]); hold on;
%     h90 = plot(z_airborne{i},d90_airborne_z{i},'^');
%     h50 = plot(z_airborne{i},d50_airborne_z{i},'o');
%     h10 = plot(z_airborne{i},d10_airborne_z{i},'v');
%     c90 = get(h90,'Color');
%     c50 = get(h50,'Color');
%     c10 = get(h10,'Color');
%     s90 = plot([0 max(z_airborne{i})],d90_surface_site(i)*[1 1],'Color',c90);
%     s50 = plot([0 max(z_airborne{i})],d50_surface_site(i)*[1 1],'Color',c50);
%     s10 = plot([0 max(z_airborne{i})],d10_surface_site(i)*[1 1],'Color',c10);
%     ylim([0 0.9]);
%     xlabel('height, $$z$$ (m)','Interpreter','Latex')
%     if i==1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     elseif i==2
%         legend([h90,h50,h10],'airborne d_{90}','airborne d_{50}','airborne d_{10}','Location','NorthEast');
%     elseif i==3
%         legend([s90,s50,s10],'surface d_{90}','surface d_{50}','surface d_{10}','Location','NorthEast');
%     end
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
%     title(Sites{i});
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9 4],'PaperPosition',[0 0 9 4],'PaperPositionMode','Manual');
% print([folder_Plots,'HeightGrainSize.png'],'-dpng');
% 
% %%plot variation in expected calibration factor with height
% figure(6); clf; 
% subplot('Position',[0.08 0.115 0.4 0.86]); hold on;
% for i = 1:N_Sites
%     for j = 1:length(zW_all{i})
%         if ~isnan(zW_all{i}{j})
%             plot(zW_all{i}{j},Cqnbar_all{i}{j},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
%         end
%     end
% end
% xlabel('Wenglor height, $$z$$ (m)','Interpreter','Latex')
% ylabel('Observed calibration factor, $$C_{qn}$$ (g m$$^{-2}$$)','Interpreter','Latex');
% text(0.012,80,'(a)','FontSize',PlotFont);
% xlim([1e-2, 1]); 
% ylim([1e-1 1e2]);
% set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont);
% 
% subplot('Position',[0.58 0.115 0.4 0.86]); hold on;
% for i = 1:N_Sites
%     plot(z_airborne{i},Cqn_airborne{i},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
% end
% xlabel('Height of grain size measurement, $$z$$ (m)','Interpreter','Latex')
% ylabel('Expected calibration factor, $$C_{qn,pred}$$ (g m$$^{-2}$$)','Interpreter','Latex');
% text(0.012,80,'(b)','FontSize',PlotFont);
% legend(Sites,'Location','SouthEast');
% xlim([1e-2, 1]); 
% ylim([1e-1 1e2]);
% set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont)
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperSize',[8 4.5],'PaperPosition',[0 0 8 4.5],'PaperPositionMode','Manual');
% print([folder_Plots,'HeightCalibration.png'],'-dpng');

%% PLOT size-conditioned Q VS u*
figure(7); clf;

Marker_bin = {'s','d','o','v'};
Color_bin = {'r','g','b','m'};

%create legend
legend_items = cell(N_d_bin,1);
for k = 1:N_d_bin
    legend_items{k} = ['d_',int2str(k),'=',num2str(d_bin_lower(k)),'-',num2str(d_bin_upper(k)),' mm'];
end

for i = 1:N_Sites
    %subplot(1,N_Sites,i); hold on;
    subplot('Position',[-0.24+0.32*i 0.13 0.27 0.82]); hold on;

    for k = 1:N_d_bin
        %plot data
        plot(tau_interval{i},Q_bin{i}(:,k),Marker_bin{k},'Color',Color_bin{k});
    end
    

    
%     %plot error bars
%     for j = 1:length(tauRe_all{i})
%         plot(ones(2,1)*tauRe_all{i}(j),Q_all{i}(j)+[-1 1]*sigma_Q_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error
%         plot(tauRe_all{i}(j)+[1 -1]*sigma_tauRe_all{i}(j),ones(2,1)*Q_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x error
%     end
    
    %format plot
    %xlim([0 0.45]);
    xlim([0 ceil(max(tau_interval{i}/0.05))*0.05]);
    %ylim([0 65]);
    xlims = xlim;
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    if i==1
        ylabel('Size-specific mass flux, $$Q_{i}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
        legend(legend_items,'Location','NorthWest');
    end
    %text(0.05*xlims(2), 62, panel_labels{i},'FontSize',9,'FontWeight','Bold');
    title(SiteNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
print([folder_Plots,'Q_d_30minute.png'],'-dpng'); %for draft
print([folder_Plots,'Q_d_30minute.tif'],'-dtiff','-r600'); %for publication