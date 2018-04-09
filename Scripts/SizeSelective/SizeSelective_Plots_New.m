%% initialize
clearvars;
close all;

%%
%%%%%%%%%%%%%%%%
% DATA SOURCES %
%%%%%%%%%%%%%%%%

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/SizeSelective/'; %folder for loading/saving analysis data
folder_LitData = '../../AnalysisData/Literature/'; %folder for loading/saving literature data
AnalysisData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis');
LitData_Path = strcat(folder_LitData,'LitAnalysis');
folder_Plots = '../../PlotOutput/SizeSelective/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(AnalysisData_Path); %load analysis data
load(LitData_Path); %load literature data
addpath(folder_Functions); %point MATLAB to location of functions

%%
%%%%%%%%%%%%%%%%%
% PLOTTING INFO %
%%%%%%%%%%%%%%%%%

%% plotting information
PlotFont = 12; %font for labels
LineWidth_Plot = 1; %width of lines
Marker_Cluster = {'s','d','o','p','h','^','v','>'}; %markers for Clusters
Color_Cluster = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]}; %colors for sites / clusters
Label_Cluster = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'}; %markers for Clusters
Marker_bin = {'o','s','d','^','v','>','<','p','h','x','*','+'}; %markers for bins

%% colors for bins by wind strength or other property
Color_bin = cell(N_bins,1);
for i = 1:N_bins
    Color_bin{i} = [0, (i-1)/(N_bins-1), 1-(i-1)/(N_bins-1)];
end
Color_ust_bin = cell(N_ust_bins,1);
for i = 1:N_ust_bins
    Color_ust_bin{i} = [(i-1)/(N_ust_bins-1), 0, 1-(i-1)/(N_ust_bins-1)];
end
Color_ustnorm_bin = cell(N_ustnorm_bins,1);
for i = 1:N_ustnorm_bins
    Color_ustnorm_bin{i} = [(i-1)/(N_ustnorm_bins-1), 0, 1-(i-1)/(N_ustnorm_bins-1)];
end
Color_taunorm_bin = cell(N_taunorm_bins,1);
for i = 1:N_taunorm_bins
    Color_taunorm_bin{i} = [(i-1)/(N_taunorm_bins-1), 0, 1-(i-1)/(N_taunorm_bins-1)];
end
Color_z_Namikas06 = cell(N_z_Namikas06,1);
for i = 1:N_z_Namikas06
    Color_z_Namikas06{i} = [(i-1)/(N_z_Namikas06-1), 0, 1-(i-1)/(N_z_Namikas06-1)];
end

%%
%%%%%%%%%
% PLOTS %
%%%%%%%%%

% %% plot variation in reference grain sizes with saltation flux
% figure(1); clf;
% for i = 1:N_Cluster
%     subplot(round(N_Cluster/2),2,i); hold on;
%     ind_usable = ind_usable_profile_Cluster{i};
%     h90_air = plot(taunorm_profile_Cluster{i}(ind_usable),d90_profilebar_airborne_Cluster{i}(ind_usable),'^');
%     h50_air = plot(taunorm_profile_Cluster{i}(ind_usable),d50_profilebar_airborne_Cluster{i}(ind_usable),'o');
%     h10_air = plot(taunorm_profile_Cluster{i}(ind_usable),d10_profilebar_airborne_Cluster{i}(ind_usable),'v');
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     h90_sfc = plot([1 4],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_sfc = plot([1 4],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_sfc = plot([1 4],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     xlim([1 4]);
%     ylim([0 0.9]);
%     text(1.05, 0.8, Label_Cluster{i},'FontSize',12);
%     if i==N_Cluster || i==N_Cluster-1
%         xlabel('shear stress ratio, $$\tau/\tau_{it}$$','Interpreter','Latex')
%     end
%     if mod(i,2) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air], 'd_{90,air}','d_{50,air}','d_{10,air}');
%     elseif i == N_Cluster
%         legend([h90_sfc, h50_sfc, h10_sfc], 'd_{90,bed}','d_{50,bed}','d_{10,bed}');
%     end
%     title(ClusterNames{i});
% 
%     set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[7 9],'PaperPosition',[0 0 7 9],'PaperPositionMode','Manual');
% print([folder_Plots,'IndexGrainSize_NormShearStress.png'],'-dpng');
% 
% 
% %% plot variation in reference grain sizes with shear velocity
% figure(2); clf;
% for i = 1:N_Cluster
%     subplot(2,round(N_Cluster/2),i); hold on;
%     h90_air = plot(ust_profile_Cluster{i},d90_profilebar_airborne_Cluster{i},'^');
%     h50_air = plot(ust_profile_Cluster{i},d50_profilebar_airborne_Cluster{i},'o');
%     h10_air = plot(ust_profile_Cluster{i},d10_profilebar_airborne_Cluster{i},'v');
%     hbar_air = plot(ust_profile_Cluster{i},dbar_profilebar_airborne_Cluster{i},'s');
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     cbar = get(hbar_air,'Color');
%     h90_bed = plot([0.2 0.6],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_bed = plot([0.2 0.6],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_bed = plot([0.2 0.6],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     hbar_bed = plot([0.2 0.6],dbar_bar_surface_Cluster(i)*[1 1],'Color',cbar);
%     xlim([0.2 0.6]);
%     ylim([0 0.9]);
%     if i>round(N_Cluster/2)
%         xlabel('shear velocity, $$u_{*}$$ (m s$$^{-1}$$)','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air,hbar_air], 'd_{90,air}','d_{50,air}','d_{10,air}','d_{bar,air}');
%     elseif i == N_Cluster
%         legend([h90_bed, h50_bed, h10_bed,hbar_bed], 'd_{90,bed}','d_{50,bed}','d_{10,bed}','d_{bar,bed}');
%     end
%     title(ClusterNames{i});
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 7],'PaperPosition',[0 0 6 7],'PaperPositionMode','Manual');
% print([folder_Plots,'ReferenceGrainSize_ShearVelocity.png'],'-dpng');
% 
% 
% %% plot variation in reference grain sizes with shear stress
% figure(3); clf;
% for i = 1:N_Cluster
%     subplot(2,round(N_Cluster/2),i); hold on;
%     
%     %plot airborne sizes
%     h90_air = plot(tau_profile_Cluster{i},d90_profilebar_airborne_Cluster{i},'^');
%     h50_air = plot(tau_profile_Cluster{i},d50_profilebar_airborne_Cluster{i},'o');
%     h10_air = plot(tau_profile_Cluster{i},d10_profilebar_airborne_Cluster{i},'v');
%     hbar_air = plot(tau_profile_Cluster{i},dbar_profilebar_airborne_Cluster{i},'s');
%     
%     %get info about airborne plots
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     cbar = get(hbar_air,'Color');
%     
%     %plot airborne uncertainties
%     for j = 1:length(tau_profile_Cluster{i})
%         plot(tau_profile_Cluster{i}(j)*[1 1],d90_profilebar_airborne_Cluster{i}(j)*[1 1]+d90_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',c90);
%         plot(tau_profile_Cluster{i}(j)*[1 1],d50_profilebar_airborne_Cluster{i}(j)*[1 1]+d50_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',c50);
%         plot(tau_profile_Cluster{i}(j)*[1 1],d10_profilebar_airborne_Cluster{i}(j)*[1 1]+d10_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',c10);
%         plot(tau_profile_Cluster{i}(j)*[1 1],dbar_profilebar_airborne_Cluster{i}(j)*[1 1]+dbar_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',cbar);
%     end
%         
%     %make surface plots
%     h90_bed = plot([0 0.45],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_bed = plot([0 0.45],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_bed = plot([0 0.45],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     hbar_bed = plot([0 0.45],dbar_bar_surface_Cluster(i)*[1 1],'Color',cbar);       
%     xlim([0 0.45]);
%     ylim([0 0.9]);
%     if i>round(N_Cluster/2)
%         xlabel('shear stress, $$\tau$$ (Pa)','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air,hbar_air], 'd_{90,air}','d_{50,air}','d_{10,air}','d_{bar,air}');
%     elseif i == N_Cluster
%         legend([h90_bed, h50_bed, h10_bed,hbar_bed], 'd_{90,bed}','d_{50,bed}','d_{10,bed}','d_{bar,bed}');
%     end
%     title(ClusterNames{i});
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[10 7],'PaperPosition',[0 0 10 7],'PaperPositionMode','Manual');
% print([folder_Plots,'ReferenceGrainSize_ShearStress.png'],'-dpng');
% 
% %% plot variation in reference grain sizes with height
% figure(4); clf;
% for i = 1:N_Cluster
%     subplot(2,round(N_Cluster/2),i); hold on;
%     h90_air = plot(z_airborne_Cluster{i},d90_airborne_Cluster{i},'^');
%     h50_air = plot(z_airborne_Cluster{i},d50_airborne_Cluster{i},'o');
%     h10_air = plot(z_airborne_Cluster{i},d10_airborne_Cluster{i},'v');
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     h90_sfc = plot([0 0.6],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_sfc = plot([0 0.6],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_sfc = plot([0 0.6],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     xlim([0 0.6]);
%     ylim([0 0.9]);
%     if i>round(N_Cluster/2)
%         xlabel('height above surface, $$z$$ (m)','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air], 'd_{90,air}','d_{50,air}','d_{10,air}');
%     elseif i == N_Cluster
%         legend([h90_sfc, h50_sfc, h10_sfc], 'd_{90,surface}','d_{50,surface}','d_{10,surface}');
%     end
%     title(ClusterNames{i});
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'ReferenceGrainSize_Height.png'],'-dpng');
% 
% 
% %% plot variation in reference grain sizes with normalized height
% figure(5); clf;
% for i = 1:N_Cluster
%     subplot(2,round(N_Cluster/2),i); hold on;
%     h90_air = plot(znorm_airborne_Cluster{i},d90_airborne_Cluster{i},'^');
%     h50_air = plot(znorm_airborne_Cluster{i},d50_airborne_Cluster{i},'o');
%     h10_air = plot(znorm_airborne_Cluster{i},d10_airborne_Cluster{i},'v');
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     h90_sfc = plot([0 7],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_sfc = plot([0 7],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_sfc = plot([0 7],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     xlim([0 7]);
%     ylim([0 0.9]);
%     if i>round(N_Cluster/2)
%         xlabel('norm. ht above sfc., $$z/z_{q}$$','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air], 'd_{90,air}','d_{50,air}','d_{10,air}');
%     elseif i == N_Cluster
%         legend([h90_sfc, h50_sfc, h10_sfc], 'd_{90,surface}','d_{50,surface}','d_{10,surface}');
%     end
%     title(ClusterNames{i});
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'ReferenceGrainSize_NormHeight.png'],'-dpng');

% %% plot variation in airborne size distributions with height
% % note - don't include top bin, as its midpoint is a bit uncertain
% figure(6); clf;
% 
% %initialize subplots
% h_subplot = gobjects(N_Cluster,1);
% 
% for i = 1:N_Namikas06
%     
%     %initialize subplot
%     if i == 1
%         h_subplot(1) = subplot('position',[0.1 0.1 0.25 0.85]); hold on;
%     elseif i == 2
%         h_subplot(2) = subplot('position',[0.41 0.1 0.25 0.85]); hold on;
%     else
%         h_subplot(3) = subplot('position',[0.72 0.1 0.25 0.85]); hold on;
%     end
%         
%     %plot airborne distribution
%     for k = 1:N_z_Namikas06
%         plot(d_Namikas06(1:end-1),...
%             dVdlogd_Namikas06{i}(k,1:end-1),...
%             ['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_z_Namikas06{k})
%     end
%     
%     %format plot
%     set(gca,'XScale','log','YScale','log','YMinorTick','On','Box','On');
%     xlim([0.06, 2]);
%     set(gca,'xtick',[0.06:0.01:0.1,0.2:0.1:2]);
%     set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','','0.7','','','1','','','','','','','','','','2'});
%     ylim([1e-3 1e1]);
% 
%     %label plot
%     htitle = title(['u_{*} = ',num2str(ust_Namikas06(i)),' m/s']);
%     set(htitle,'Position',[0.35,3.5]); %set title below edge of box to accommodate second axis
%     xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
%     if mod(i,3) == 1
%         ylabel('Normalized volume size distr., $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
%     end
%     
%     %create legend
%     legend_items = cell(N_z_Namikas06,1);
%     for j = 1:N_z_Namikas06
%         %legend_items{j} = [num2str(z_bottom_Namikas06(j),'%10.2f'),' < z \leq ',num2str(z_top_Namikas06(j),'%10.2f')];
%         legend_items{j} = [num2str(znorm_bottom_Namikas06(j),'%10.2f'),' < z/z_{q} \leq ',num2str(znorm_top_Namikas06(j),'%10.2f')];
%     end
%     if i == 1
%         h_legend = legend(legend_items,'Location','SouthWest');
%         set(h_legend,'FontSize',8);
%     end
% end
% 
% %print plot - landscape
% set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
% print([folder_Plots,'GSD_Namikas06_z.png'],'-dpng');


%% plot variation in mean airborne grain sizes (normalized by mean surface)
% versus height (normalized by z_q)
figure(7); clf;

%initialize subplots
h_subplot = gobjects(8,1);

% Namikas
h_subplot(1) = subplot('position',[0.1 0.7 0.25 0.24]); hold on;

%plot airborne profiles
for i = 1:N_Namikas06
    plot(znorm_Namikas06,dbar_airborne_Namikas06{i}./...
        mean(dbar_bar_surface_Cluster(3:6)),[Marker_bin{i},'-']);
end
%for normalization, use mean surface grain size from all of the Oceano
%clusters

%format plot
htitle = title('Namikas 2006');
set(gca,'YMinorTick','On','Box','On');
ylabel('Norm. mean grain diam., $$\bar{d}_{air}/\bar{d}_{bed}$$ (mm)','Interpreter','Latex')
xlim([0 6]);
ylim([0.5 1]);

%create legend
legend_items = cell(N_Namikas06,1);
for i = 1:N_Namikas06
    legend_items{i} = ['u_{*} = ',num2str(ust_Namikas06(i),'%10.2f'),' m/s'];
end
h_legend = legend(legend_items,'Location','NorthEast');
set(h_legend,'FontSize',8);

%% Farrell
h_subplot(2) = subplot('position',[0.41 0.7 0.25 0.24]); hold on;

%plot airborne profiles
for i = 1:N_ustbins_Farrell12
    plot(znorm_airborne_ustbin_Farrell12{i},dbar_airborne_ustbin_Farrell12{i}./...
        mean(dbar_bar_surface_Cluster(1)),[Marker_bin{i},'-']);
end
%for normalization, use mean surface grain size from Jericoacoara

%format plot
htitle = title('Farrell 2012');
set(gca,'YMinorTick','On','Box','On');
xlim([0 6]);
ylim([0.5 1]);

%create legend
for i = 1:N_ustbins_Farrell12
    legend_items{i} = [num2str(ust_lower_Farrell12(i),'%10.2f'),'\leq u_{*} \leq',num2str(ust_upper_Farrell12(i),'%10.2f'),' m/s'];
end
h_legend = legend(legend_items,'Location','NorthEast');
set(h_legend,'FontSize',8);


%% our data
for i = 1:N_Cluster
    if i == 1
        h_subplot(3) = subplot('position',[0.1 0.38 0.25 0.24]); hold on;
    elseif i == 2
        h_subplot(4) = subplot('position',[0.41 0.38 0.25 0.24]); hold on;
    elseif i == 3 %make this subplot bigger to accommodate legend
        %h_subplot(5) = subplot('position',[0.72 0.38 0.25 0.24]); hold on;
        h_subplot(5) = subplot('position',[0.72 0.38 0.25 0.385]); hold on;
    elseif i == 4
        h_subplot(6) = subplot('position',[0.1 0.06 0.25 0.24]); hold on;
    elseif i == 5
        h_subplot(7) = subplot('position',[0.41 0.06 0.25 0.24]); hold on;
    else
        h_subplot(8) = subplot('position',[0.72 0.06 0.25 0.24]); hold on;
    end

    %plot airborne profiles
    for j = 1:N_taunorm_bins
        ind_plot = find(~isnan(dbar_profile_airborne_taunorm_Cluster{i}(j,:)));
        if ~isempty(ind_plot) %plot only profiles with data
            plot(znorm_profile_airborne_taunorm_Cluster{i}(j,ind_plot),...
                dbar_profile_airborne_taunorm_Cluster{i}(j,ind_plot)./...
                dbar_bar_surface_Cluster(i),[Marker_bin{j},'-']);
        else %for those profiles with no data, create a dummy plot so that legend renders properly
            plot([0 0],[0 0],[Marker_bin{j},'-']); 
        end
    end

    %format plot
    htitle = title(ClusterNames{i});
    set(gca,'YMinorTick','On','Box','On');
    if mod(i,3) == 1
        ylabel('Norm. mean grain diam., $$\bar{d}_{air}/\bar{d}_{bed}$$ (mm)','Interpreter','Latex')
    end
    if i>=2*round(N_Cluster/2)-2
        xlabel('Normalized height, $$z/z_q$$','Interpreter','Latex')
    end
    xlim([0 6]);
    ylim([0.5 1]);

    %create legend
    if i == 3
        for j = 1:N_taunorm_bins
            legend_items{j} = [num2str(taunorm_min_bins(j),'%10.1f'),'\leq \tau/\tau_{it} \leq',num2str(taunorm_max_bins(j),'%10.1f')];
        end
        h_legend = legend(legend_items,'Location','NorthOutside');
        set(h_legend,'FontSize',8);
    end
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[9 9],'PaperPosition',[0 0 9 9],'PaperPositionMode','Manual');
print([folder_Plots,'dbar_profile_comparison.png'],'-dpng');