%% initialize
clearvars;
close all;

%%
%%%%%%%%%%%%%%%%
% DATA SOURCES %
%%%%%%%%%%%%%%%%

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/SizeSelective/'; %folder for saving analysis data
AnalysisData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis_Additional');
folder_Plots = '../../PlotOutput/SizeSelective/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(AnalysisData_Path); %load additional analysis data
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

%% colors for bins by wind strength
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

%%
%%%%%%%%%
% PLOTS %
%%%%%%%%%

%% plot variation in reference zqnorm with saltation flux - dref values
figure(1); clf;
for i = 1:N_Cluster
    subplot(round(N_Cluster/2),2,i); hold on;
    ind_usable = ind_usable_profile_Cluster{i};
    h10_air = plot(taunorm_profile_Cluster{i}(ind_usable),zqinorm_d10_Cluster{i}(ind_usable),'v');
    h50_air = plot(taunorm_profile_Cluster{i}(ind_usable),zqinorm_d50_Cluster{i}(ind_usable),'o');
    h90_air = plot(taunorm_profile_Cluster{i}(ind_usable),zqinorm_d90_Cluster{i}(ind_usable),'^');
    c10 = get(h10_air,'Color');
    c50 = get(h50_air,'Color');
    c90 = get(h90_air,'Color');
    xlim([1 4]);
    ylim([0 400]);
    ylims = ylim;
    text(1.1, ylims(1)+range(ylims)*0.94, Label_Cluster{i},'FontSize',12);
    if i==N_Cluster || i==N_Cluster-1
        xlabel('Normalized shear stress, $$\tau/\tau_{it}$$','Interpreter','Latex')
    end
    if mod(i,2) == 1
        ylabel('Norm. size-sel. salt. ht., $$z_{q,i}/d_{i}$$','Interpreter','Latex');
    end
    if i == 2
        legend([h10_air, h50_air, h90_air], '<d_{10,air}>','<d_{50,air}>','<d_{90,air}>','Location','West');
    end
    title(ClusterNames{i});
    set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 9],'PaperPosition',[0 0 8 9],'PaperPositionMode','Manual');
print([folder_Plots,'SaltationHeight_IndexGrainSize_NormShearStress.png'],'-dpng');


%% PLOT normalized size-conditioned zq VS tau
figure(2); clf;

for i = 1:N_Cluster
    
    %initialize subplot
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot data
    for k = 1:N_bins
        ind_plot = find(~isnan(zqinorm_bin_Cluster{i}(:,k)));
        if ~isempty(ind_plot)
            plot(taunorm_zqinorm_bin_Cluster{i}(ind_plot),zqinorm_bin_Cluster{i}(ind_plot,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
        end
    end
    
    %plot error bars
    for k = 1:N_bins
        N_tau = length(tau_profile_Cluster{i});
        for j = 1:N_tau
            plot(tau_profile_Cluster{i}(j)*[1 1],1000*(zqi_Cluster{i}(j,k)*[1 1]+sigma_zqi_Cluster{i}(j,k)*[-1 1])./d_bin_mid_Cluster(i,k),'Color',Color_bin{k});
        end
    end
    
    %plot fit
    for k = 1:N_bins
        if ~isempty(tau_zqi_fit{i}{k})
            plot(tau_zqi_fit{i}{k},1000*zqi_fit{i}{k}./d_bin_mid_Cluster(i,k),'Color',Color_bin{k}); %plot fit
        end
    end   
    
    %format plot
    xlim([0 0.5]);
    ylim([0 400]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    if i>round(N_Cluster/2)
        xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    end
    ylabel('Norm. size-sel. salt. ht., $$z_{q,i}/d_{i}$$','Interpreter','Latex');
    text(0.02, 380, Label_Cluster{i},'FontSize',9);
    
    %create legend
    legend_items = cell(N_bins,1);
    for k = 1:N_bins
        legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
    end
    h_legend = legend(legend_items,'Location','EastOutside');
    set(h_legend,'FontSize',6);

    title(ClusterNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_tau_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');




%% Plot ust-conditioned airborne versus surface size distributions
figure(3); clf;

%initialize subplots
h_subplot = gobjects(N_Cluster,1);

for i = 1:N_Cluster
        
    %initialize subplot
    if i < round(N_Cluster/2)
        h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i); hold on;
    elseif i == round(N_Cluster/2)
        h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i:(i+1)); hold on;
    else
        h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i+1); hold on;
    end
        
    %plot surface distribution
    plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
    
    %plot airborne distribution
    for k = 1:N_ust_bins
        plot(d_ustbar_airborne_mid_Cluster{i},dVdlogd_ustbar_airborne_Cluster{i}(k,:),['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_ust_bin{k})
    end
    
    %format plot
    set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlim([0.06, 3]);
    ylim([1e-4 1e1]);
    text(0.08, 6, Label_Cluster{i},'FontSize',12);

    %label plot
    title(ClusterNames{i});
    if i>round(N_Cluster/2)
        xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    end
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('Normalized volume, $$\frac{dV}{d\textrm{log}d}$$','Interpreter','Latex');
    end
        
    %create legend
    legend_items = cell(N_ust_bins+1,1);
    legend_items{1} = 'Surface';
    for j = 1:N_ust_bins
        legend_items{j+1} = ['u_{*} = ',num2str(ust_bin_min(j)),'-',num2str(ust_bin_max(j)),' m/s'];
    end
    if i == round(N_Cluster/2)
        h_legend = legend(legend_items,'Location','EastOutside');
        set(h_legend,'FontSize',12);
    end    
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[10 6],'PaperPosition',[0 0 10 6],'PaperPositionMode','Manual');
print([folder_Plots,'GSD_ust.png'],'-dpng');

% add bins to subplots
for i = 1:N_Cluster
    subplot(h_subplot(i))
    ylims = ylim;
    for j = 1:N_bins
        d_plot = [[1 1]*d_bin_lower_Cluster(i,j),[1 1]*d_bin_upper_Cluster(i,j)];
        y_plot = [ylims([2,1]),ylims];
        plot(d_plot,y_plot,'-','Color',Color_bin{j});
    end
    if i == round(N_Cluster/2)
        h_legend = legend(legend_items,'Location','EastOutside');
        set(h_legend,'FontSize',12);
    end
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[10 6],'PaperPosition',[0 0 10 6],'PaperPositionMode','Manual');
print([folder_Plots,'GSD_ust_bins_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% Plot ustnorm-conditioned airborne versus surface size distributions
figure(4); clf;

%initialize subplots
h_subplot = gobjects(N_Cluster,1);

for i = 1:N_Cluster

    %initialize subplot
    if N_Cluster == 4 %defined subplot sizes for four clusters
        if i == 1
            h_subplot(1) = subplot('position',[0.10 0.56 0.4 0.4]); hold on;
        elseif i == 2
            h_subplot(2) = subplot('position',[0.58 0.56 0.4 0.4]); hold on;
        elseif i == 3
            h_subplot(3) = subplot('position',[0.10 0.08 0.4 0.4]); hold on;
        else
            h_subplot(4) = subplot('position',[0.58 0.08 0.4 0.4]); hold on;
        end
    else %otherwise, automated subplot sizes
        if i < round(N_Cluster/2)
            h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i); hold on;
        elseif i == round(N_Cluster/2)
            h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i:(i+1)); hold on;
        else
            h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i+1); hold on;
        end
    end
        
    %plot surface distribution
    plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
    
    %plot airborne distribution
    for k = 1:N_ustnorm_bins
        plot(d_ustnormbar_airborne_mid_Cluster{i},dVdlogd_ustnormbar_airborne_Cluster{i}(k,:),['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_ustnorm_bin{k})
    end
    
    %format plot
    set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlim([0.06, 2]);
    set(gca,'xtick',[0.06:0.01:0.1, 0.2:0.1:1, 2]);
    set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','0.6','','','','1','2'});
    ylim([1e-4 1e1]);
    text(0.07, 6, Label_Cluster{i},'FontSize',12);

    %label plot
    title(ClusterNames{i});
    if i>round(N_Cluster/2)
        xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    end
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('Normalized volume, $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
    end
        
    %create legend
    legend_items = cell(N_ustnorm_bins+1,1);
    legend_items{1} = 'Surface';
    for j = 1:N_ustnorm_bins
        if j == 1
            legend_items{j+1} = ['u_{*}/u_{*,th} \leq ',num2str(ustnorm_bin_max(j),'%10.2f')];       
        else
            legend_items{j+1} = [num2str(ustnorm_bin_min(j),'%10.2f'),' < u_{*}/u_{*,th} \leq ',num2str(ustnorm_bin_max(j),'%10.2f')];
        end
    end
    if i == 1
        h_legend = legend(legend_items,'Location','SouthWest');
        set(h_legend,'FontSize',10);
    end    
%     if i == round(N_Cluster/2)
%         h_legend = legend(legend_items,'Location','EastOutside');
%         set(h_legend,'FontSize',12);
%     end    
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
print([folder_Plots,'GSD_ustnorm.png'],'-dpng');

% add bins to subplots
for i = 1:N_Cluster
    subplot(h_subplot(i))
    ylims = ylim;
    for j = 1:N_bins
        d_plot = [[1 1]*d_bin_lower_Cluster(i,j),[1 1]*d_bin_upper_Cluster(i,j)];
        y_plot = [ylims([2,1]),ylims];
        plot(d_plot,y_plot,'-','Color',Color_bin{j});
    end
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[10 6],'PaperPosition',[0 0 10 6],'PaperPositionMode','Manual');
print([folder_Plots,'GSD_ustnorm_bins_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');



%% PLOT size-conditioned normalized Q VS u* - dhat bins
figure(5); clf;

for i = 1:N_Cluster
    
    %initialize subplot
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot data
    for k = 1:N_bins
        plot(tau_profile_Cluster{i},Qhat_Cluster{i}(:,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
    end
    
    %plot error bars
    for k = 1:N_bins
        N_tau = length(tau_profile_Cluster{i});
        for j = 1:N_tau
            plot(tau_profile_Cluster{i}(j)*[1 1],Qhat_Cluster{i}(j,k)*[1 1]+sigma_Qhat_Cluster{i}(j,k)*[-1 1],'Color',Color_bin{k});
        end
    end
    
    %plot fit
    for k = 1:N_bins
        plot([tauth_Qhat_tau_fit{i}(k); tau_Qhat_fit{i}{k}],[0; Qhat_fit{i}{k}],'Color',Color_bin{k}); %plot fit
    end   
    
    %format plot
    xlim([0 0.5]);
    ylims = ylim;
    ylim([0 ylims(2)]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    if i>round(N_Cluster/2)
        xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    end
    ylabel('Bed fraction norm. size-spec. flux, $$\hat{Q_{i}}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
    text(0.05, ylims(2)*0.95, Label_Cluster{i},'FontSize',9);
    
    %create legend
    legend_items = cell(N_bins,1);
    for k = 1:N_bins
        legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
    end
    h_legend = legend(legend_items,'Location','EastOutside');
    set(h_legend,'FontSize',6);

    title(ClusterNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
print([folder_Plots,'Qhat_tau_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT Qhat,i versus tau fit values
figure(6); clf;

%plot C values and error bars
subplot(1,2,1); hold on;
for i = 1:N_Cluster
    plot(dhat_bin_mid_Cluster(i,:),C_Qhat_tau_fit{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
end
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],C_Qhat_tau_fit{i}(k)+[-1 1]*sigma_C_Qhat_tau_fit{i}(k),'Color',Color_Cluster{i}); %error bars
    end
end

%format plot
xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
ylabel('flux coefficient, $$C$$ (s)','Interpreter','Latex');
legend(ClusterNames,'Location','SouthWest');
set(gca,'Yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
xlims = xlim;
ylims = ylim;
text(xlims(1)+0.1*range(xlims), ylims(1)+10^(0.95*log10(range(ylims))), '(a)','FontSize',9);

%plot tauth values
subplot(1,2,2); hold on;
for i = 1:N_Cluster
    plot(dhat_bin_mid_Cluster(i,:),tauth_Qhat_tau_fit{i},Marker_Cluster{i},'Color',Color_Cluster{i});
end
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],tauth_Qhat_tau_fit{i}(k)+[-1 1]*sigma_tauth_Qhat_tau_fit{i}(k),'Color',Color_Cluster{i}); %error bars
    end
end

%format plot
xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
ylabel('threshold, $$\tau_{th}$$ (Pa)','Interpreter','Latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlims = xlim;
ylims = ylim;
text(xlims(1)+0.1*range(xlims), ylims(1)+0.99*range(ylims), '(b)','FontSize',9);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
print([folder_Plots,'C_tauth_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');

%plot only positive tauth values
ylim([0 ylims(2)]);
print([folder_Plots,'C_tauth_dhat_positive_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% plot f_airborne / f_surface versus d - binned
figure(7); clf;

for i = 1:N_Cluster
    
    %initialize subplot
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot data
    for k = 1:N_bins
        plot(tau_profile_Cluster{i},f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
    end
    
    %plot fit
    for k = 1:N_bins
        plot(tau_fratio_fit{i}{k},fratio_fit{i}{k},'Color',Color_bin{k}); %plot fit
    end   
    
    %format plot
    set(gca,'YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlim([0 0.5]);
    ylims = ylim;
    if i>round(N_Cluster/2)
        xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    end
    ylabel('Airborne / surface fraction, $$f_{air}/f_{sfc}$$','Interpreter','Latex');
    text(0.05, ylims(2)*0.95, Label_Cluster{i},'FontSize',9);
    
    %create legend
    legend_items = cell(N_bins,1);
    for k = 1:N_bins
        legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
    end
    h_legend = legend(legend_items,'Location','EastOutside');
    set(h_legend,'FontSize',6);

    title(ClusterNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
print([folder_Plots,'fratio_tau_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% plot f_airborne / f_surface versus d
figure(8); clf;

for i = 1:N_Cluster
    
    %initialize subplot
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot data
    for k = 1:N_bins
        plot(tau_profile_Cluster{i},f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
    end
    
    %plot fit
    for k = 1:N_bins
        plot(tau_fratio_fit{i}{k},fratio_fit{i}{k},'Color',Color_bin{k}); %plot fit
    end   
    
    %format plot
    set(gca,'YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlim([0 0.5]);
    ylims = ylim;
    xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    if i==1
        ylabel('Airborne / surface fraction, $$f_{air}/f_{bed}$$','Interpreter','Latex');
    end
    text(0.05, ylims(2)*0.95, Label_Cluster{i},'FontSize',9);
    
    %create legend
    legend_items = cell(N_bins,1);
    for k = 1:N_bins
        legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
    end
    h_legend = legend(legend_items,'Location','EastOutside');
    set(h_legend,'FontSize',6);

    title(ClusterNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
print([folder_Plots,'fratio_tau_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT mean fratio / fsurface - unbinned
figure(9); clf; hold on;

%plot fratio values
for i = 1:N_Cluster
    if strcmp(dref_type,'d50') 
        plot(d_surface_mid_Cluster{i}./d50_bar_surface_Cluster(i),dV_bar_airborne_Cluster{i}./dV_bar_surface_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    elseif strcmp(dref_type,'dmodal')
        plot(d_surface_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),dV_bar_airborne_Cluster{i}./dV_bar_surface_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    end
end

%format plot
if strcmp(dref_type,'d50')==1
    xlabel('Normalized grain size, $$d / d_{50}$$','Interpreter','Latex');
elseif strcmp(dref_type,'dmodal')==1
    xlabel('Normalized grain size, $$d / d_{modal}$$','Interpreter','Latex');
end
ylabel('Mean airborne / bed fraction, $$\langle f_{air}/f_{bed} \rangle$$','Interpreter','Latex');
legend(ClusterNames,'Location','SouthWest');
set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
print([folder_Plots,'fratiobar_tau_unbinned_',dref_type,'.png'],'-dpng');


%% PLOT mean fratio / fsurface - unbinned - limited d 
figure(10); clf; hold on;

%plot fratio values
for i = 1:N_Cluster
    ind_plot = find(d_surface_lower_Cluster{i} >= d_min & d_surface_upper_Cluster{i} <= d_max); %get only values in range of trustworthy d
    plot(d_surface_lower_Cluster{i}(ind_plot),dV_bar_airborne_Cluster{i}(ind_plot)./dV_bar_surface_Cluster{i}(ind_plot),Marker_Cluster{i},'Color',Color_Cluster{i}); %values
end

%format plot
xlabel('Grain size, $$d$$ (mm)','Interpreter','Latex');
ylabel('Mean airborne / bed fraction, $$\langle f_{air}/f_{bed} \rangle$$','Interpreter','Latex');
legend(ClusterNames,'Location','SouthWest');
set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
xlim([d_min d_max]);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
print([folder_Plots,'fratiobar_d_unbinned_limited_',dref_type,'.png'],'-dpng');



%% PLOT mean fratio / fsurface - binned
figure(11); clf; hold on;

%plot fratio values
for i = 1:N_Cluster
    plot(d_bin_mid_Cluster(i,:),f_bar_airborne_Cluster{i}./f_surface_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
end

%format plot
xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
ylabel('Mean airborne / surface fraction, $$\langle f_{air}/f_{sfc} \rangle$$','Interpreter','Latex');
legend(ClusterNames,'Location','SouthWest');
set(gca,'YScale','log','XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
print([folder_Plots,'fratiobar_tau_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT normalized size-conditioned zq VS u*
figure(12); clf;

for i = 1:N_Cluster
    
    %initialize subplot
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot data
    for k = 1:N_bins
        plot(tau_profile_Cluster{i},1000*zqi_Cluster{i}(:,k)./d_bin_mid_Cluster(i,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
    end
    
    %plot error bars
    for k = 1:N_bins
        N_tau = length(tau_profile_Cluster{i});
        for j = 1:N_tau
            plot(tau_profile_Cluster{i}(j)*[1 1],1000*(zqi_Cluster{i}(j,k)*[1 1]+sigma_zqi_Cluster{i}(j,k)*[-1 1])./d_bin_mid_Cluster(i,k),'Color',Color_bin{k});
        end
    end
    
    %plot fit
    for k = 1:N_bins
        if ~isempty(tau_zqi_fit{i}{k})
            plot(tau_zqi_fit{i}{k},1000*zqi_fit{i}{k}./d_bin_mid_Cluster(i,k),'Color',Color_bin{k}); %plot fit
        end
    end   
    
    %format plot
    xlim([0 0.5]);
    ylim([0 400]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    if i>round(N_Cluster/2)
        xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    end
    ylabel('Norm. size-sel. salt. ht., $$z_{q,i}/d_{i}$$','Interpreter','Latex');
    text(0.02, 380, Label_Cluster{i},'FontSize',9);
    
    %create legend
    legend_items = cell(N_bins,1);
    for k = 1:N_bins
        legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
    end
    h_legend = legend(legend_items,'Location','EastOutside');
    set(h_legend,'FontSize',6);

    title(ClusterNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_tau_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT mean size-conditioned zq VS dhat
figure(13); clf; hold on;

%plot data
for i = 1:N_Cluster
    plot(dhat_bin_mid_Cluster(i,:),zqinorm_bar_Cluster{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],zqinorm_bar_Cluster{i}(k)*[1 1]+zqinorm_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end
   
%format plot
xlim([0.28 2]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
set(gca,'XTick',[0.3:0.1:1,2],'XTickLabel',{'0.3','0.4','0.5','','0.7','','','1','2'});
if strcmp(dref_type,'dmodal')
    xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
else
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
end
ylabel('Mean normalized size-selective saltation height, $$\langle z_{q,i} \rangle/d_{i}$$','Interpreter','Latex');

%create legend
h_legend = legend(ClusterNames,'Location','NorthEast');
set(h_legend,'FontSize',10);

%inset plot with dimensional values
axes('Position',[.23 .2 .35 .3]); hold on;

%plot data
for i = 1:N_Cluster
    plot(d_bin_mid_Cluster(i,:),zqi_bar_Cluster{i}',[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(d_bin_mid_Cluster(i,k)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end

%format plot
%ylim([-0.1 0.15]);
xlim([0.12 1.0]);
xlabel('grain size, $$d$$ (mm)','Interpreter','Latex');
ylabel('Mean size-sel. salt. ht., $$\langle z_{q,i} \rangle$$ (m)','Interpreter','Latex');
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT mean size-conditioned zq VS d
figure(14); clf; hold on;

%plot data
for i = 1:N_Cluster
    plot(d_bin_mid_Cluster(i,:),zqi_bar_Cluster{i}',[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(d_bin_mid_Cluster(i,k)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end

%format plot
%ylim([-0.1 0.15]);
xlim([0.12 1.0]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
xlabel('Grain size, $$d$$ (mm)','Interpreter','Latex');
ylabel('Mean size-selective saltation height, $$\langle z_{q,i} \rangle$$ (m)','Interpreter','Latex');

%create legend
h_legend = legend(ClusterNames,'Location','SouthWest');
set(h_legend,'FontSize',10);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_d_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT slope of fair/fbed VS tau/tauth and slope of zq/d VS tau/tauth
figure(15); clf;

% PLOT slope of fair/fbed
subplot(2,1,1); hold on;

%plot data
for i = 1:N_Cluster
%    plot(dhat_bin_mid_Cluster(i,:),b_fratio_taunorm_fit{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
    plot(dhat_bin_mid_Cluster(i,:),b_fratio_taunorm_fit{i}./fratio_bar_Cluster{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
%        plot(dhat_bin_mid_Cluster(i,k)*[1 1],b_fratio_taunorm_fit{i}(k)*[1 1]+sigma_b_fratio_taunorm_fit{i}(k)*[-1 1],'Color',Color_Cluster{i});
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],(b_fratio_taunorm_fit{i}(k)*[1 1]+sigma_b_fratio_taunorm_fit{i}(k)*[-1 1])./fratio_bar_Cluster{i}(k),'Color',Color_Cluster{i});
    end
end
   
%format plot
xlim([0.4 2]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
set(gca,'XTick',[0.4:0.1:1,2],'XTickLabel',{'0.4','0.5','','0.7','','','1','2'});
if strcmp(dref_type,'dmodal')
    xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
else
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
end
%ylabel('Airborne fraction vs. shear stress trend, $$\frac{\textrm{d}f_{air}/f_{bed}}{\textrm{d}\tau/\tau_{th}}$$','Interpreter','Latex');
ylabel('Norm air/bed frac vs stress trend, $$(\frac{\textrm{d}f_{air}/f_{bed}}{\textrm{d}\tau/\tau_{th}})/\langle f_{air}/f_{bed} \rangle$$','Interpreter','Latex');

% PLOT slope of zq/d
subplot(2,1,2); hold on;

%plot data
for i = 1:N_Cluster
%    plot(dhat_bin_mid_Cluster(i,:),b_zqinorm_taunorm_fit{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
    plot(dhat_bin_mid_Cluster(i,:),b_zqinorm_taunorm_fit{i}./zqinorm_bar_Cluster{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
%        plot(dhat_bin_mid_Cluster(i,k)*[1 1],b_zqinorm_taunorm_fit{i}(k)*[1 1]+sigma_b_zqinorm_taunorm_fit{i}(k)*[-1 1],'Color',Color_Cluster{i});
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],(b_zqinorm_taunorm_fit{i}(k)*[1 1]+sigma_b_zqinorm_taunorm_fit{i}(k)*[-1 1])./zqinorm_bar_Cluster{i}(k),'Color',Color_Cluster{i});
    end
end
   
%format plot
xlim([0.4 2]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
set(gca,'XTick',[0.4:0.1:1,2],'XTickLabel',{'0.4','0.5','','0.7','','','1','2'});
if strcmp(dref_type,'dmodal')
    xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
else
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
end
%ylabel('Saltation height vs. shear stress trend, $$\frac{\textrm{d}z_{q,i}/d_{i}}{\textrm{d}\tau/\tau_{th}}$$','Interpreter','Latex');
ylabel('Norm salt ht vs stress trend, $$(\frac{\textrm{d}z_{q,i}/d_{i}}{\textrm{d}\tau/\tau_{th}})/\langle z_{q,i} \rangle$$','Interpreter','Latex');

%create legend
h_legend = legend(ClusterNames,'Location','NorthWest');
set(h_legend,'FontSize',10);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 8],'PaperPosition',[0 0 6 8],'PaperPositionMode','Manual');
print([folder_Plots,'fratio_zqinorm_taunorm_d_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');





%% plot variation in z of BSNEs
figure(16); clf;
for i = 1:N_Cluster
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot z
    for j = 1:length(z_profile_Cluster{i})
        plot(j*ones(size(z_profile_Cluster{i}{j})),z_profile_Cluster{i}{j},'ko');
    end
    
    if i>round(N_Cluster/2)
        xlabel('number of profile','Interpreter','Latex')
    end
    
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('BSNE $$z$$ (m)','Interpreter','Latex');
    end

    title(ClusterNames{i});
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[10 7],'PaperPosition',[0 0 10 7],'PaperPositionMode','Manual');
print([folder_Plots,'z_profile_BSNE.png'],'-dpng');


%% plot variation in z of BSNEs
figure(17); clf;
for i = 1:N_Cluster
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot z - usable profiles
    ind_usable = ind_usable_profile_Cluster{i};
    z_plot = z_profile_Cluster{i}(ind_usable);
    Time_plot = Time_profile_Cluster{i}(ind_usable);
    N_plot = length(z_plot);
    for j = 1:N_plot
        Time_profile = []; for k = 1:length(z_plot{j}); Time_profile = [Time_profile; Time_plot(j)]; end; %get times for profile
        h1 = plot(Time_profile,z_plot{j},'bo-');
    end
    
    %plot z - unusable profiles
    ind_unusable = setdiff(1:length(z_profile_Cluster{i}),ind_usable);
    z_plot = z_profile_Cluster{i}(ind_unusable);
    Time_plot = Time_profile_Cluster{i}(ind_unusable);
    N_plot = length(z_plot);
    for j = 1:N_plot
        Time_profile = []; for k = 1:length(z_plot{j}); Time_profile = [Time_profile; Time_plot(j)]; end; %get times for profile
        h2 = plot(Time_profile,z_plot{j},'rx-');
    end
        
    if i>round(N_Cluster/2)
        xlabel('number of profile','Interpreter','Latex')
    end
    
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('BSNE $$z$$ (m)','Interpreter','Latex');
    end

    if i == 1
        legend([h1 h2],{'usable','unusable'},'Location','North');
    end
    
    title(ClusterNames{i});
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[10 7],'PaperPosition',[0 0 10 7],'PaperPositionMode','Manual');
print([folder_Plots,'z_profile_BSNE.png'],'-dpng');