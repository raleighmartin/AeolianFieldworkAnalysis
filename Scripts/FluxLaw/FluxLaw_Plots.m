%% PLOTS FOR ANALYSIS OF SALTATION FLUX VERSUS SHEAR STRESS RELATIONSHIP

%%
%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%

%initialize
clearvars;
close all;

%parameters
g = 9.8; %gravity (m/s^2)
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)

%% folders for loading data and functions
folder_LoadData = '../../AnalysisData/FluxLaw/'; %folder for retrieving data for this analysis
folder_Functions = '../Functions/'; %folder with functions

%% paths for loading and printing plots - restricted
LoadData_Path = strcat(folder_LoadData,'FluxLawCalcs_30min_Restricted'); %path for 30 minute data
folder_Plots = '../../PlotOutput/FluxLaw/'; %folder for plots

% %% paths for loading and printing plots - restricted - alternate data
% LoadData_Path = strcat(folder_LoadData,'FluxLawCalcs_30min_Restricted_alt'); %path for 30 minute data
% folder_Plots = '../../PlotOutput/FluxLaw/AltAnemometer/'; %folder for plots

%load data
load(LoadData_Path);
addpath(folder_Functions); %point MATLAB to location of functions

%set info for plotting
Markers_Field = {'s','d','o'};
Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
Markers_Lit = {'<','>','^'};
Colors_Lit = {[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
%Markers_WindTunnel = {'+','x','*'};
Markers_WindTunnel = {'<','>','^'};
Colors_WindTunnel = {[0.48, 0.429, 0.372], [0.3835, 0.7095, 0.5605],[0.3975, 0.4656, 0.7445]};
panel_labels = {'A','B','C'};

%set info about sizing
PlotFont = 12;
MarkerSize_plot = 6;
LineWidth_plot = 0.5;

%% Generate values for plotting Q versus tau fits

%initialize values for plotting fits
tau_linearfit_plot = cell(N_Sites,1); %values for plotting linear fits
Q_linearfit_plot = cell(N_Sites,1); %values for plotting linear fits
ust_threehalvesfit_plot = cell(N_Sites,1); %values for plotting 3/2 fits
tau_threehalvesfit_plot = cell(N_Sites,1); %values for plotting 3/2 fits
Q_threehalvesfit_plot = cell(N_Sites,1); %values for plotting 3/2 fits

% generate values for plotting fits
for i = 1:N_Sites
    %get values for plotting fits
    tau_linearfit_plot{i} = linspace(tauit_linearfit_all(i),max(tau_fit_all{i}),50);
    Q_linearfit_plot{i} = C_linearfit_all(i)*(tau_linearfit_plot{i}-tauit_linearfit_all(i));
    tau_threehalvesfit_plot{i} = linspace(tauit_threehalvesfit_all(i),max(tau_fit_all{i}),50);
    ust_threehalvesfit_plot{i} = sqrt(tau_threehalvesfit_plot{i}/rho_a(i));
    Q_threehalvesfit_plot{i} = C_threehalvesfit_all(i)*...
        ust_threehalvesfit_plot{i}.*...
        (tau_threehalvesfit_plot{i}-tauit_threehalvesfit_all(i));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMARY PLOTS FOR PAPER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT ZQ VERSUS UST
figure(1); clf; 
legend_items = [SiteNames;LitNames];

%PANEL A - dimensional
%subplot(2,7,1:7); hold on;
subplot('Position',[0.1 0.6 0.88 0.38]); hold on;

%plot binned Field data
for i = 1:N_Sites
    plot(ust_zqustfit_all{i},zq_zqustfit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot);
end
%plot Literature data
for i = 1:N_Lit
    plot(ust_Lit{i},zq_Lit{i},Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_plot);
end

%plot error bars separately
for i = 1:N_Sites %plot binned Field data
    for j = 1:length(ust_zqustfit_all{i})
        plot(ones(2,1)*ust_zqustfit_all{i}(j), zq_zqustfit_all{i}(j)+[-1 1]*sigma_zq_zqustfit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error-bars
    end
end
for i = 1:N_Sites %plot Literature data
    for j = 1:length(ust_Lit{i})
        plot(ones(2,1)*ust_Lit{i}(j), zq_Lit{i}(j)+[-1 1]*sigma_zq_Lit{i}(j),'Color',Colors_Lit{i},'LineWidth',LineWidth_plot); %y error-bars
    end
end

%organize plot
xlim([0.25 0.65]);
ylim([0 0.15]);
legend(legend_items,'Location','EastOutside');
text(0.67,0.118,'Field measurements:','FontSize',PlotFont)
text(0.255, 0.14,'A','FontSize',9,'FontWeight','Bold');
set(gca,'XMinorTick','On','XScale','log','YMinorTick','On','Box','On');
set(gca,'XMinorTick','On');
xlabel('Wind shear velocity, $$u_{*}$$ (ms$$^{-1}$$)','Interpreter','Latex');
ylabel('Saltation layer height, $$z_q$$ (m)','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%PANEL B - fit values
%subplot(2,7,8:9); hold on;
subplot('Position',[0.1 0.08 0.33 0.42]); hold on;

%plot Field data
for i = 1:N_Sites
    plot(d50_Site(i),slope_zqustfit_all(i),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot); %values
    plot(ones(2,1)*d50_Site(i),slope_zqustfit_all(i)+[-1 1]*sigma_slope_zqustfit_all(i),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %error bars
end
%plot Literature data
for i = 1:2 %neglect Farrell (2012) data
    plot(d50_Lit(i),slope_zqustfit_Lit_all(i),Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_plot); %values
    plot(ones(2,1)*d50_Lit(i),slope_zqustfit_Lit_all(i)+[-1 1]*sigma_slope_zqustfit_Lit_all(i),'Color',Colors_Lit{i},'LineWidth',LineWidth_plot); %error bars
end

%plot 0 line
plot([0.2 0.6],[0 0],'k--','LineWidth',LineWidth_plot);

%organize plot
xlim([0.2 0.6]);
ylim([-0.2 0.2]);
text(0.21, 0.175,'B','FontSize',9,'FontWeight','Bold');
set(gca,'XMinorTick','On','XScale','log','YMinorTick','On','Box','On');
xlabel('Particle diameter, $$d_{50}$$ (mm)','Interpreter','Latex');
ylabel('$$z_q$$-$$u_{*}$$ fit slope, $$b$$ (m/ms$$^{-1}$$)','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%PANEL B2 - fit values for Farrell et al. (2012), without d50
%subplot(2,7,10); hold on;
subplot('Position',[0.44 0.08 0.05 0.42]); hold on;

%plot best fit values
plot([0 0],slope_zqustfit_Lit_all(3),Markers_Lit{3},'Color',Colors_Lit{3},'MarkerSize',MarkerSize_plot); %values
plot([0 0],slope_zqustfit_Lit_all(3)+[-1 1]*sigma_slope_zqustfit_Lit_all(3),'Color',Colors_Lit{3},'LineWidth',LineWidth_plot); %error bars

%plot 0 line
plot([-1 1],[0 0],'k--','LineWidth',LineWidth_plot);

%organize plot
xlim([-1 1]);
ylim([-0.2 0.2]);
set(gca,'XTick',[],'XTickLabel',{''},'YTickLabel',{''},'YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

%PANEL C - dimensionless saltation heights
%subplot(2,7,12:14); hold on;
subplot('Position',[0.61 0.08 0.38 0.42]); hold on;

%plot dimensionless heights versus d50 for wind tunnel data
for i = 1:length(WindTunnelNames)
    plot(d50_WindTunnel{i}, zqnorm_WindTunnel{i}, Markers_WindTunnel{i},'Color',Colors_WindTunnel{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot,'MarkerFaceColor',Colors_WindTunnel{i}); %values
end
%plot dimensionless heights versus d50 for Field data
for i = 1:N_Sites
    plot(d50_Site(i), zqnorm_bar_all(i), Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot); %values
    plot(ones(2,1)*d50_Site(i), zqnorm_bar_all(i)+[-1 1]*sigma_zqnorm_bar_all(i),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %error bars
end
%plot dimensionless heights versus d50 for lit data
for i = 1:2 %only plot first two entries, ignorning Farrell (2012) with no d50
    plot(d50_Lit(i), zqnorm_bar_Lit_all(i), Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_plot); %values
    plot(ones(2,1)*d50_Lit(i), zqnorm_bar_Lit_all(i)+[-1 1]*sigma_zqnorm_bar_Lit_all(i),'Color',Colors_Lit{i},'LineWidth',LineWidth_plot); %error bars
end
%organize plot
xlim([0.2 0.7]);
ylim([0 250]);
legend(WindTunnelNames,'Location','NorthOutside'); %add legend
text(0.235,380,'Wind tunnel measurements:','FontSize',PlotFont)
text(0.205, 230,'C','FontSize',9,'FontWeight','Bold');
set(gca,'XMinorTick','On','XScale','log','YMinorTick','On','Box','On');
xlabel('Particle diameter, $$d_{50}$$ (mm)','Interpreter','Latex');
ylabel('Dimensionless height, $$\langle z_q \rangle/d_{50}$$','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'zq_ust.png'],'-dpng'); %for draft
print([folder_Plots,'zq_ust.tif'],'-dtiff','-r600'); %for publication


%% PLOT Q VERSUS TAU AS SINGLE PLOT
figure(2); clf; hold on;

%plotting
for i = 1:N_Sites
    plot(tau_fit_all{i},Q_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot); %plot binned values
    plot(tau_linearfit_plot{i},Q_linearfit_plot{i},'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %plot linear fit
    plot(tau_threehalvesfit_plot{i},Q_threehalvesfit_plot{i},'--','Color',Colors_Field{i},'LineWidth',LineWidth_plot); %plot 3/2 fit
end

%generate of legend items
legend_items = cell(N_Sites*3,1);
for i = 1:N_Sites
    legend_items{3*i-2} = SiteNames{i}; %add to list of legend items
    legend_items{3*i-1} = 'linear fit';
    legend_items{3*i} = 'nonlinear 3/2 fit';
end

%plot error bars
for i = 1:N_Sites
    for j = 1:length(tau_fit_all{i})
        plot(ones(2,1)*tau_fit_all{i}(j),Q_fit_all{i}(j)+[-1 1]*sigma_Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error
        plot(tau_fit_all{i}(j)+[1 -1]*sigma_tau_fit_all{i}(j),ones(2,1)*Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x error
    end
end

%format plot
xlim([0 0.45]);
ylim([0 65]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('Wind shear stress, $$\tau_{i}$$ (Pa)','Interpreter','Latex');
ylabel('Saltation mass flux, $$Q_{i}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf,'PaperUnits','inches','PaperSize',[7.3 5],'PaperPosition',[0 0 7.3 5],'PaperPositionMode','Manual');
print([folder_Plots,'Q_tau.png'],'-dpng'); %for draft
print([folder_Plots,'Q_tau.tif'],'-dtiff','-r600'); %for publication


%% PLOT Q VERSUS TAU AS SUBPLOTS
figure(2); clf;

for i = 1:N_Sites
    %subplot(N_Sites,1,i); hold on;
    subplot('Position',[0.15 1.03-0.32*i 0.8 0.25]); hold on;
    
    %plot binned values
    plot(tau_fit_all{i},Q_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);
    
    %plot fit lines
    plot(tau_linearfit_plot{i},Q_linearfit_plot{i},'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %plot linear fit
    plot(tau_threehalvesfit_plot{i},Q_threehalvesfit_plot{i},'--','Color',Colors_Field{i},'LineWidth',LineWidth_plot); %plot 3/2 fit

    %plot error bars
    for j = 1:length(tau_fit_all{i})
        plot(ones(2,1)*tau_fit_all{i}(j),Q_fit_all{i}(j)+[-1 1]*sigma_Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error
        plot(tau_fit_all{i}(j)+[1 -1]*sigma_tau_fit_all{i}(j),ones(2,1)*Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x error
    end

    %format plot
    xlim([0 0.45]);
    ylim([0 ceil(max(Q_threehalvesfit_plot{i}/5))*5]);
%     xlim([0 ceil(max(tau_fit_all{i}/0.05))*0.05]);
%     ylim([0 65]);
    xlims = xlim; ylims = ylim; %get x and y limit information
    text(0.02*xlims(2),0.05*ylims(2),panel_labels{i},'FontSize',9.5,'FontWeight','Bold');
    set(gca,'YTick',0:10:ylims(2));
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    set(gca, 'LooseInset', get(gca,'TightInset'))
    if i==3
        xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    end
    ylabel('Saltation flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
    title(SiteNames{i});
%     legend('measurements','linear fit','nonlinear 3/2 fit','Location','NorthWest');
    legend('meas.','linear fit','3/2 fit','Location','NorthWest');
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[3.5 8],'PaperPosition',[0 0 3.5 8],'PaperPositionMode','Manual');
print([folder_Plots,'Q_tau_panels.png'],'-dpng'); %for draft
print([folder_Plots,'Q_tau_panels.tif'],'-dtiff','-r600'); %for publication


%% PLOT NORMALIZED Q [Q/(ustit/g)] VERSUS TAUEX BINNED DATA
figure(3); clf; hold on;

%plot binned values
for i = 1:N_Sites
    plot(tauex_fit_all{i},Qnorm_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot);
end

%plot error bars
for i = 1:N_Sites
    for j = 1:length(tauex_fit_all{i})
        plot(ones(2,1)*tauex_fit_all{i}(j),Qnorm_fit_all{i}(j)+[-1 1]*sigma_Qnorm_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error
        plot(tauex_fit_all{i}(j)+[1 -1]*sigma_tauex_fit_all{i}(j),ones(2,1)*Qnorm_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x error
    end
end

%plot fit values
tauex_fit = [0 0.35];
for i = 1:N_Sites %for individual sites
    plot(tauex_fit, CQ_all(i)*tauex_fit,'--','Color',Colors_Field{i},'LineWidth',LineWidth_plot);
end
plot(tauex_fit, CQ_bar_all*tauex_fit,'k-','LineWidth',LineWidth_plot*3); %for all sites

xlim([0 0.35]);
ylim([0 3]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('Excess shear stress, $$\tau_{ex} = \tau - \tau_{it}$$ (Pa)','Interpreter','Latex');
ylabel('Normalized saltation flux, $$\hat{Q}=\frac{Q}{g/u_{*,it}}$$ (Pa)','Interpreter','Latex');
legend(SiteNames,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf,'PaperUnits','inches','PaperSize',[5 4.5],'PaperPosition',[0 0 5 4.5],'PaperPositionMode','Manual');
print([folder_Plots,'Q_tauex_norm.png'],'-dpng'); %for draft
print([folder_Plots,'Q_tauex_norm.tif'],'-dtiff','-r600'); %for publication

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY PLOTS FOR FLUX LAW PAPER %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT THETA MINUS MEAN THETA VERSUS TAU
figure(11); clf;

for i = 1:N_Sites
    %subplot(1,N_Sites,i); hold on;
    subplot('Position',[-0.22+0.315*i 0.13 0.25 0.82]); hold on;
    
    %plot theta values
    plot(tauRe_all{i}(Q_all{i}>0),theta_adjusted_all{i}(Q_all{i}>0),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);
    plot(tauRe_all{i}(Q_all{i}==0),theta_adjusted_all{i}(Q_all{i}==0),Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);  
    
    %format plot
    %xlim([0 0.4]);
    xlim([0 ceil(max(tauRe_all{i}/0.05))*0.05]);
    ylim([-110 80]);
    xlims = xlim; %get x limit information
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    if i==1
        ylabel('Wind angle, $$\theta$$ ($$^{\circ}$$)','Interpreter','Latex');
    end
    text(0.09*xlims(2), 73, panel_labels{i},'FontSize',9,'FontWeight','Bold');
    title(SiteNames{i});
    legend('saltation','no saltation','Location','SouthEast');
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 4],'PaperPosition',[0 0 7.3 4],'PaperPositionMode','Manual');
print([folder_Plots,'theta_tau.png'],'-dpng'); %for draft
print([folder_Plots,'theta_tau.tif'],'-dtiff','-r600'); %for publication

%% PLOT z/L VERSUS TAU
figure(12); clf;

for i = 1:N_Sites
    %perform plot
    %subplot(1,N_Sites,i); hold on;
    subplot('Position',[-0.23+0.32*i 0.13 0.25 0.82]); hold on;

    plot(tauRe_all{i}(Q_all{i}>0),zL_all{i}(Q_all{i}>0),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);
    plot(tauRe_all{i}(Q_all{i}==0),zL_all{i}(Q_all{i}==0),Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);
    
    %format plot
    legend('saltation','no saltation','Location','SouthEast');
%    xlim([0 0.4]);
    xlim([0 ceil(max(tauRe_all{i}/0.05))*0.05]);
    ylim([-0.5 0.1]);
    xlims = xlim; %get x limit information
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    if i==1
        ylabel('Stability parameter, $$z/L$$','Interpreter','Latex');
    end
    text(0.05*xlims(2), 0.08, panel_labels{i},'FontSize',9,'FontWeight','Bold');
    title(SiteNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 4],'PaperPosition',[0 0 7.3 4],'PaperPositionMode','Manual');
print([folder_Plots,'zL_tau.png'],'-dpng'); %for draft
print([folder_Plots,'zL_tau.tif'],'-dtiff','-r600'); %for publication

%% PLOT 30-MINUTE ZQ VS U*
figure(14); clf;

for i = 1:N_Sites
    %zq plot
    %subplot(1,N_Sites,i); hold on;
    subplot('Position',[-0.23+0.32*i 0.13 0.25 0.82]); hold on;
    
    %get values to plot
    ind_plot = find(Q_all{i}>0);
    ust_plot = ustRe_all{i}(ind_plot);
    sigma_ust_plot = sigma_ustRe_all{i}(ind_plot);
    zq_plot = zq_all{i}(ind_plot);
    sigma_zq_plot = sigma_zq_all{i}(ind_plot);
    
    %plot data
    plot(ust_plot,zq_plot,Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot);

    %plot error bars
    for j = 1:length(ust_plot)
        plot(ones(2,1)*ust_plot(j),zq_plot(j)+[-1 1]*sigma_zq_plot(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y-error
        plot(ust_plot(j)+[-1 1]*sigma_ust_plot(j),ones(2,1)*zq_plot(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x-error
    end
    
    %format plot
    %xlim([0.25 0.6]);
    xlim([floor(min(ust_plot)/0.05)*0.05 ceil(max(ust_plot)/0.05)*0.05]);
    ylim([0 0.14]);
    xlims = xlim; %get x limit information
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','log');
    xlabel('Shear velocity, $$u_{*}$$ (ms$$^{-1}$$)','Interpreter','Latex');
    if i==1
        ylabel('Saltation layer height, $$z_q$$ (m)','Interpreter','Latex');
    end
    text(0.03*range(xlims)+xlims(1), 0.135, panel_labels{i},'FontSize',9,'FontWeight','Bold');
    title(SiteNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 4],'PaperPosition',[0 0 7.3 4],'PaperPositionMode','Manual');
print([folder_Plots,'zq_ust_30minute.png'],'-dpng'); %for draft
print([folder_Plots,'zq_ust_30minute.tif'],'-dtiff','-r600'); %for publication


%% PLOT 30-MINUTE Q VS TAU
figure(15); clf;

for i = 1:N_Sites
    %subplot(1,N_Sites,i); hold on;
    subplot('Position',[-0.24+0.32*i 0.13 0.27 0.82]); hold on;

    %plot data
    plot(tauRe_all{i},Q_all{i},Markers_Field{i},'MarkerSize',MarkerSize_plot,'Color',Colors_Field{i});

    %plot error bars
    for j = 1:length(tauRe_all{i})
        plot(ones(2,1)*tauRe_all{i}(j),Q_all{i}(j)+[-1 1]*sigma_Q_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error
        plot(tauRe_all{i}(j)+[1 -1]*sigma_tauRe_all{i}(j),ones(2,1)*Q_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x error
    end
    
    %format plot
    %xlim([0 0.45]);
    xlim([0 ceil(max(tauRe_all{i}/0.05))*0.05]);
    ylim([0 65]);
    xlims = xlim;
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    if i==1
        ylabel('Saltation mass flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
    end
    text(0.05*xlims(2), 62, panel_labels{i},'FontSize',9,'FontWeight','Bold');
    title(SiteNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 4],'PaperPosition',[0 0 7.3 4],'PaperPositionMode','Manual');
print([folder_Plots,'Q_tau_30minute.png'],'-dpng'); %for draft
print([folder_Plots,'Q_tau_30minute.tif'],'-dtiff','-r600'); %for publication

%% PLOT TYPICAL STANDARD DEVIATIONS FOR Q's IN BINS
figure(16); clf; hold on;

%plot standard deviations
for i = 1:N_Sites
    plot(Q_bin_avg_stdmed_all{i}, Q_bin_std_stdmed_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot,'LineWidth',LineWidth_plot);
end
    
%plot median standard deviations
for i = 1:N_Sites
    plot(xlim,[Q_bin_stdmed_all(i) Q_bin_stdmed_all(i)],'Color',Colors_Field{i},'LineWidth',LineWidth_plot);
end

%format plot
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('Saltation mass flux, $$Q_i$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
ylabel('Saltation flux std. dev., $$SD_{Q_i}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
legend(SiteNames,'Location','NorthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf,'PaperUnits','inches','PaperSize',[5 4.5],'PaperPosition',[0 0 5 4.5],'PaperPositionMode','Manual');
print([folder_Plots,'Q_std.png'],'-dpng'); %for draft
print([folder_Plots,'Q_std.tif'],'-dtiff','-r600'); %for publication

%% PLOT Q VERSUS TAUEX BINNED DATA FOR ARXIV PAPER
figure(22); clf; hold on;

%plot binned values
for i = 1:N_Sites
    plot(tauex_fit_all{i},Q_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot);
end

%plot error bars
for i = 1:N_Sites
    for j = 1:length(tauex_fit_all{i})
        plot(ones(2,1)*tauex_fit_all{i}(j),Q_fit_all{i}(j)+[-1 1]*sigma_Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %y error
        plot(tauex_fit_all{i}(j)+[1 -1]*sigma_tauex_fit_all{i}(j),ones(2,1)*Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_plot); %x error
    end
end

%plot fit values
tauex_fit = [0 0.35];
for i = 1:N_Sites %for individual sites
    plot(tauex_fit, 1e3*ustit_linearfit_all(i)/g*CQ_all(i)*tauex_fit,'-','Color',Colors_Field{i},'LineWidth',LineWidth_plot);
end

xlim([0 0.35]);
ylim([0 65]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('Excess shear stress, $$\tau_{ex} = \tau - \tau_{it}$$ (Pa)','Interpreter','Latex');
ylabel('Saltation mass flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
legend(SiteNames,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf,'PaperUnits','inches','PaperSize',[6 4.5],'PaperPosition',[0 0 6 4.5],'PaperPositionMode','Manual');
print([folder_Plots,'Q_tauex.png'],'-dpng'); %for draft
print([folder_Plots,'Q_tauex.tif'],'-dtiff','-r600'); %for publication

%% PLOT CQ VERSUS TAUEX
figure(23); clf; hold on;

%plot binned values
for i = 1:N_Sites
    ind_plot = find(tauex_fit_all{i}>2*sigma_tauit_linearfit_all(i));
    plot(tauratio_fit_all{i}(ind_plot),CQ_fit_all{i}(ind_plot),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_plot);
    %errorbar(tauex_fit_all{i},CQ_fit_all{i},sigma_CQ_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
end

%plot error bars
for i = 1:N_Sites
    ind_plot = find(tauex_fit_all{i}>2*sigma_tauit_linearfit_all(i));
    for j = 1:length(ind_plot)
        plot(ones(2,1)*tauratio_fit_all{i}(ind_plot(j)),CQ_fit_all{i}(ind_plot(j))+[-1, 1]*sigma_CQ_fit_all{i}(ind_plot(j)),'Color',Colors_Field{i},'LineWidth',LineWidth_plot);
    end
end

%plot fit values
tauratio_fit = [1 4];
for i = 1:N_Sites
    plot(tauratio_fit, ones(2,1)*CQ_all(i),'Color',Colors_Field{i});
end

xlim([1 4]);
ylim([0 10]);
set(gca,'XMinorTick','On','YMinorTick','On','XScale','Log','Box','On');
xlabel('Dimensionless shear stress, $$\tau/\tau_{it}$$','Interpreter','Latex');
ylabel('Dimensionless saltation flux, $$\hat{C_{Q}}$$','Interpreter','Latex');
legend(SiteNames,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf,'PaperUnits','inches','PaperSize',[6 4.5],'PaperPosition',[0 0 6 4.5],'PaperPositionMode','Manual');
print([folder_Plots,'CQ_tauratio.png'],'-dpng'); %for draft
print([folder_Plots,'CQ_tauratio.tif'],'-dtiff','-r600'); %for publication