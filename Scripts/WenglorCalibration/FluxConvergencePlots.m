% %% INITIALIZATION
% %initialize
% clearvars;
% close all;
% 
% %% INFO FOR PLOTTING
% Markers_Field = {'s','d','o'};
% Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
% PlotFont = 14;
% 
% %% LOAD DATA AND FUNCTIONS
% %folders for loading data, saving data, and functions
% folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
% folder_Plots = '../../PlotOutput/WenglorCalibration/'; %folder containing plot output
% folder_Functions = '../Functions/'; %folder with functions
% 
% %paths for loading data
% LoadData_Path = strcat(folder_LoadData,'DataFullSubwindowAnalysis_30min_Restricted'); %path for 30 minute data
% 
% %load data
% load(LoadData_Path); %load window data
% 
% %load functions
% addpath(folder_Functions); %point MATLAB to location of functions

%%
%%%%%%%%%%%%%%%%%
% SUMMARY PLOTS %
%%%%%%%%%%%%%%%%%

%% plot median relative difference between Qsum and Qfit
figure(1); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, 100*median_Qrel_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$|Q_{fit}-Q_{sum}|/Q_{sum}$$ ($$\%$$)','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'Qsum_Qfit_T.png'],'-dpng'); %for draft


%% plot median Chi2nu of q(z) profile fit
figure(2); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, median_Chi2nu_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$\chi^{2}_{\nu}$$','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'zqfit_Chi2nu_T.png'],'-dpng'); %for draft


%% plot median relative uncertainy on q(z) values
figure(3); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, median_sigmaq_rel_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$\sigma_{q}/{q}$$','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'sigmaq_rel_T.png'],'-dpng'); %for draft


%% plot percent NaN for q(z) profile fits
figure(4); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, 100*fNaN_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('percent NaN for profile fit','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'zqfit_PercentNaN_T.png'],'-dpng'); %for draft


%% plot parameters for Q versus u^2 fit
figure(5); clf;

subplot('Position',[0.15 0.71 0.8 0.25]); hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, Chi2nu_Q_u2_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
%xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('$$\chi^{2}_{\nu}$$ for $$Q$$ vs $$u^2$$ fit','Interpreter','LaTeX');
legend(SiteNames,'Location','SouthWest');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log','YScale','Log');
set(gca,'FontSize',PlotFont);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% plot uth
subplot('Position',[0.15 0.39 0.8 0.25]); hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, uth_Q_u2_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
%xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('$$u_{th}$$ (m s$$^{-1}$$) for $$Q$$ vs $$u^2$$ fit','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% plot C
subplot('Position',[0.15 0.07 0.8 0.25]); hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, C_Q_u2_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('$$C$$ (g m$$^{-2}$$) for $$Q$$ vs $$u^u$$ fit','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[4 9],'PaperPosition',[0 0 5 9],'PaperPositionMode','Manual');
print([folder_Plots,'Q_u_fit_T.png'],'-dpng'); %for draft


%% plot standard deviation of zq values
figure(6); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s,std_zq_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
    plot(T_subwindow_s,std_zq_fW1_positive_analysis{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i})
    plot(T_subwindow_s,std_zq_fW1_nooutliers_analysis{i},['--',Markers_Field{i}],'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i})
end

% create legend
legend_items = cell(N_Sites*3,1);
for i = 1:N_Sites
    legend_items{i*3-2} = Sites{i};
    legend_items{i*3-1} = [Sites{i},' filtered'];
    legend_items{i*3} = [Sites{i},' no outliers'];
end 

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('std. dev. in $$z_q$$ from fit (m)','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'yscale','log','xscale','log');
set(gca,'FontSize',PlotFont);
h_legend = legend(legend_items,'Location','NorthEast');
set(h_legend,'FontSize',7);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 4],'PaperPosition',[0 0 7 4],'PaperPositionMode','Manual');
print([folder_Plots,'std_zq_T.png'],'-dpng'); %for draft

%% plot median of zq values
figure(7); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s,median_zq_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$z_q$$ from fit (m)','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'xscale','log');
set(gca,'FontSize',PlotFont);
legend(Sites,'Location','East');
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 4],'PaperPosition',[0 0 7 4],'PaperPositionMode','Manual');
print([folder_Plots,'median_zq_T.png'],'-dpng'); %for draft

%% plot median fW values
figure(8); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s,median_fW_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$N_{used}/N_{total}$$','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'xscale','log');
set(gca,'FontSize',PlotFont);
legend(Sites,'Location','East');
ylim([0 1]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 4],'PaperPosition',[0 0 7 4],'PaperPositionMode','Manual');
print([folder_Plots,'median_fW_T.png'],'-dpng'); %for draft

%% plot fraction of fW = 1 values
figure(9); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s,fraction_fW_full_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
end

% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('fraction $$N_{used}/N_{total} = 1$$','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'xscale','log');
set(gca,'FontSize',PlotFont);
legend(SiteNames,'Location','NorthEast');
ylim([0 1]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 4],'PaperPosition',[0 0 7 4],'PaperPositionMode','Manual');
print([folder_Plots,'fraction_fW_1_T.png'],'-dpng'); %for draft

% %%
% %%%%%%%%%%%%%%%%%%%%
% % TIME-SCALE PLOTS %
% %%%%%%%%%%%%%%%%%%%%
% 
% %% plot comparison of flux estimates
% for m = 1:N_T_subwindow
%     
%     % initialize figure for calculation timescale
%     figure(10+m); clf; hold on;
%     
%     % plot comparison of flux estimates
%     for i = 1:N_Sites
%         plot(Qsum_analysis{i}{m},Qfit_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
%     end
%     
%     % determine maximum Q
%     Qmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         Qmax(i) = ceil(max(max(Qsum_analysis{i}{m}),max(Qfit_analysis{i}{m}))/10)*10;
%     end
%     Qmax = max(Qmax);
%     
%     % plot 1-1 line
%     plot([0 Qmax],[0 Qmax],'k--');
%     
%     % format plot
%     xlim([0 Qmax]);
%     ylim([0 Qmax]);
%     xlabel('$$Q_{sum}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
%     ylabel('$$Q_{fit}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     legend(SiteNames,'Location','NorthWest');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'Q_comparison_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% end
% 
% 
% %% plot comparison of relative difference in flux estimates
% for m = 1:N_T_subwindow
%     
%     % initialize figure for calculation timescale
%     figure(10+N_T_subwindow+m); clf; hold on;
%     
%     % plot comparison of flux estimates
%     for i = 1:N_Sites
%         ind_nonzero = find(Qsum_analysis{i}{m}>=0);
%         plot(Qsum_analysis{i}{m}(ind_nonzero),100*Qrel_analysis{i}{m}(ind_nonzero),Markers_Field{i},'Color',Colors_Field{i})
%     end
%     
%     % determine maximum Q
%     Qmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         Qmax(i) = ceil(max(Qsum_analysis{i}{m}/10))*10;
%     end
%     Qmax = max(Qmax); 
%         
%     % format plot
%     xlim([0 Qmax]);
%     xlabel('$$Q_{sum}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
%     ylabel('$$|Q_{fit}-Q_{sum}|/Q_{sum}$$ ($$\%$$)','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     legend(SiteNames,'Location','NorthEast');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'Q_comparison_relative_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% end
% 
% 
% %% plot Q versus u for different timescales
% for m = 1:N_T_subwindow
%     
%     % initialize figure for calculation timescale
%     figure(10+2*N_T_subwindow+m); clf; hold on;
%   
%     % determine maximum u
%     umax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         umax(i) = ceil(max(ubar_analysis{i}{m}/10))*10;
%     end
%     umax = max(umax);
%     
%     % plot comparison of flux estimates
%     for i = 1:N_Sites
%         plot(ubar_analysis{i}{m},Qsum_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
%     end
% 
%     % plot fit
%     for i = 1:N_Sites
%         uth = uth_Q_u2_analysis{i}(m);
%         C = C_Q_u2_analysis{i}(m);
%         u_fit = linspace(uth,umax,100);
%         Q_fit = C*(u_fit.^2-uth.^2);
%         plot(u_fit,Q_fit,'Color',Colors_Field{i},'LineWidth',2)
%     end
%         
%     % determine maximum Q
%     Qmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         Qmax(i) = ceil(max(Qsum_analysis{i}{m}/10))*10;
%     end
%     Qmax = max(Qmax);  
%        
%     % format plot
%     ylim([0 Qmax]);
%     xlabel('$$u$$ (m s$$^{-1}$$)','Interpreter','LaTeX');
%     ylabel('$$Q_{sum}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     legend(SiteNames,'Location','NorthWest');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'Q_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% end
% 
% %% plot zq versus u for different timescales
% % create legend
% legend_items = cell(N_Sites*2,1);
% for i = 1:N_Sites
%     legend_items{i*2-1} = Sites{i};
%     legend_items{i*2} = [Sites{i},' filtered'];
% end 
% 
% for m = 1:N_T_subwindow
%     
%     % initialize figure for calculation timescale
%     figure(10+3*N_T_subwindow+m); clf; hold on;
%   
%     % determine maximum u
%     umax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         umax(i) = ceil(max(ubar_analysis{i}{m}/10))*10;
%     end
%     umax = max(umax);
%     
%     % plot comparison of flux height estimates
%     for i = 1:N_Sites
%         plot(ubar_analysis{i}{m},zq_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
%         plot(ubar_analysis{i}{m}(fW_analysis{i}{m}==1),zq_analysis{i}{m}(fW_analysis{i}{m}==1),Markers_Field{i},'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i})
%     end
%        
%     % determine minimum and maximum zq for plotting
%     zqmin = zeros(N_Sites,1); %compute for each site, then take min of all sites
%     zqmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         zqmin(i) = floor(min(zq_analysis{i}{m}(zq_analysis{i}{m}>0))*100)/100;
%         zqmax(i) = ceil(max(zq_analysis{i}{m}*100))/100;
%     end
%     zqmin = min(zqmin);
%     zqmax = max(zqmax);
%        
%     % format plot
%     ylim([0 zqmax]);
%     xlabel('$$u$$ (m s$$^{-1}$$)','Interpreter','LaTeX');
%     ylabel('$$z_{q}$$ (m)','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     h_legend = legend(legend_items,'Location','NorthEast');
%     set(h_legend,'FontSize',8);
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'zq_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
%     
%     %change plot to logscale, eliminate unfiltered values
%     ylim([zqmin zqmax]);
%     set(gca,'yscale','log');
%     
%     % plot comparison of flux height estimates
%     cla;
%     for i = 1:N_Sites
%         plot(ubar_analysis{i}{m}(fW_analysis{i}{m}==1),zq_analysis{i}{m}(fW_analysis{i}{m}==1),Markers_Field{i},'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i})
%     end
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'zq_u_',num2str(T_subwindow_s(m)),'_log.png'],'-dpng'); %for draft
% end
% 
% %% plot profile for zqmax - Oceano, T = 30 min, unfiltered
% 
% %get values for plotting
% i = 3; m = N_T_subwindow;
% ind_zqmax = find(zq_analysis{i}{m}==max(zq_analysis{i}{m})); %get index for max zq
% qbar = qbar_analysis{i}{m}{ind_zqmax}; %get qbar profile
% sigma_qbar = sigma_qbar_analysis{i}{m}{ind_zqmax}; %get sigma_qbar profile
% zW = zW_analysis{i}{m}{ind_zqmax}; %get zW profile
% sigma_zW = sigma_zW_analysis{i}{m}{ind_zqmax}; %get sigma_zW profile
% 
% %create figure
% figure(11+4*N_T_subwindow); clf; hold on;
% plot(zW,qbar,'bo');
% for j = 1:length(zW)
%     plot(zW(j)+[-1 1]*sigma_zW(j),qbar(j)*[1 1],'b');
%     plot(zW(j)*[1 1],qbar(j)+[-1 1]*sigma_qbar(j),'b');
% end
% 
% %organize figure
% xlabel('$$z$$ (m)','Interpreter','LaTeX');
% ylabel('$$q$$ (g m$$^{-1}$$ s$$^{-2}$$)','Interpreter','LaTeX');
% title([Sites{i},', ',datestr(StartTime_analysis{i}{m}(ind_zqmax),'yyyy-mm-dd HH:MM'),', T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
% print([folder_Plots,'q_z_profile_',Sites{i},'_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% 
% %% plot profile for zqmax - Oceano, T = 2 s, unfiltered
% 
% %get values for plotting
% i = 3; m = find(T_subwindow_s == 2);
% ind_calc = find(fW_analysis{i}{m}==1);
% ind_zqmax = find(zq_analysis{i}{m}==max(zq_analysis{i}{m}(ind_calc))); %get index for max zq
% qbar = qbar_analysis{i}{m}{ind_zqmax}; %get qbar profile
% sigma_qbar = sigma_qbar_analysis{i}{m}{ind_zqmax}; %get sigma_qbar profile
% zW = zW_analysis{i}{m}{ind_zqmax}; %get zW profile
% sigma_zW = sigma_zW_analysis{i}{m}{ind_zqmax}; %get sigma_zW profile
% 
% %create figure
% figure(12+4*N_T_subwindow); clf; hold on;
% plot(zW,qbar,'bo');
% for j = 1:length(zW)
%     plot(zW(j)+[-1 1]*sigma_zW(j),qbar(j)*[1 1],'b');
%     plot(zW(j)*[1 1],qbar(j)+[-1 1]*sigma_qbar(j),'b');
% end
% 
% %organize figure
% xlabel('$$z$$ (m)','Interpreter','LaTeX');
% ylabel('$$q$$ (g m$$^{-1}$$ s$$^{-2}$$)','Interpreter','LaTeX');
% title([Sites{i},', ',datestr(StartTime_analysis{i}{m}(ind_zqmax),'yyyy-mm-dd HH:MM'),', T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
% print([folder_Plots,'q_z_profile_',Sites{i},'_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% 
% 
% %% plot fW (fraction of Wenglors used) versus u for different timescales
% for m = 1:N_T_subwindow
%     
%     % initialize figure for calculation timescale
%     figure(12+4*N_T_subwindow+m); clf; hold on;
%   
%     % determine maximum u
%     umax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         umax(i) = ceil(max(ubar_analysis{i}{m}/10))*10;
%     end
%     umax = max(umax);
%     
%     % plot comparison of N_used / N_total
%     for i = 1:N_Sites
%         plot(ubar_analysis{i}{m},fW_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
%     end
%        
%     % format plot
%     ylim([0 1]);
%     xlabel('$$u$$ (m/s)','Interpreter','LaTeX');
%     ylabel('$$N_{used} / N_{total}$$','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     legend(SiteNames,'Location','NorthWest');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
%     print([folder_Plots,'fW_used_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% end