%% INITIALIZATION
%initialize
clearvars;
close all;

%% INFO FOR PLOTTING
Markers_Field = {'s','d','o'};
%Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.125]*0.8};
PlotFont = 14;

%% LOAD DATA AND FUNCTIONS
%folders for loading data, saving data, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Plots = '../../PlotOutput/Methods/'; %folder containing plot output
folder_Functions = '../Functions/'; %folder with functions

%paths for loading data
LoadData_Path = strcat(folder_LoadData,'DataFullSubwindowAnalysis_30min_Restricted'); %path for 30 minute data

%load data
load(LoadData_Path); %load window data

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS - PROFILE FITTING %
% %%%%%%%%%%%%%%%%%%%%%%%%%

%% figure on Wenglor heights with zeros
figure(1); clf;

% boxplot of fraction of unused Wenglor heights in analysis
subplot('Position',[0.07 0.15 0.5 0.8]); hold on;

%plot symbols
for i = 1:N_Sites
    plot(T_subwindow_s,1-median_f_zW_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
end
xlim([0.9*min(T_subwindow_s) max(T_subwindow_s)]);
set(gca,'xscale','log');
xTicks = get(gca, 'xtick');
xTickLabels = get(gca, 'xticklabel');

%generate and plot boxplot
for i = 1:N_Sites
    %create vector for boxplot
    f_boxplot = [];
    T_boxplot = [];
    for j = 1:N_T_subwindow
        ind_positive_Qsum = find(Qsum_all_analysis{i}{j}>0);
        f = 1-f_zW_analysis{i}{j}(ind_positive_Qsum);
        f_boxplot = [f_boxplot; f];
        T_boxplot = [T_boxplot; T_subwindow_s(j)*ones(size(f))];
    end
        
    %create boxplot
    boxplot(f_boxplot,T_boxplot,'Colors',Colors_Field{i},'Symbol','','Positions',T_subwindow_s,'Widths',T_subwindow_s/4);
end

% format plot
xlabel('sampling time scale, $$\Delta t$$ (s)','Interpreter','LaTeX');
ylabel('frac. of HF sensor hts. with zero flux');
set(gca,'XTick',xTicks,'XTickLabel',{'1','10','100','1000'},'XMinorTick','On');
set(gca,'YMinorTick','On','Box','On');
set(gca,'xscale','log');
set(gca,'FontSize',PlotFont);
legend(SiteNames,'Location','NorthEast');
ylim([0 1]);
xlim([0.9*min(T_subwindow_s) max(T_subwindow_s)]);
text(0.8*mean(T_subwindow_s(1:2)),0.95,'(a)','FontSize',PlotFont);

% % plot median fraction of unused Wenglor heights in analysis
% subplot('Position',[0.07 0.15 0.5 0.8]); hold on;
% for i = 1:N_Sites
% %    plot(T_subwindow_s,1-median_f_zW_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
%     errorbar(T_subwindow_s,1-median_f_zW_analysis{i},...
%         (upperquartile_f_zW_analysis{i}-median_f_zW_analysis{i}),...
%         (median_f_zW_analysis{i}-lowerquartile_f_zW_analysis{i}),...
%         Markers_Field{i},'Color',Colors_Field{i})
% end
% 
% % format plot
% xlabel('sampling time scale, $$\Delta t$$ (s)','Interpreter','LaTeX');
% ylabel('frac. of HF sensor hts. with zero flux');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'xscale','log');
% set(gca,'FontSize',PlotFont);
% %legend(SiteNames,'Location','NorthEast');
% ylim([0 1]);
% xlim([0.9*min(T_subwindow_s) max(T_subwindow_s)]);
% text(0.8*mean(T_subwindow_s(1:2)),0.95,'(a)','FontSize',PlotFont);

% plot example of number of Wenglors with zeros
subplot('Position',[0.65 0.15 0.34 0.8]); hold on;
for i = 1:N_Sites
    plot(u_f0_analysis,f0_analysis{i}, Markers_Field{i},'Color',Colors_Field{i})
end
% format plot
xlabel('wind speed, $$u$$ (m s$$^{-1}$$)','Interpreter','LaTeX');
ylabel('frac. of time intervals with zero flux');
set(gca,'YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);
legend(SiteNames,'Location','NorthEast');
ylim([0 1]);
text(u_f0_analysis(1)+0.25,0.95,'(b)','FontSize',PlotFont);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[10 4],'PaperPosition',[0 0 10 4],'PaperPositionMode','Manual');
print([folder_Plots,'f_zW_T.png'],'-dpng');
print([folder_Plots,'f_zW_T.tif'],'-dtiff');

% 
% 
% %% plot fraction of profiles with full number of Wenglor heights used
% figure(2); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s,fraction_f_zW_full_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
% end
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('fraction of Wenglor profiles that are full');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'xscale','log');
% set(gca,'FontSize',PlotFont);
% legend(SiteNames,'Location','SouthEast');
% ylim([0 1]);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'fraction_f_zW_full_T.png'],'-dpng');
% 
% 
% %% plot fraction of intervals with outlier values for q(z) fits
% figure(3); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, f_outlier_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
%     plot(T_subwindow_s, f_full_outlier_analysis{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i});
% end
% 
% % create legend
% legend_items = cell(N_Sites*2,1);
% for i = 1:N_Sites
%     legend_items{i*2-1} = [SiteNames{i},' all profiles'];
%     legend_items{i*2} = [SiteNames{i},' full profiles'];
% end 
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('fraction of profiles with outlier $$z_q$$','Interpreter','LaTeX');
% h_legend = legend(legend_items,'Location','SouthEast');
% set(h_legend,'FontSize',10);
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On')
% set(gca,'XScale','Log','YScale','Log');
% set(gca,'FontSize',PlotFont);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[7.5 4],'PaperPosition',[0 0 7.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'f_zq_outlier_T.png'],'-dpng');
% 
% %% plot median Chi2nu for q(z) profile fit
% figure(4); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, median_Chi2nu_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
% end
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('median $$\chi^{2}_{\nu}$$ for profile fit','Interpreter','LaTeX');
% legend(SiteNames,'Location','SouthWest');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
% set(gca,'FontSize',PlotFont);
% ylims = ylim;
% ylim([0 ylims(2)]);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
% print([folder_Plots,'Chi2nu_profilefit_T.png'],'-dpng');
% 
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS - SALTATION HEIGHT %
% %%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot median of zq values
figure(5); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s,median_zq_real_analysis{i},[Markers_Field{i},'-'],'Color',Colors_Field{i})
    %plot(T_subwindow_s,mean_zq_real_analysis{i},[Markers_Field{i},'-'],'Color',Colors_Field{i})
    plot(T_subwindow_s,median_zq_full_analysis{i},[Markers_Field{i},'--'],'Color',Colors_Field{i})
end

% include error bars as median of sigma_zq values
for i = 1:N_Sites
    for m = 1:N_T_subwindow
        plot(T_subwindow_s(m)*[1 1],median_zq_real_analysis{i}(m)+median_sigma_zq_real_analysis{i}(m)*[-1 1],'Color',Colors_Field{i})
        %plot(T_subwindow_s(m)*[1 1],mean_zq_real_analysis{i}(m)+mean_sigma_zq_real_analysis{i}(m)*[-1 1],'Color',Colors_Field{i})
        plot(T_subwindow_s(m)*[1 1],median_zq_full_analysis{i}(m)+median_sigma_zq_full_analysis{i}(m)*[-1 1],'Color',Colors_Field{i})
    end
end

% create legend
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    legend_items{i*2-1} = [SiteNames{i},' all profiles'];
    legend_items{i*2} = [SiteNames{i},' full profiles'];
end 

% format plot
xlabel('subsampling time scale, $$\Delta t$$ (s)','Interpreter','LaTeX');
ylabel('median saltation height from exp. profile fit, $$z_q \pm \sigma_{z_q}$$ (m)','Interpreter','LaTeX');
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont,'LooseInset', get(gca,'TightInset'));
legend(legend_items,'Location','SouthEast','FontSize',11);
xlim([0.9*min(T_subwindow_s) max(T_subwindow_s)]);
ylim([0 0.14]);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 6],'PaperPosition',[0 0 7 6],'PaperPositionMode','Manual');
print([folder_Plots,'median_zq_T.png'],'-dpng');
print([folder_Plots,'median_zq_T.tif'],'-dtiff');

% 
% 
% %% plot median of sigma_zq values
% figure(6); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s,median_sigma_zq_full_analysis{i}./median_zq_full_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
% end
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('median $$\sigma_{z_q}/z_q$$ from fit','Interpreter','LaTeX');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'xscale','log','yscale','log');
% set(gca,'FontSize',PlotFont);
% legend(SiteNames,'Location','East');
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
% print([folder_Plots,'median_sigma_zq_T.png'],'-dpng');
% 
% 
% %% plot standard deviation of zq values
% figure(7); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s,std_zq_real_analysis{i},Markers_Field{i},'Color',Colors_Field{i})
%     plot(T_subwindow_s,std_zq_full_analysis{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i})
%     plot(T_subwindow_s,std_zq_nooutlier_analysis{i},['--',Markers_Field{i}],'Color',Colors_Field{i})
%     plot(T_subwindow_s,std_zq_full_nooutlier_analysis{i},['--',Markers_Field{i}],'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i})
% end
% 
% % create legend
% legend_items = cell(N_Sites*4,1);
% for i = 1:N_Sites
%     legend_items{i*4-3} = [SiteNames{i},' all profiles'];
%     legend_items{i*4-2} = [SiteNames{i},' full profiles'];
%     legend_items{i*4-1} = [SiteNames{i},' no outlier profiles'];
%     legend_items{i*4} = [SiteNames{i},' full profiles, no outliers'];
% end 
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('std. dev. in $$z_q$$ from fit (m)','Interpreter','LaTeX');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% set(gca,'yscale','log','xscale','log');
% set(gca,'FontSize',PlotFont);
% h_legend = legend(legend_items,'Location','NorthEast');
% set(h_legend,'FontSize',7);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 6],'PaperPosition',[0 0 6 6],'PaperPositionMode','Manual');
% print([folder_Plots,'std_zq_T.png'],'-dpng');
% 
% %%
% % %%%%%%%%%%%%%%%%%%%%
% % PLOTS - TOTAL FLUX %
% % %%%%%%%%%%%%%%%%%%%%
% 
% %% plot median relative difference between Qsum and Qfit
% figure(8); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, median_Qrel_real_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
% end
% 
% % format plot
% xlabel('subsampling time scale, $$\Delta t$$ (s)','Interpreter','LaTeX');
% ylabel('median $$(Q_{fit}-Q_{sum})/Q_{sum}$$','Interpreter','LaTeX');
% legend(SiteNames,'Location','NorthEast');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log','LooseInset', get(gca,'TightInset'));
% set(gca,'FontSize',PlotFont);
% %ylims = ylim;
% %ylim([0 ylims(2)]);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'Qsum_Qfit_T.png'],'-dpng');
% 
% %% plot median ratio between sigma_Qsum and sigma_Qfit
% figure(9); clf; hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, median_sigma_Qratio_real_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
% end
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('median $$\sigma_{Q_{fit}}/\sigma_{Q_{sum}}$$','Interpreter','LaTeX');
% legend(SiteNames,'Location','South');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
% set(gca,'FontSize',PlotFont);
% ylims = ylim;
% ylim([0 ylims(2)]);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'sigma_Qfit_Qsum_ratio_T.png'],'-dpng');

%% plot mean of Q values
figure(10); clf; hold on;

%subplot(2,1,1); hold on;
subplot('Position',[0.12 0.48 0.86 0.5]); hold on;

%for i = 1
for i = 1:N_Sites
    plot(T_subwindow_s,mean_Qfit_calc_analysis{i},[Markers_Field{i},'-'],'Color',Colors_Field{i})
    %plot(T_subwindow_s,mean_Qsum_calc_analysis{i},[Markers_Field{i},'--'],'Color',Colors_Field{i})
    plot(T_subwindow_s,mean_Qsum_all_analysis{i},[Markers_Field{i},'--'],'Color',Colors_Field{i})
end

% % include error bars as mean of sigma_zq values
% for i = 1:N_Sites
%     for m = 1:N_T_subwindow
%         plot(T_subwindow_s(m)*[1 1],mean_Qfit_calc_analysis{i}(m)+SD_Qfit_calc_analysis{i}(m)*[-1 1],'Color',Colors_Field{i})
%         %plot(T_subwindow_s(m)*[1 1],mean_Qsum_calc_analysis{i}(m)+SD_Qsum_calc_analysis{i}(m)*[-1 1],'Color',Colors_Field{i})
%         plot(T_subwindow_s(m)*[1 1],mean_Qsum_all_analysis{i}(m)+SD_Qsum_all_analysis{i}(m)*[-1 1],'Color',Colors_Field{i})
%     end
% end

% create legend
legend_items = cell(N_Sites*2,1);
for i = 1:N_Sites
    legend_items{i*2-1} = [SiteNames{i},' Q_{fit}'];
    %legend_items{i*3-1} = [SiteNames{i},' Q_{sum}'];
    legend_items{i*2} = [SiteNames{i},' Q_{sum}'];
end 

% format plot
%xlabel('subsampling time scale, $$\Delta t$$ (s)','Interpreter','LaTeX');
ylabel('mean saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont,'LooseInset', get(gca,'TightInset'));
legend(legend_items,'Location','East','FontSize',11);
xlim([0.9*min(T_subwindow_s) max(T_subwindow_s)]);
text(1.5e0,34,'(a)','FontSize',PlotFont);
% ylims = ylim;
% ylim([0 ylims(2)]);


% plot uncertainty
%subplot(2,1,2); hold on;
%subplot('Position',[0.15 0.39 0.8 0.25]); hold on;
subplot('Position',[0.12 0.08 0.86 0.32]); hold on;

%for i = 1
for i = 1:N_Sites
    plot(T_subwindow_s,mean_sigma_Qfit_calc_analysis{i},[Markers_Field{i},'-'],'Color',Colors_Field{i})
    %plot(T_subwindow_s,mean_sigma_Qsum_calc_analysis{i},[Markers_Field{i},'--'],'Color',Colors_Field{i})
    plot(T_subwindow_s,mean_sigma_Qsum_positive_analysis{i},[Markers_Field{i},'--'],'Color',Colors_Field{i})
end

% format plot
xlabel('subsampling time scale, $$\Delta t$$ (s)','Interpreter','LaTeX');
ylabel('mean flux uncertainty, $$\sigma_{Q}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
set(gca,'xscale','log','yscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont,'LooseInset', get(gca,'TightInset'));
%legend(legend_items,'Location','East','FontSize',11);
xlim([0.9*min(T_subwindow_s) max(T_subwindow_s)]);
ylim([1e-1 1e4]);
text(1.5e0,5e3,'(b)','FontSize',PlotFont);

% print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 8],'PaperPosition',[0 0 6 8],'PaperPositionMode','Manual');
print([folder_Plots,'median_Q_T.png'],'-dpng');
print([folder_Plots,'median_Q_T.tif'],'-dtiff');

% %%
% % %%%%%%%%%%%%%%%%%%
% % PLOTS - FLUX LAW %
% % %%%%%%%%%%%%%%%%%%
% 
% %% plot parameters for Q versus u^2 fit
% figure(10); clf;
% 
% subplot('Position',[0.15 0.71 0.8 0.25]); hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, Chi2nu_Q_u2_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
% end
% 
% % format plot
% ylabel('$$\chi^{2}_{\nu}$$ for $$Q$$ vs $$u^2$$ fit','Interpreter','LaTeX');
% legend(SiteNames,'Location','SouthWest');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log','YScale','Log');
% set(gca,'FontSize',PlotFont);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % plot uth
% subplot('Position',[0.15 0.39 0.8 0.25]); hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, uth_Q_u2_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
% end
% 
% % format plot
% ylabel('$$u_{th}$$ (m s$$^{-1}$$) for $$Q$$ vs $$u^2$$ fit','Interpreter','LaTeX');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
% set(gca,'FontSize',PlotFont);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % plot C
% subplot('Position',[0.15 0.07 0.8 0.25]); hold on;
% for i = 1:N_Sites
%     plot(T_subwindow_s, C_Q_u2_analysis{i},Markers_Field{i},'Color',Colors_Field{i});
% end
% 
% % format plot
% xlabel('averaging time scale, $$T_{HF}$$ (s)','Interpreter','LaTeX');
% ylabel('$$C$$ (g m$$^{-2}$$) for $$Q$$ vs $$u^u$$ fit','Interpreter','LaTeX');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
% set(gca,'FontSize',PlotFont);
% xlim([min(T_subwindow_s) max(T_subwindow_s)]);
% 
% % print plot
% set(gcf,'PaperUnits','inches','PaperSize',[4 9],'PaperPosition',[0 0 5 9],'PaperPositionMode','Manual');
% print([folder_Plots,'Q_u_fit_T.png'],'-dpng');
% 
%
% % %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % INDIVIDUAL PLOTS BY TIME-SCALE %
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%     set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
%     print([folder_Plots,'Q_comparison_',num2str(T_subwindow_s(m)),'.png'],'-dpng');
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
%         plot(Qsum_analysis{i}{m},Qrel_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
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
%     ylabel('$$|Q_{fit}-Q_{sum}|/Q_{sum}$$','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     legend(SiteNames,'Location','SouthEast');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'YScale','log');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
%     print([folder_Plots,'Q_comparison_relative_',num2str(T_subwindow_s(m)),'.png'],'-dpng');
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
%         umax(i) = ceil(max(ubar_full_notoutlier_analysis{i}{m}/10))*10;
%     end
%     umax = max(umax);
%     
%     % plot comparison of flux estimates
%     for i = 1:N_Sites
%         plot(ubar_full_notoutlier_analysis{i}{m},Qsum_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
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
%     set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
%     print([folder_Plots,'Q_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng');
% end
% 
% 
% %% plot zq versus u for different timescales
% 
% % create legend
% legend_items = cell(N_Sites*2,1);
% for i = 1:N_Sites
%     legend_items{i*2-1} = [SiteNames{i},' all profiles'];
%     legend_items{i*2} = [SiteNames{i},' full profiles'];
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
%         plot(ubar_analysis{i}{m}, zq_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i});
%         plot(ubar_full_analysis{i}{m},full_zq_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i},'MarkerFaceColor',Colors_Field{i});
%     end
%        
%     % determine minimum and maximum zq for plotting
%     zqmin = zeros(N_Sites,1); %compute for each site, then take min of all sites
%     zqmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         zqmin(i) = floor(min(zq_analysis{i}{m}(zq_analysis{i}{m}>0))*100)/100;
%         zqmax(i) = ceil(max(zq_analysis{i}{m}*100))/100;
% %         zqmin(i) = floor(min(full_notoutlier_zq_analysis{i}{m}(full_notoutlier_zq_analysis{i}{m}>0))*100)/100;
% %         zqmax(i) = ceil(max(full_notoutlier_zq_analysis{i}{m}*100))/100;
%     end
%     zqmin = min(zqmin);
%     zqmax = max(zqmax);
%        
%     % format plot
%     %ylim([0 zqmax]);
%     xlabel('$$u$$ (m s$$^{-1}$$)','Interpreter','LaTeX');
%     ylabel('$$z_{q}$$ (m)','Interpreter','LaTeX');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     h_legend = legend(legend_items,'Location','NorthEast');
%     set(h_legend,'FontSize',10);
%     %legend(SiteNames,'Location','SouthEast');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
%     print([folder_Plots,'zq_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
% end
%
% 
% %% plot fW (fraction of Wenglors used) versus u for different timescales
% for m = 1:N_T_subwindow
%     
%     % initialize figure for calculation timescale
%     figure(10+4*N_T_subwindow+m); clf; hold on;
%   
%     % determine maximum u
%     umax = zeros(N_Sites,1); %compute for each site, then take max of all sites
%     for i = 1:N_Sites
%         umax(i) = ceil(max(ubar_all_analysis{i}{m}/10))*10;
%     end
%     umax = max(umax);
%     
%     % plot comparison of N_used / N_total
%     for i = 1:N_Sites
%         plot(ubar_all_analysis{i}{m},f_zW_analysis{i}{m},Markers_Field{i},'Color',Colors_Field{i})
%     end
%        
%     % format plot
%     ylim([0 1]);
%     xlabel('mean wind speed, $$u$$ (m/s)','Interpreter','LaTeX');
%     ylabel('fraction of Wenglor heights used in profile fit');
%     title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
%     legend(SiteNames,'Location','NorthWest');
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     set(gca,'FontSize',PlotFont);
%     
%     % print plot
%     set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
%     print([folder_Plots,'fW_used_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng');
% end