clearvars;

%load flux stress window data
folder_WindowData = '../../AnalysisData/Windowing/'; %folder for loading window data
load(strcat(folder_WindowData,'FluxLawWindows_Analysis'));
N_Sites = length(Sites);

%set plotting info
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots
PlotFont = 14;
Markers_Field = {'s','d','o','<','>'};
Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};

%determine z0 based on fQ = 0 windows
ind_notransport = find(fQ_all{3}<0.05);
z0 = mean(zs_all{3}(ind_notransport));
sigma_z0 = std(zs_all{3}(ind_notransport))/sqrt(length(ind_notransport));

%plot zs versus fQ
figure(1); clf;
for i = 1:N_Sites
    hold on;
    ind_plot = intersect(find(fQ_all{i}>=0.05),find(fQ_all{i}<=0.95));
    plot(fQ_all{i}(ind_plot),zs_all{i}(ind_plot),Markers_Field{i},'Color',Colors_Field{i})
end
set(gca,'yscale','log');

xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective roughness ht., $$z_{s}$$','interpreter','latex');
legend(SiteNames,'Location','SouthEast');

set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'zs_fQ.png'],'-dpng');