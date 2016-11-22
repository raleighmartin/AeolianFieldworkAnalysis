%% initialize
clearvars;

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/Thresholds/'; %folder for threshold data
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data

%% load data
load(strcat(folder_AnalysisData,'ThresholdAnalysisData')); %load analysis windows
load(strcat(folder_GrainSizeData,'MeanGrainSize')); %load grain size data

%% tauit and tauft for d50
tauit_d50 = [0.072; 0.073; 0.054];
sigma_tauit_d50 = [0.006; 0.004; 0.011];
tauft_d50 = [0.143; 0.138; 0.104];
sigma_tauft_d50 = [0.010; 0.008; 0.018];

%% plotting info
PlotFont = 14;
PlotMarkers_Site = {'s','d','o'};
PlotColors_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
PlotMarkers_Sensitivity = {'s','d','o','<','>'};
PlotColors_Sensitivity = {[0.6473 0.7456 0.4188],[0.2116 0.1898 0.5777],[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};

%% Get info about sites
N_Sites = length(Sites);

%additional plot for alternative values
figure(1); clf; hold on;
for i = 1:N_Sites
    plot([0 0]+i/2,tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'-k','LineWidth',2); %plot intercept SE dummy values for legend
    plot(4+[0 0]+i/2,tauft_d50(i)+sigma_tauft_d50(i)*[-1 1],'-.k','LineWidth',2); %plot ft SE dummy values for legend
    plot(4+[0 0]+i/2,tauit_d50(i)+sigma_tauit_d50(i)*[-1 1],'--k','LineWidth',2); %plot it SE dummy values for legend
    plot([0 0]+i/2,tauit_intercept(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i}); %plot intercept average values
    plot(4+[0 0]+i/2,tauft_d50(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i}); %plot ft average values
    plot(4+[0 0]+i/2,tauit_d50(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i}); %plot it average values
    plot([0 0]+i/2,tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',2); %plot intercept SE values
    plot(4+[0 0]+i/2,tauft_d50(i)+sigma_tauft_d50(i)*[-1 1],'-.','Color',PlotColors_Site{i},'LineWidth',2); %plot ft SE values
    plot(4+[0 0]+i/2,tauit_d50(i)+sigma_tauit_d50(i)*[-1 1],'--','Color',PlotColors_Site{i},'LineWidth',2); %plot it SE values
end
ylims = ylim; %get min and max y-values

%organize plot
legend('\tau_{th,flux}','\tau_{ft,d_{50}}','\tau_{it,d_{50}}','Location','SouthWest');
xlim([0 6]);
ylabel('alternative threshold stress estimates (Pa)','interpreter','latex');
set(gca,'YTickLabel',{''},'YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

