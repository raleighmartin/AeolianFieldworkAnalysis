%% initialize
clearvars;

%% plot calibration coefficient for sample day
Site_SamplePlot = 'Oceano';
Date_SamplePlot = datetime(2015,6,3);

%% plotting information
PlotFont = 10;
LineWidth_Plot = 1;
Marker_Site = {'s','d','o'};
Color_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/BSNE/'; %folder for mean grain size data
folder_SaltationData = '../../AnalysisData/Windowing/'; %folder for saltation flux data
SaltationFluxData_Path = strcat(folder_SaltationData,'DataWindowCalcs_30min_Restricted'); %path for loading saltation data
folder_Plots = '../../PlotOutput/Methods/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%load BSNE flux data
BSNEData_Path = strcat(folder_ProcessedData,'FluxBSNE_',Site_SamplePlot);
load(BSNEData_Path);

%load Wenglor saltation flux data
load(SaltationFluxData_Path);

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% get data for BSNEs for interval
ind_BSNE = find([FluxBSNE.Date]==Date_SamplePlot);
StartTime_BSNE = [FluxBSNE(ind_BSNE).StartTime];
EndTime_BSNE = [FluxBSNE(ind_BSNE).EndTime];
MidTime_BSNE = mean([StartTime_BSNE; EndTime_BSNE]);
N_BSNE = length(ind_BSNE);
Q_BSNE = zeros(N_BSNE,1);
sigma_Q_BSNE = zeros(N_BSNE,1);
for i=1:N_BSNE
    Q_BSNE(i) = FluxBSNE(ind_BSNE(i)).Q.Q;
    sigma_Q_BSNE(i) = FluxBSNE(ind_BSNE(i)).Q.sigma_Q;
end

%% get data for Wenglors for interval
ind_Site = find(strcmp(SiteNames,Site_SamplePlot));
ind_Time = zeros(N_BSNE,1);
for i=1:N_BSNE
    ind_Time(i) = find(StartTimes_all{ind_Site}>=StartTime_BSNE(i), 1 );
end
N_zW = length(zW_all{ind_Site}{ind_Time(1)});

zW_matrix = zeros(N_BSNE,N_zW);
Cqn_matrix = zeros(N_BSNE,N_zW);
for i = 1:N_BSNE
    zW_matrix(i,:) = zW_all{ind_Site}{ind_Time(i)};
    Cqn_matrix(i,:) = Cqnbar_all{ind_Site}{ind_Time(i)}';
end

%% generate plot
figure(1); clf; hold on;

subplot(3,1,1:2); hold on;

Time_Plot = [StartTime_BSNE, EndTime_BSNE];
[Time_Plot,ind_Plot] = sort(Time_Plot);

for i = 1:N_BSNE
    Cqn_Plot = ([Cqn_matrix(i,:), Cqn_matrix(i,:)]);
    Cqn_Plot = Cqn_Plot(ind_Plot);
    plot(Time_Plot,Cqn_Plot);
    %plot(MidTime_BSNE,Cqn_matrix(i,:));
end
ylabel('HF sensor calibration factor, $$C_{qn,i}$$ (g m$$^{-2}$$)','Interpreter','Latex');

legend_items = cell(N_zW,1);
for i = 1:N_zW
    legend_items{i} = ['z_{HF,',int2str(i),'} = ',num2str(round(zW_matrix(1,i)*100,1)),' cm'];
end
legend(legend_items,'Location','NorthWest');
text(StartTime_BSNE(2),8.5,'(a)','FontSize',PlotFont);
set(gca,'yscale','log','box','on','FontSize',PlotFont);
ylim([7e-1 1e1]);

subplot(3,1,3); hold on;
plot(MidTime_BSNE,Q_BSNE,'ko');

for i = 1:N_BSNE
    plot([MidTime_BSNE(i) MidTime_BSNE(i)],Q_BSNE(i)*[1 1]+sigma_Q_BSNE(i)*[-1 1],'k');
    plot([StartTime_BSNE(i) EndTime_BSNE(i)],Q_BSNE(i)*[1 1],'k');
end
ylabel('LF trap saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
text(StartTime_BSNE(1),35,'(b)','FontSize',PlotFont);
set(gca,'box','on','FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[5 7],'PaperPosition',[0 0 5 7],'PaperPositionMode','Manual');
print([folder_Plots,'CalibrationCoefficient_SamplePlot.png'],'-dpng');
