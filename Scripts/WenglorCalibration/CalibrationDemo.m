%% INITIALIZATION
%initialize
clearvars;
close all;

%% PARAMETERS AND INPUTS
%info about site and time intervals for plots
Site_Demo = 'Jericoacoara';
StartTime_Demo_BSNE = datetime(2014,11,14,15,26,0); %start time for BSNE profile demo
EndTime_Demo_BSNE = datetime(2014,11,14,16,26,0); %end time for BSNE profile demo
%StartTime_Demo_Wenglor = datetime(2014,11,14,15,26,0); %start time for Wenglor profile demo
%EndTime_Demo_Wenglor = datetime(2014,11,14,15,56,0); %end time for Wenglor profile demo
StartTime_Demo_Subwindow = datetime(2014,11,14,15,26,0); %start time for Wenglor profile demo
EndTime_Demo_Subwindow = datetime(2014,11,14,15,26,10); %end time for Wenglor profile demo

%folders for loading data, plotting, and functions
folder_LoadData_BSNE = '../../AnalysisData/BSNE/'; %folder for retrieving BSNE data
folder_LoadData_Wenglor = '../../AnalysisData/Windowing/'; %folder for retrieving Wenglor data
folder_LoadData_Subwindow = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Plots = '../../PlotOutput/Methods/'; %folder for plot output
folder_Functions = '../Functions/'; %folder with functions

%paths for loading data
LoadData_BSNE_Path = strcat(folder_LoadData_BSNE,'FluxBSNE'); %path for retrieving BSNE data
LoadData_Wenglor_Path = strcat(folder_LoadData_Wenglor,'DataWindowCalcs_30min_Restricted'); %path for retrieving Wenglor data
LoadData_Subwindow_Path = strcat(folder_LoadData_Subwindow,'DataFullSubwindowCalcs_30min_Restricted'); %path for 30 minute data

%plotting info
plot_fontsize = 13;
plot_markersize = 8;
plot_xsize = 5;
plot_ysize = 4;

%% LOAD DATA AND FUNCTIONS
%load data
load(LoadData_BSNE_Path); %load BSNE data
load(LoadData_Wenglor_Path); %load Wenglor data
load(LoadData_Subwindow_Path); %load subwindow data

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% EXTRACT RELEVANT DATA FOR PLOTS
%get indices for specific time intervals - BSNE
ind_Site = find(strcmp(Sites,Site_Demo)); %get index for demo site
ind_BSNE = find([FluxBSNE{ind_Site}.StartTime]==StartTime_Demo_BSNE); %get index of BSNE start time
% ind_Wenglor = find(StartTimes_all{ind_Site}==StartTime_Demo_Wenglor); %get index of BSNE start time

%get indices for specific time intervals - Wenglor
ind_T = find(T_subwindow == EndTime_Demo_Subwindow-StartTime_Demo_Subwindow);
ind_Subwindow = find(StartTime_Demo_Subwindow == StartTime_subwindow{ind_Site}{ind_T});

%extract height and flux information for BSNE profile plot
z_BSNE = FluxBSNE{ind_Site}(ind_BSNE).z.z; %get BSNE heights
sigma_z_BSNE = FluxBSNE{ind_Site}(ind_BSNE).z.sigma_z; %get BSNE height uncertainty
q_BSNE = FluxBSNE{ind_Site}(ind_BSNE).qz.qz; %get BSNE fluxes
sigma_q_BSNE = FluxBSNE{ind_Site}(ind_BSNE).qz.sigma; %get BSNE flux uncertainties
qpred_BSNE = FluxBSNE{ind_Site}(ind_BSNE).qz.qz_pred; %get predicted BSNE fluxes

%sort BSNE profile information by height
[~,ind_sort] = sort(z_BSNE);
z_BSNE = z_BSNE(ind_sort);
sigma_z_BSNE = sigma_z_BSNE(ind_sort);
q_BSNE = q_BSNE(ind_sort);
sigma_q_BSNE = sigma_q_BSNE(ind_sort);
qpred_BSNE = qpred_BSNE(ind_sort);

%extract fit information for BSNE profile plot
zq_BSNE = FluxBSNE{ind_Site}(ind_BSNE).z.zq; %get BSNE fitted saltation height
sigma_zq_BSNE = FluxBSNE{ind_Site}(ind_BSNE).z.sigma_zq; %get BSNE fitted saltation height
q0_BSNE = FluxBSNE{ind_Site}(ind_BSNE).qz.q0; %get BSNE fitted saltation scaling
sigma_q0_BSNE = FluxBSNE{ind_Site}(ind_BSNE).qz.sigma_q0; %get BSNE fitted saltation scaling
sigma2_q0zq_BSNE = FluxBSNE{ind_Site}(ind_BSNE).Q.sigma2_q0zq; %get BSNE fit parameter uncertainty

%extract height information for Wenglor calibration plots
z_Wenglor = zW_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get Wenglor heights
sigma_z_Wenglor = sigma_zW_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get Wenglor height uncertainties

% %extract height information for Wenglor calibration plots
% z_Wenglor = zW_all{ind_Site}{ind_Wenglor}; %get Wenglor heights
% sigma_z_Wenglor = sigma_zW_all{ind_Site}{ind_Wenglor}; %get Wenglor height uncertainties

%compute predicted Wenglor flux and its uncertainty for Wenglor calibration plot
qpred_calibration_Wenglor = q0_BSNE*exp(-z_Wenglor/zq_BSNE); %predicted Wenglor fluxes at BSNE heights 
sigma_qpred_calibration_Wenglor = qpred_calibration_Wenglor.*... %uncertainty in predicted Wenglor fluxes at BSNE heights
    sqrt((sigma_q0_BSNE/q0_BSNE)^2+...
    (sigma_zq_BSNE*z_Wenglor/zq_BSNE^2).^2+...
    2*sigma2_q0zq_BSNE*z_Wenglor/(q0_BSNE*zq_BSNE^2));

%extract calibration coefficients for Wenglor calibration plots
Cqn_Wenglor = Cqnbar_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get calibration coefficients
sigma_Cqn_Wenglor = sigma_Cqnbar_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get calibration coefficient uncertainties

% %extract calibration coefficients for Wenglor calibration plots
% Cqn_Wenglor = Cqnbar_all{ind_Site}{ind_Wenglor}; %get calibration coefficients
% sigma_Cqn_Wenglor = sigma_Cqnbar_all{ind_Site}{ind_Wenglor}; %get calibration coefficient uncertainties

%compute Wenglor counts for Wenglor calibration plot
n_calibration_Wenglor = qpred_calibration_Wenglor./Cqn_Wenglor; %get counts rates for Wenglor calibration

%extract counts rates for calibrated flux profile
n_Wenglor = nbar_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get Wenglor counts rate
sigma_n_Wenglor = sigma_nbar_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get calibration coefficient uncertainties

% %extract counts rates for calibrated flux profile
% n_Wenglor = nbar_all{ind_Site}{ind_Wenglor}; %get Wenglor counts rate
% sigma_n_Wenglor = sigma_nbar_all{ind_Site}{ind_Wenglor}; %get calibration coefficient uncertainties

%extract Wenglor flux information for calibrated flux profile
q_Wenglor = qbar_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get Wenglor fluxes
sigma_q_Wenglor = sigma_qbar_subwindow{ind_Site}{ind_T}{ind_Subwindow}; %get Wenglor flux uncertainties

% %extract Wenglor flux information for calibrated flux profile
% q_Wenglor = qbar_all{ind_Site}{ind_Wenglor}; %get Wenglor fluxes
% sigma_q_Wenglor = sigma_qbar_all{ind_Site}{ind_Wenglor}; %get Wenglor flux uncertainties

%extract fit information for Wenglor profile plot
zq_Wenglor = zq_subwindow{ind_Site}{ind_T}(ind_Subwindow); %get Wenglor fitted saltation height
Q_Wenglor_fit = Qfit_subwindow{ind_Site}{ind_T}(ind_Subwindow); %get Wenglor total flux - fitting
Q_Wenglor_sum = Qsum_subwindow{ind_Site}{ind_T}(ind_Subwindow); %get Wenglor total flux - summation

% %extract fit information for Wenglor profile plot
% zq_Wenglor = zq_all{ind_Site}(ind_Wenglor); %get Wenglor fitted saltation height
% Q_Wenglor = Q_all{ind_Site}(ind_Wenglor); %get Wenglor total flux

%get predicted Wenglor flux profile
q0_Wenglor = Q_Wenglor_fit./zq_Wenglor; %compute Wenglor fitted saltation scaling
% q0_Wenglor = Q_Wenglor./zq_Wenglor; %compute Wenglor fitted saltation scaling
qpred_Wenglor = q0_Wenglor*exp(-z_Wenglor/zq_Wenglor); %predicted Wenglor fluxes at Wenglor heights 
z_linearplot = linspace(0,max(max(z_Wenglor),max(z_BSNE)),100); %more z values for plot on linear scale
qpred_linearplot = q0_Wenglor*exp(-z_linearplot/zq_Wenglor); %predicted Wenglor fluxes at wide range of heights 

%compute delta_z weights for additive flux calculations
z1_Qsum = [0 sqrt(z_Wenglor(1:end-1).*z_Wenglor(2:end))];
z2_Qsum = [sqrt(z_Wenglor(1:end-1).*z_Wenglor(2:end)) Inf];
q0_Qsum = q_Wenglor./(exp(-z_Wenglor/zq_BSNE));
deltaQ = q0_Qsum.*zq_BSNE.*(exp(-z1_Qsum/zq_BSNE)-exp(-z2_Qsum/zq_BSNE));
deltaz = deltaQ./q_Wenglor; %equivalent deltaz for summation
%Q_Wenglor_sum = sum(deltaQ);

%% GET PLOTTING LIMITS

%get ranges of z values for plotting
z_lims = [0, ceil(max(max(z_BSNE),max(z_Wenglor))*10)/10]; %z limits

%get ranges of q values for plotting
q_min = min(min(min(qpred_BSNE),min(q_BSNE)),min(min(qpred_calibration_Wenglor),min(q_Wenglor))); %mininum q across plots
q_max = max(max(max(qpred_BSNE),max(q_BSNE)),max(max(qpred_calibration_Wenglor),max(q_Wenglor))); %maximum q across plots
q_min_scale = 10^floor(log10(q_min)); %scale (decade) for min q
q_max_scale = 10^floor(log10(q_max)); %scale (decade) for max q
q_lims = [floor(q_min/q_min_scale)*q_min_scale,ceil(q_max/q_max_scale)*q_max_scale]; %q limits

%get ranges of n values for plotting
n_min = min(min(n_Wenglor),min(n_calibration_Wenglor)); %mininum n across plots
n_max = max(max(n_Wenglor),max(n_calibration_Wenglor)); %maximum n across plots
n_min_scale = 10^floor(log10(n_min)); %scale (decade) for min n
n_max_scale = 10^floor(log10(n_max)); %scale (decade) for max n
n_lims = [floor(n_min/n_min_scale)*n_min_scale,ceil(n_max/n_max_scale)*n_max_scale]; %n limits

%get ranges of Cqn values for plotting
Cqn_min = min(Cqn_Wenglor); %mininum Cqn
Cqn_max = max(Cqn_Wenglor); %maximum Cqn
Cqn_min_scale = 10^floor(log10(Cqn_min)); %scale (decade) for min Cqn
Cqn_max_scale = 10^floor(log10(Cqn_max)); %scale (decade) for max Cqn
Cqn_lims = [floor(Cqn_min/Cqn_min_scale)*Cqn_min_scale,ceil(Cqn_max/Cqn_max_scale)*Cqn_max_scale]; %Cqn limits

%% PLOT BSNE profile
figure(1); clf;
subplot('Position',[0.07 0.17 0.35 0.76]); hold on;

%BSNE fluxes with error bars
h1 = plot(q_BSNE,z_BSNE,'ko','MarkerSize',plot_markersize);
for i = 1:length(z_BSNE)
    plot(q_BSNE(i)*[1 1],z_BSNE(i)+sigma_z_BSNE(i)*[-1 1],'k');
    plot(q_BSNE(i)+sigma_q_BSNE(i)*[-1 1],z_BSNE(i)*[1 1],'k');
end

%BSNE profile fit
h2 = plot(qpred_BSNE,z_BSNE,'k');

%predicted Wenglor fluxes with error bars
h3 = plot(qpred_calibration_Wenglor,z_Wenglor,'bs','MarkerSize',plot_markersize);
for i = 1:length(z_Wenglor)
    plot(qpred_calibration_Wenglor(i)*[1 1],z_Wenglor(i)+sigma_z_Wenglor(i)*[-1 1],'b');
    plot(qpred_calibration_Wenglor(i)+sigma_qpred_calibration_Wenglor(i)*[-1 1],z_Wenglor(i)*[1 1],'b');
end

%annotate plot
ylabel('height above bed surface, $$z$$','Interpreter','Latex');
xlabel('height-specific saltation flux, $$q$$ (gm$$^{-2}$$s$$^{-1}$$)','Interpreter','Latex');
h_legend = legend([h1 h2 h3],'LF trap obs, q_{LF,i}','LF profile fit, q_{LF,exp}','HF sensor pred, q_{pred,HF,i}','Location','SouthWest');
set(h_legend,'FontSize',plot_fontsize-3);
text(q_lims(2)-0.3*(q_lims(2)/q_lims(1)),z_lims(1)+0.95*diff(z_lims),'(a)','FontSize',plot_fontsize)

%format plot
title([Site_Demo,', ',datestr(StartTime_Demo_BSNE,'yyyy-mm-dd HH:MM'),'-',datestr(EndTime_Demo_BSNE,'HH:MM')])
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
set(gcf,'PaperUnits','inches','PaperSize',[plot_xsize plot_ysize],'PaperPosition',[0 0 plot_xsize plot_ysize],'PaperPositionMode','Manual');
xlim(q_lims);
ylim(z_lims);

% %print plot
% print([folder_Plots,'CalibrationDemo_BSNEFlux.png'],'-dpng'); %for draft
% %print([folder_Plots,'CalibrationDemo_BSNEFlux.tif'],'-dtiff','-r600'); %for publication

%% PLOT Wenglor counts profile - for calibration
subplot('Position',[0.47 0.17 0.23 0.76]); hold on;

%predicted Wenglor fluxes with error bars
plot(n_calibration_Wenglor,z_Wenglor,'bs','MarkerSize',plot_markersize);

%annotate plot
xlabel('HF counts rate, $$n_{cal,HF,i}$$ (s$$^{-1}$$)','Interpreter','Latex');
text(n_lims(2)-0.9*(n_lims(2)/n_lims(1)),z_lims(1)+0.95*diff(z_lims),'(b)','FontSize',plot_fontsize)

%format plot
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
%xlim(n_lims);
xlim([1 n_lims(2)]);
ylim(z_lims);

%% PLOT Wenglor calibration coefficients - for calibration
subplot('Position',[0.75 0.17 0.2 0.76]); hold on;

%predicted Wenglor fluxes with error bars
h1 = plot(Cqn_Wenglor,z_Wenglor,'rd','MarkerSize',plot_markersize);
for i = 1:length(z_Wenglor)
    plot(Cqn_Wenglor(i)+sigma_Cqn_Wenglor(i)*[-1 1],z_Wenglor(i)*[1 1],'r');
end

%annotate plot
xlabel('HF calibration, $$C_{qn,i}$$ (gm$$^{-2}$$)','Interpreter','Latex');
text(Cqn_lims(2)-0.9*(Cqn_lims(2)/Cqn_lims(1)),z_lims(1)+0.95*diff(z_lims),'(c)','FontSize',plot_fontsize)

%format plot
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
set(gcf,'PaperUnits','inches','PaperSize',[plot_xsize plot_ysize],'PaperPosition',[0 0 plot_xsize plot_ysize],'PaperPositionMode','Manual');
xlim(Cqn_lims);
xlim([5 Cqn_lims(2)]);
ylim(z_lims);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 4]);
print([folder_Plots,'CalibrationDemo.png'],'-dpng');


%% PLOT Wenglor calibration coefficients -- for flux prediction
figure(2); clf;
subplot('Position',[0.07 0.17 0.23 0.76]); hold on;

%predicted Wenglor fluxes with error bars
h1 = plot(Cqn_Wenglor,z_Wenglor,'rd','MarkerSize',plot_markersize);
for i = 1:length(z_Wenglor)
    plot(Cqn_Wenglor(i)+sigma_Cqn_Wenglor(i)*[-1 1],z_Wenglor(i)*[1 1],'r','LineWidth',1);
end

%annotate plot
xlabel('HF calibration, $$C_{qn,i}$$ (gm$$^{-2}$$)','Interpreter','Latex');
ylabel('height above bed surface, $$z$$','Interpreter','Latex');
% title([Site_Demo,', ',datestr(StartTime_Demo_Wenglor,'yyyy-mm-dd HH:MM'),'-',datestr(EndTime_Demo_Wenglor,'HH:MM')])
text(Cqn_lims(2)-0.7*(Cqn_lims(2)/Cqn_lims(1)),0.51,'(a)')
%text(Cqn_lims(2)-0.7*(Cqn_lims(2)/Cqn_lims(1)),z_lims(1)+0.95*diff(z_lims),'(a)')

%format plot
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
xlim(Cqn_lims);
ylim([0 0.55]);


%% PLOT Wenglor counts profile - for flux prediction
subplot('Position',[0.35 0.17 0.25 0.76]); hold on;

%predicted Wenglor fluxes with error bars
plot(n_Wenglor,z_Wenglor,'bs','MarkerSize',plot_markersize*2);
for i = 1:length(z_Wenglor)
    plot(n_Wenglor(i)+sigma_n_Wenglor(i)*[-1 1],z_Wenglor(i)*[1 1],'b','LineWidth',1);
end

%annotate plot
title([Site_Demo,', ',datestr(StartTime_Demo_Subwindow,'yyyy-mm-dd HH:MM:SS'),'-',datestr(EndTime_Demo_Subwindow,'HH:MM:SS')])
xlabel('HF counts rate, $$n_{i}$$ (s$$^{-1}$$)','Interpreter','Latex');
%text(n_lims(2)-0.7*(n_lims(2)/n_lims(1)),z_lims(1)+0.95*diff(z_lims),'(b)')
text(n_lims(2)-0.7*(n_lims(2)/n_lims(1)),0.51,'(b)')

%format plot
set(gca,'xscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
xlim(n_lims);
ylim([0 0.55]);


%% PLOT calibrated Wenglor fluxes
subplot('Position',[0.65 0.17 0.3 0.76]); hold on;

%Wenglor weighted sum
for i = 1:length(z_Wenglor)
    if i == 1
        h3 = plot([0 q_Wenglor(i) q_Wenglor(i) 0],[0 0 z2_Qsum(i) z2_Qsum(i)],'r','LineWidth',2);
        %h3 = area([q_Wenglor(i) q_Wenglor(i)],[z2_Qsum(i)-deltaz(i) z2_Qsum(i)],'FaceColor','r');
    elseif i == length(z_Wenglor)
        plot([0 q_Wenglor(i) q_Wenglor(i) 0],[z1_Qsum(i) z1_Qsum(i) z1_Qsum(i)+deltaz(i) z1_Qsum(i)+deltaz(i)],'r','LineWidth',2);
        %area([q_Wenglor(i) q_Wenglor(i)],[z1_Qsum(i) z1_Qsum(i)+deltaz(i)],'FaceColor','r');
    else
        plot([0 q_Wenglor(i) q_Wenglor(i) 0],[z1_Qsum(i) z1_Qsum(i) z2_Qsum(i) z2_Qsum(i)],'r','LineWidth',2);
        %area([q_Wenglor(i) q_Wenglor(i)],[z1_Qsum(i) z2_Qsum(i)],'FaceColor','r');
    end
end

%predicted Wenglor fluxes with error bars
h1 = plot(q_Wenglor,z_Wenglor,'bs','MarkerSize',plot_markersize*2,'LineWidth',1);
for i = 1:length(z_Wenglor)
    plot(q_Wenglor(i)+sigma_q_Wenglor(i)*[-1 1],z_Wenglor(i)*[1 1],'b','LineWidth',1);
end

%Wenglor profile fit
%h2 = plot(z_Wenglor,qpred_Wenglor,'b');
h2 = plot(qpred_linearplot,z_linearplot,'k','LineWidth',2);

%annotate plot
xlabel('Calibrated HF flux, $$q_{i}$$ (gm$$^{-2}$$s$$^{-1}$$)','Interpreter','Latex');
%text(q_lims(1)+0.9*diff(q_lims),z_lims(1)+0.95*diff(z_lims),'(c)')
text(q_lims(1)+0.9*diff(q_lims),0.51,'(c)')
legend([h1 h2 h3],'HF sensor cal. flux','HF profile fit','HF profile sum','Location','North');

%plotting info
plot_fontsize = 14;
plot_xsize = 5;
plot_ysize = 4;

%format plot
%set(gca,'yscale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','FontSize',plot_fontsize);
xlim(q_lims);
ylim([0 0.55]);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 4]);
print([folder_Plots,'CalibrationPrediction.png'],'-dpng'); 