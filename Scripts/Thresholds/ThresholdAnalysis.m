%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% Site-specific parameter values
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
z0_Site = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m)
tauit_intercept = [0.137, 0.106, 0.086];
sigma_tauit_intercept = [0.015, 0.020, 0.008];

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,1); %sampling interval for analysis

%% diurnal info
diurnalrange_starthour = [12; 14; 16];
diurnalrange_endhour = [14; 16; 18];

%% date info
daterange_startdate = [datetime(2015,5,15);datetime(2015,5,23);datetime(2015,6,1)];
daterange_enddate = [datetime(2015,5,19);datetime(2015,5,28);datetime(2015,6,4)];

%% plotting info
PlotFont = 14;
PlotMarkers = {'s','d','o','<','>'};
PlotColors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};

%% load data 
folder_LoadData = '../../AnalysisData/Thresholds/'; %folder for retrieving processed data
LoadData_Path = strcat(folder_LoadData,'ThresholdAnalysisWindows'); %path for loading time window data
load(LoadData_Path);

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% specify where to save data and plots
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
SaveData_Path = strcat(folder_SaveData,'ThresholdAnalysisData'); %path for saving output data
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMARY ANALYSIS BY SITE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize list of all values by site
fQ_bin_avg_all = cell(N_Sites,1); %transport activity - average
fQ_bin_SE_all = cell(N_Sites,1); %transport activity - standard error
zs_log10_bin_avg_all = cell(N_Sites,1); %roughness height - average
zs_log10_bin_SE_all = cell(N_Sites,1); %roughness height - standard error
uth_bin_avg_all = cell(N_Sites,1); %TFEM threshold wind - average
uth_bin_SE_all = cell(N_Sites,1); %TFEM threshold wind - standard error
ustth_bin_avg_all = cell(N_Sites,1); %TFEM threshold shear velocity - average
ustth_bin_SE_all = cell(N_Sites,1); %TFEM threshold shear velocity - standard error
tauth_bin_avg_all = cell(N_Sites,1); %TFEM threshold stress - average
tauth_bin_SE_all = cell(N_Sites,1); %TFEM threshold stress - standard error
tauft_all = zeros(N_Sites,1); %inferred fluid threshold stress - best fit
sigma_tauft_all = zeros(N_Sites,1); %inferred fluid threshold stress - uncertainty
tauit_all = zeros(N_Sites,1); %inferred impact threshold stress - best fit
sigma_tauit_all = zeros(N_Sites,1); %inferred impact threshold stress - uncertainty
ustitftratio_all = zeros(N_Sites,1); %ratio of impact and fluid threshold shear velocities - best fit
sigma_ustitftratio_all = zeros(N_Sites,1); %ratio of impact and fluid threshold shear velocities - uncertainty
Chi2nu_all = zeros(N_Sites,1); %normalized Chi2 for threshold versus fQ fit

%% perform main analysis by site
for i = 1:length(Sites)
    
    %% get parameter values
    rho = rho_Site(i); %get air density
    z0 = z0_Site(i); %get roughness length
    
    %% apply threshold binning analysis function
    [fQ_bin_avg,fQ_bin_SE,...
        zs_log10_bin_avg,zs_log10_bin_SE,...
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio,...
        Chi2nu] = ...
    ThresholdBinning(rho,z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Sites{i},...
        fQ_all, uth_all, zU_all, zs_all,...
        starttime_all,timeofday_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);
        
    %% assign values by site
    fQ_bin_avg_all{i} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_all{i} = fQ_bin_SE; %transport activity - standard error
    zs_log10_bin_avg_all{i} = zs_log10_bin_avg; %roughness height - average
    zs_log10_bin_SE_all{i} = zs_log10_bin_SE; %roughness height - standard error
    uth_bin_avg_all{i} = uth_bin_avg; %TFEM threshold wind - average
    uth_bin_SE_all{i} = uth_bin_SE; %TFEM threshold wind - standard error
    ustth_bin_avg_all{i} = ustth_bin_avg; %TFEM threshold shear velocity - average
    ustth_bin_SE_all{i} = ustth_bin_SE; %TFEM threshold shear velocity - standard error
    tauth_bin_avg_all{i} = tauth_bin_avg; %TFEM threshold stress - average
    tauth_bin_SE_all{i} = tauth_bin_SE; %TFEM threshold stress - standard error
    tauft_all(i) = tauft; %inferred fluid threshold stress - best fit
    sigma_tauft_all(i) = sigma_tauft; %inferred fluid threshold stress - uncertainty
    tauit_all(i) = tauit; %inferred impact threshold stress - best fit
    sigma_tauit_all(i) = sigma_tauit; %inferred impact threshold stress - uncertainty
    ustitftratio_all(i) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_all(i) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
    Chi2nu_all(i) = Chi2nu; %normalized Chi2 for threshold versus fQ fit
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform sensitivity analyses for Oceano %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Site_sensitivityanalysis = 'Oceano'; %get site name
ind_Site = find(strcmp(Sites,Site_sensitivityanalysis)); %get index for site
rho = rho_Site(ind_Site); %get air density
z0 = z0_Site(ind_Site); %get roughness length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensitivity analysis for roughness height, get comparison for fQ = 0 at Oceano
m = find(Deltat_all==Deltat_analysis);
s = find(deltat_all==deltat_analysis);
ind_fQ0 = find(fQ_all{ind_Site}{m,s}==0);
z0_fQ0 = zs_all{ind_Site}{m,s}(ind_fQ0);
z0_log10_avg = mean(log10(z0_fQ0));
z0_log10_SE = std(log10(z0_fQ0))/sqrt(length(z0_fQ0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensitivity analysis for measurement interval
N_Deltat = length(Deltat_all);

%% initialize list of all values by Deltat
fQ_bin_avg_Deltat = cell(N_Deltat,1); %transport activity - average
fQ_bin_SE_Deltat = cell(N_Deltat,1); %transport activity - standard error
zs_log10_bin_avg_Deltat = cell(N_Deltat,1); %roughness height - average
zs_log10_bin_SE_Deltat = cell(N_Deltat,1); %roughness height - standard error
uth_bin_avg_Deltat = cell(N_Deltat,1); %TFEM threshold wind - average
uth_bin_SE_Deltat = cell(N_Deltat,1); %TFEM threshold wind - standard error
ustth_bin_avg_Deltat = cell(N_Deltat,1); %TFEM threshold shear velocity - average
ustth_bin_SE_Deltat = cell(N_Deltat,1); %TFEM threshold shear velocity - standard error
tauth_bin_avg_Deltat = cell(N_Deltat,1); %TFEM threshold stress - average
tauth_bin_SE_Deltat = cell(N_Deltat,1); %TFEM threshold stress - standard error
tauft_Deltat = zeros(N_Deltat,1); %inferred fluid threshold stress - best fit
sigma_tauft_Deltat = zeros(N_Deltat,1); %inferred fluid threshold stress - uncertainty
tauit_Deltat = zeros(N_Deltat,1); %inferred impact threshold stress - best fit
sigma_tauit_Deltat = zeros(N_Deltat,1); %inferred impact threshold stress - uncertainty
ustitftratio_Deltat = zeros(N_Deltat,1); %ratio of impact and fluid threshold shear velocities - best fit
sigma_ustitftratio_Deltat = zeros(N_Deltat,1); %ratio of impact and fluid threshold shear velocities - uncertainty

%% perform analysis by measurement interval
for m = 1:N_Deltat
  
    %apply threshold binning analysis function
    [fQ_bin_avg,fQ_bin_SE,...
        zs_log10_bin_avg,zs_log10_bin_SE,...
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,...
        Deltat_all,deltat_all,...
        Deltat_all(m), deltat_analysis,...
        Sites,Site_sensitivityanalysis,...
        fQ_all, uth_all, zU_all, zs_all,...
        starttime_all,timeofday_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);

    %% assign values by Deltat
    fQ_bin_avg_Deltat{m} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_Deltat{m} = fQ_bin_SE; %transport activity - standard error
    zs_log10_bin_avg_Deltat{m} = zs_log10_bin_avg; %roughness height - average
    zs_log10_bin_SE_Deltat{m} = zs_log10_bin_SE; %roughness height - standard error
    uth_bin_avg_Deltat{m} = uth_bin_avg; %TFEM threshold wind - average
    uth_bin_SE_Deltat{m} = uth_bin_SE; %TFEM threshold wind - standard error
    ustth_bin_avg_Deltat{m} = ustth_bin_avg; %TFEM threshold shear velocity - average
    ustth_bin_SE_Deltat{m} = ustth_bin_SE; %TFEM threshold shear velocity - standard error
    tauth_bin_avg_Deltat{m} = tauth_bin_avg; %TFEM threshold stress - average
    tauth_bin_SE_Deltat{m} = tauth_bin_SE; %TFEM threshold stress - standard error
    tauft_Deltat(m) = tauft; %inferred fluid threshold stress - best fit
    sigma_tauft_Deltat(m) = sigma_tauft; %inferred fluid threshold stress - uncertainty
    tauit_Deltat(m) = tauit; %inferred impact threshold stress - best fit
    sigma_tauit_Deltat(m) = sigma_tauit; %inferred impact threshold stress - uncertainty
    ustitftratio_Deltat(m) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_Deltat(m) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensitivity analysis by sampling interval
N_deltat = length(deltat_all);

%% initialize list of all values by sampling interval
fQ_bin_avg_deltat = cell(N_deltat,1); %transport activity - average
fQ_bin_SE_deltat = cell(N_deltat,1); %transport activity - standard error
zs_log10_bin_avg_deltat = cell(N_deltat,1); %roughness height - average
zs_log10_bin_SE_deltat = cell(N_deltat,1); %roughness height - standard error
uth_bin_avg_deltat = cell(N_deltat,1); %TFEM threshold wind - average
uth_bin_SE_deltat = cell(N_deltat,1); %TFEM threshold wind - standard error
ustth_bin_avg_deltat = cell(N_deltat,1); %TFEM threshold shear velocity - average
ustth_bin_SE_deltat = cell(N_deltat,1); %TFEM threshold shear velocity - standard error
tauth_bin_avg_deltat = cell(N_deltat,1); %TFEM threshold stress - average
tauth_bin_SE_deltat = cell(N_deltat,1); %TFEM threshold stress - standard error
tauft_deltat = zeros(N_deltat,1); %inferred fluid threshold stress - best fit
sigma_tauft_deltat = zeros(N_deltat,1); %inferred fluid threshold stress - uncertainty
tauit_deltat = zeros(N_deltat,1); %inferred impact threshold stress - best fit
sigma_tauit_deltat = zeros(N_deltat,1); %inferred impact threshold stress - uncertainty
ustitftratio_deltat = zeros(N_deltat,1); %ratio of impact and fluid threshold shear velocities - best fit
sigma_ustitftratio_deltat = zeros(N_deltat,1); %ratio of impact and fluid threshold shear velocities - uncertainty

%% perform analysis by measurement interval
for s = 1:N_deltat
  
    %apply threshold binning analysis function
    [fQ_bin_avg,fQ_bin_SE,...
        zs_log10_bin_avg,zs_log10_bin_SE,...
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_all(s),...
        Sites,Site_sensitivityanalysis,...
        fQ_all, uth_all, zU_all, zs_all,...
        starttime_all,timeofday_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);

    %% assign values by sampling interval
    fQ_bin_avg_deltat{s} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_deltat{s} = fQ_bin_SE; %transport activity - standard error
    zs_log10_bin_avg_deltat{s} = zs_log10_bin_avg; %roughness height - average
    zs_log10_bin_SE_deltat{s} = zs_log10_bin_SE; %roughness height - standard error
    uth_bin_avg_deltat{s} = uth_bin_avg; %TFEM threshold wind - average
    uth_bin_SE_deltat{s} = uth_bin_SE; %TFEM threshold wind - standard error
    ustth_bin_avg_deltat{s} = ustth_bin_avg; %TFEM threshold shear velocity - average
    ustth_bin_SE_deltat{s} = ustth_bin_SE; %TFEM threshold shear velocity - standard error
    tauth_bin_avg_deltat{s} = tauth_bin_avg; %TFEM threshold stress - average
    tauth_bin_SE_deltat{s} = tauth_bin_SE; %TFEM threshold stress - standard error
    tauft_deltat(s) = tauft; %inferred fluid threshold stress - best fit
    sigma_tauft_deltat(s) = sigma_tauft; %inferred fluid threshold stress - uncertainty
    tauit_deltat(s) = tauit; %inferred impact threshold stress - best fit
    sigma_tauit_deltat(s) = sigma_tauit; %inferred impact threshold stress - uncertainty
    ustitftratio_deltat(s) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_deltat(s) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensitivity analysis by diurnal range
N_diurnalrange = length(diurnalrange_starthour);

%% initialize list of all values by sampling interval
fQ_bin_avg_diurnalrange = cell(N_diurnalrange,1); %transport activity - average
fQ_bin_SE_diurnalrange = cell(N_diurnalrange,1); %transport activity - standard error
zs_log10_bin_avg_diurnalrange = cell(N_diurnalrange,1); %roughness height - average
zs_log10_bin_SE_diurnalrange = cell(N_diurnalrange,1); %roughness height - standard error
uth_bin_avg_diurnalrange = cell(N_diurnalrange,1); %TFEM threshold wind - average
uth_bin_SE_diurnalrange = cell(N_diurnalrange,1); %TFEM threshold wind - standard error
ustth_bin_avg_diurnalrange = cell(N_diurnalrange,1); %TFEM threshold shear velocity - average
ustth_bin_SE_diurnalrange = cell(N_diurnalrange,1); %TFEM threshold shear velocity - standard error
tauth_bin_avg_diurnalrange = cell(N_diurnalrange,1); %TFEM threshold stress - average
tauth_bin_SE_diurnalrange = cell(N_diurnalrange,1); %TFEM threshold stress - standard error
tauft_diurnalrange = zeros(N_diurnalrange,1); %inferred fluid threshold stress - best fit
sigma_tauft_diurnalrange = zeros(N_diurnalrange,1); %inferred fluid threshold stress - uncertainty
tauit_diurnalrange = zeros(N_diurnalrange,1); %inferred impact threshold stress - best fit
sigma_tauit_diurnalrange = zeros(N_diurnalrange,1); %inferred impact threshold stress - uncertainty
ustitftratio_diurnalrange = zeros(N_diurnalrange,1); %ratio of impact and fluid threshold shear velocities - best fit
sigma_ustitftratio_diurnalrange = zeros(N_diurnalrange,1); %ratio of impact and fluid threshold shear velocities - uncertainty

%% perform analysis by diurnal range
for d = 1:N_diurnalrange
  
    %apply threshold binning analysis function
    [fQ_bin_avg,fQ_bin_SE,...
        zs_log10_bin_avg,zs_log10_bin_SE,...
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Site_sensitivityanalysis,...
        fQ_all, uth_all, zU_all, zs_all,...
        starttime_all,timeofday_all,...
        diurnalrange_starthour(d),diurnalrange_endhour(d));%,...
        %startdate_analysis,enddate_analysis);

    %% assign values by sampling interval
    fQ_bin_avg_diurnalrange{d} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_diurnalrange{d} = fQ_bin_SE; %transport activity - standard error
    zs_log10_bin_avg_diurnalrange{d} = zs_log10_bin_avg; %roughness height - average
    zs_log10_bin_SE_diurnalrange{d} = zs_log10_bin_SE; %roughness height - standard error
    uth_bin_avg_diurnalrange{d} = uth_bin_avg; %TFEM threshold wind - average
    uth_bin_SE_diurnalrange{d} = uth_bin_SE; %TFEM threshold wind - standard error
    ustth_bin_avg_diurnalrange{d} = ustth_bin_avg; %TFEM threshold shear velocity - average
    ustth_bin_SE_diurnalrange{d} = ustth_bin_SE; %TFEM threshold shear velocity - standard error
    tauth_bin_avg_diurnalrange{d} = tauth_bin_avg; %TFEM threshold stress - average
    tauth_bin_SE_diurnalrange{d} = tauth_bin_SE; %TFEM threshold stress - standard error
    tauft_diurnalrange(d) = tauft; %inferred fluid threshold stress - best fit
    sigma_tauft_diurnalrange(d) = sigma_tauft; %inferred fluid threshold stress - uncertainty
    tauit_diurnalrange(d) = tauit; %inferred impact threshold stress - best fit
    sigma_tauit_diurnalrange(d) = sigma_tauit; %inferred impact threshold stress - uncertainty
    ustitftratio_diurnalrange(d) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_diurnalrange(d) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensitivity analysis by date range
N_daterange = length(daterange_startdate);

%% initialize list of all values by sampling interval
fQ_bin_avg_daterange = cell(N_daterange,1); %transport activity - average
fQ_bin_SE_daterange = cell(N_daterange,1); %transport activity - standard error
zs_log10_bin_avg_daterange = cell(N_daterange,1); %roughness height - average
zs_log10_bin_SE_daterange = cell(N_daterange,1); %roughness height - standard error
uth_bin_avg_daterange = cell(N_daterange,1); %TFEM threshold wind - average
uth_bin_SE_daterange = cell(N_daterange,1); %TFEM threshold wind - standard error
ustth_bin_avg_daterange = cell(N_daterange,1); %TFEM threshold shear velocity - average
ustth_bin_SE_daterange = cell(N_daterange,1); %TFEM threshold shear velocity - standard error
tauth_bin_avg_daterange = cell(N_daterange,1); %TFEM threshold stress - average
tauth_bin_SE_daterange = cell(N_daterange,1); %TFEM threshold stress - standard error
tauft_daterange = zeros(N_daterange,1); %inferred fluid threshold stress - best fit
sigma_tauft_daterange = zeros(N_daterange,1); %inferred fluid threshold stress - uncertainty
tauit_daterange = zeros(N_daterange,1); %inferred impact threshold stress - best fit
sigma_tauit_daterange = zeros(N_daterange,1); %inferred impact threshold stress - uncertainty
ustitftratio_daterange = zeros(N_daterange,1); %ratio of impact and fluid threshold shear velocities - best fit
sigma_ustitftratio_daterange = zeros(N_daterange,1); %ratio of impact and fluid threshold shear velocities - uncertainty

%% perform analysis by date range
for d = 1:N_daterange
  
    %apply threshold binning analysis function
    [fQ_bin_avg,fQ_bin_SE,...
        zs_log10_bin_avg,zs_log10_bin_SE,...
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Site_sensitivityanalysis,...
        fQ_all, uth_all, zU_all, zs_all,...
        starttime_all,timeofday_all,...
        0,24,...
        daterange_startdate(d),daterange_enddate(d));

    %% assign values by sampling interval
    fQ_bin_avg_daterange{d} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_daterange{d} = fQ_bin_SE; %transport activity - standard error
    zs_log10_bin_avg_daterange{d} = zs_log10_bin_avg; %roughness height - average
    zs_log10_bin_SE_daterange{d} = zs_log10_bin_SE; %roughness height - standard error
    uth_bin_avg_daterange{d} = uth_bin_avg; %TFEM threshold wind - average
    uth_bin_SE_daterange{d} = uth_bin_SE; %TFEM threshold wind - standard error
    ustth_bin_avg_daterange{d} = ustth_bin_avg; %TFEM threshold shear velocity - average
    ustth_bin_SE_daterange{d} = ustth_bin_SE; %TFEM threshold shear velocity - standard error
    tauth_bin_avg_daterange{d} = tauth_bin_avg; %TFEM threshold stress - average
    tauth_bin_SE_daterange{d} = tauth_bin_SE; %TFEM threshold stress - standard error
    tauft_daterange(d) = tauft; %inferred fluid threshold stress - best fit
    sigma_tauft_daterange(d) = sigma_tauft; %inferred fluid threshold stress - uncertainty
    tauit_daterange(d) = tauit; %inferred impact threshold stress - best fit
    sigma_tauit_daterange(d) = sigma_tauit; %inferred impact threshold stress - uncertainty
    ustitftratio_daterange(d) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_daterange(d) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
end

%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%
save(SaveData_Path,...
    'fQ_bin_avg_*','fQ_bin_SE_*',...
    'tauft_*','sigma_tauft_*',...
    'tauit_*','sigma_tauit_*',...
    'ustitftratio_*','sigma_ustitftratio_*',...
    'Chi2nu_*');

%%%%%%%%%%%%
% PLOTTING %
%%%%%%%%%%%%

%% plot tau_th versus fQ for sites
figure(1); clf; hold on; %initialize plot

%plot average values
for i = 1:N_Sites
    plot(fQ_bin_avg_all{i},tauth_bin_avg_all{i},PlotMarkers{i},'Color',PlotColors{i});
end

%plot SE values
for i = 1:N_Sites
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],tauth_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{i});
        plot(fQ_bin_avg_all{i}(k)*[1 1],tauth_bin_avg_all{i}(k)+tauth_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{i});
    end
end

%plot fit lines
for i = 1:N_Sites
    plot([0 1],[tauft_all(i) tauit_all(i)],'Color',PlotColors{i});
end

%plot fit intercept values
for i = 1:N_Sites
    plot([0 0],tauft_all(i)+sigma_tauft_all(i)*[-1 1],'Color',PlotColors{i},'LineWidth',2);
    plot([0.999 0.999],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'Color',PlotColors{i},'LineWidth',2);
end

%annotate plot
legend(SiteNames,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity.png'],'-dpng');

%% plot zs versus fQ
figure(11); clf; hold on; %initialize plot

%plot average values
for i = 1:N_Sites
    plot(fQ_bin_avg_all{i},10.^zs_log10_bin_avg_all{i},PlotMarkers{i},'Color',PlotColors{i});
end

%plot SE values
for i = 1:N_Sites
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],10.^(zs_log10_bin_avg_all{i}(k)*[1 1]),'Color',PlotColors{i});
        plot(fQ_bin_avg_all{i}(k)*[1 1],10.^(zs_log10_bin_avg_all{i}(k)+zs_log10_bin_SE_all{i}(k)*[-1 1]),'Color',PlotColors{i});
    end
end

% %plot z0 for Oceano
% plot(0,10.^z0_log10_avg,PlotMarkers{ind_Site},'Color',PlotColors{ind_Site});
% plot([0 0],10.^(z0_log10_avg+z0_log10_SE*[-1 1]),'Color',PlotColors{ind_Site});

%annotate plot
legend(SiteNames,'Location','SouthWest');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective roughness length, $$z_{s}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'YScale','log');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'roughness_activity.png'],'-dpng');


%% plot tau_th versus fQ for measurement interval sensitivity analysis
figure(12); clf; hold on; %initialize plot

%plot average values
for m = 1:N_Deltat
    plot(fQ_bin_avg_Deltat{m},tauth_bin_avg_Deltat{m},PlotMarkers{m},'Color',PlotColors{m});
end

%plot SE values
for m = 1:N_Deltat
    N_fQ_bins = length(fQ_bin_avg_Deltat{m});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_Deltat{m}(k)+fQ_bin_SE_Deltat{m}(k)*[-1 1],tauth_bin_avg_Deltat{m}(k)*[1 1],'Color',PlotColors{m});
        plot(fQ_bin_avg_Deltat{m}(k)*[1 1],tauth_bin_avg_Deltat{m}(k)+tauth_bin_SE_Deltat{m}(k)*[-1 1],'Color',PlotColors{m});
    end
end

%plot fit lines
for m = 1:N_Deltat
    plot([0 1],[tauft_Deltat(m) tauit_Deltat(m)],'Color',PlotColors{m});
end

%plot fit intercept values
for m = 1:N_Deltat
    plot([0 0],tauft_Deltat(m)+sigma_tauft_Deltat(m)*[-1 1],'Color',PlotColors{m},'LineWidth',2);
    plot([0.999 0.999],tauit_Deltat(m)+sigma_tauit_Deltat(m)*[-1 1],'Color',PlotColors{m},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_Deltat,1);
for m = 1:N_Deltat
    legend_items{m} = ['\Delta t = ',num2str(minutes(Deltat_all(m))),' min.'];
end
legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity_measurementinterval.png'],'-dpng');

%% plot tau_th versus fQ for sampling interval sensitivity analysis
figure(13); clf; hold on; %initialize plot

%plot average values
for s = 1:N_deltat
    plot(fQ_bin_avg_deltat{s},tauth_bin_avg_deltat{s},PlotMarkers{s},'Color',PlotColors{s});
end

%plot SE values
for s = 1:N_deltat
    N_fQ_bins = length(fQ_bin_avg_deltat{s});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_deltat{s}(k)+fQ_bin_SE_deltat{s}(k)*[-1 1],tauth_bin_avg_deltat{s}(k)*[1 1],'Color',PlotColors{s});
        plot(fQ_bin_avg_deltat{s}(k)*[1 1],tauth_bin_avg_deltat{s}(k)+tauth_bin_SE_deltat{s}(k)*[-1 1],'Color',PlotColors{s});
    end
end

%plot fit lines
for s = 1:N_deltat
    plot([0 1],[tauft_deltat(s) tauit_deltat(s)],'Color',PlotColors{s});
end

%plot fit intercept values
for s = 1:N_deltat
    plot([0 0],tauft_deltat(s)+sigma_tauft_deltat(s)*[-1 1],'Color',PlotColors{s},'LineWidth',2);
    plot([0.999 0.999],tauit_deltat(s)+sigma_tauit_deltat(s)*[-1 1],'Color',PlotColors{s},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_deltat,1);
for s = 1:N_deltat
    legend_items{s} = ['\delta{t} = ',num2str(seconds(deltat_all(s)),2),' s'];
end
legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity_samplinginterval.png'],'-dpng');

%% plot tau_th versus fQ for diurnal period sensitivity analysis
figure(14); clf; hold on; %initialize plot

%plot average values
for d = 1:N_diurnalrange
    plot(fQ_bin_avg_diurnalrange{d},tauth_bin_avg_diurnalrange{d},PlotMarkers{d},'Color',PlotColors{d});
end

%plot SE values
for d = 1:N_diurnalrange
    N_fQ_bins = length(fQ_bin_avg_diurnalrange{d});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_diurnalrange{d}(k)+fQ_bin_SE_diurnalrange{d}(k)*[-1 1],tauth_bin_avg_diurnalrange{d}(k)*[1 1],'Color',PlotColors{d});
        plot(fQ_bin_avg_diurnalrange{d}(k)*[1 1],tauth_bin_avg_diurnalrange{d}(k)+tauth_bin_SE_diurnalrange{d}(k)*[-1 1],'Color',PlotColors{d});
    end
end

%plot fit lines
for d = 1:N_diurnalrange
    plot([0 1],[tauft_diurnalrange(d) tauit_diurnalrange(d)],'Color',PlotColors{d});
end

%plot fit intercept values
for d = 1:N_diurnalrange
    plot([0 0],tauft_diurnalrange(d)+sigma_tauft_diurnalrange(d)*[-1 1],'Color',PlotColors{d},'LineWidth',2);
    plot([0.999 0.999],tauit_diurnalrange(d)+sigma_tauit_diurnalrange(d)*[-1 1],'Color',PlotColors{d},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_diurnalrange,1);
for d = 1:N_diurnalrange
    legend_items{d} = [num2str(diurnalrange_starthour(d)),'-',num2str(diurnalrange_endhour(d)),'h'];
end

legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity_diurnalrange.png'],'-dpng');

%% plot tau_th versus fQ for date period sensitivity analysis
figure(15); clf; hold on; %initialize plot

%plot average values
for d = 1:N_daterange
    plot(fQ_bin_avg_daterange{d},tauth_bin_avg_daterange{d},PlotMarkers{d},'Color',PlotColors{d});
end

%plot SE values
for d = 1:N_daterange
    N_fQ_bins = length(fQ_bin_avg_daterange{d});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_daterange{d}(k)+fQ_bin_SE_daterange{d}(k)*[-1 1],tauth_bin_avg_daterange{d}(k)*[1 1],'Color',PlotColors{d});
        plot(fQ_bin_avg_daterange{d}(k)*[1 1],tauth_bin_avg_daterange{d}(k)+tauth_bin_SE_daterange{d}(k)*[-1 1],'Color',PlotColors{d});
    end
end

%plot fit lines
for d = 1:N_daterange
    plot([0 1],[tauft_daterange(d) tauit_daterange(d)],'Color',PlotColors{d});
end

%plot fit intercept values
for d = 1:N_daterange
    plot([0 0],tauft_daterange(d)+sigma_tauft_daterange(d)*[-1 1],'Color',PlotColors{d},'LineWidth',2);
    plot([0.999 0.999],tauit_daterange(d)+sigma_tauit_daterange(d)*[-1 1],'Color',PlotColors{d},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_daterange,1);
for d = 1:N_daterange
    legend_items{d} = [datestr(daterange_startdate(d),'mmm dd'),'-',datestr(daterange_enddate(d),'dd')];
end

legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity_daterange.png'],'-dpng');