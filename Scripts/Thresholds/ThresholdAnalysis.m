%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% Parameters
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,1); %sampling interval for analysis

%% filtering info
theta_max = 20; %maximum absolute wind angle for calcs
zL_max = 0.2; %maximum absolute stability value for calcs

%% diurnal info
diurnalrange_starthour = [12; 14; 16];
diurnalrange_endhour = [14; 16; 18];

%% date info
daterange_startdate = [datetime(2015,5,15);datetime(2015,5,23);datetime(2015,6,1)];
daterange_enddate = [datetime(2015,5,19);datetime(2015,5,28);datetime(2015,6,4)];

%% plotting info
PlotFont = 14;
PlotMarkers_Site = {'s','d','o'};
PlotColors_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
PlotMarkers_Sensitivity = {'s','d','o','<','>'};
PlotColors_Sensitivity = {[0.6473 0.7456 0.4188],[0.2116 0.1898 0.5777],[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};

%% folders for loading and saving data
folder_LoadFluxLawData = '../../AnalysisData/FluxLaw/'; %folder for loading 30 minute data
folder_LoadRoughnessData = '../../AnalysisData/Roughness/'; %folder for loading roughness data
folder_LoadSubwindowData = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
folder_LoadWindowData = '../../AnalysisData/Windowing/'; %folder for loading 30 minute data
folder_Functions = '../Functions/'; %folder with functions
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% paths for loading and saving data
LoadFluxLawData_Path = strcat(folder_LoadFluxLawData,'FluxLawCalcs_30min_Restricted'); %path for 30 minute data
LoadRoughnessData_Path = strcat(folder_LoadRoughnessData,'RoughnessCalcs_30min_Restricted'); %path for 30 minute data
LoadSubwindowData_Path = strcat(folder_LoadSubwindowData,'DataSubwindowCalcs_30min_Unrestricted'); %path for loading time window data
LoadWindowData_Path = strcat(folder_LoadWindowData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'ThresholdAnalysisData'); %path for saving output data

%% load data and rename variables as necessary
load(LoadFluxLawData_Path); %load 30-minute values
load(LoadRoughnessData_Path); %load Roughness values
load(LoadSubwindowData_Path); %load primary data
load(LoadWindowData_Path); %load 30-minute values
addpath(folder_Functions); %point MATLAB to location of functions

%% get info about measurement intervals and sampling intervals
N_Deltat = length(Deltat_all); %number of measurement intervals
N_deltat = length(deltat_all); %number of sampling intervals

%% get z0 for Sites
z0_Site = z0Re_Q_fit;
sigma_z0_Site = 10.^(sigma_z0Re_ln_Q_fit);

%% get flux law intercept values for Sites
tauit_intercept = tauit_linearfit_all;
sigma_tauit_intercept = tauit_sigma_linearfit_all;

%% get predominant wind direction for each site from mean of 30-minute windows with full transport
theta_Site = zeros(N_Sites,1); %predominant wind for each site
for i = 1:N_Sites
    %30-minute calculation
    theta_Site(i) = mean(theta_all{i}(fQ_all{i}==1)); %wind direction for site
end

%% get z/L from associated 30-minute window
zL_subwindow_all = cell(N_Sites,1);
for i = 1:N_Sites
    zL_subwindow_all{i} = cell(N_Deltat,N_deltat);
    for m = 1:N_Deltat
        for s = 1:N_Deltat
            zL_subwindow_all{i}{m,s} = zL_all{i}(ind_window_subwindow_all{i}{m,s});
        end
    end
end

%% adjust subwindow winds to theta_Site
theta_subwindow_adjusted_all = cell(N_Sites,1);
for i = 1:N_Sites
    theta_subwindow_adjusted_all{i} = cell(N_Deltat,N_deltat);
    for m = 1:N_Deltat
        for s = 1:N_Deltat
            theta_subwindow_adjusted_all{i}{m,s} = theta_subwindow_all{i}{m,s}-theta_Site(i);
        end
    end
end

%% for Jeri, convert fQ NaN values to zeros, assuming that times with no flux data are fQ = 0 (not valid for RG and not necessary for Oceano)
ind_Site = find(strcmp(SiteNames,'Jericoacoara'));
for m = 1:length(Deltat_all)
    for s = 1:length(deltat_all)
        fQ_subwindow_all{ind_Site}{m,s}(isnan(fQ_subwindow_all{ind_Site}{m,s}))=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMARY ANALYSIS BY SITE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize list of all values by site
fQ_bin_avg_all = cell(N_Sites,1); %transport activity - average
fQ_bin_SE_all = cell(N_Sites,1); %transport activity - standard error
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
    sigma_z0 = sigma_z0_Site(i); %get uncertainty of roughness length
    
    %% apply threshold binning analysis function
    [fQ_bin_avg,fQ_bin_SE,...
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio,...
        Chi2nu] = ...
    ThresholdBinning(rho,z0,sigma_z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Sites{i},...
        fQ_subwindow_all, uth_subwindow_all, zU_subwindow_all,...
        theta_subwindow_all,theta_max,zL_subwindow_all,zL_max,...
        StartTime_subwindow_all,timeofday_subwindow_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);
        
    %% assign values by site
    fQ_bin_avg_all{i} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_all{i} = fQ_bin_SE; %transport activity - standard error
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
sigma_z0 = sigma_z0_Site(ind_Site); %get uncertainty of roughness length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% sensitivity analysis for measurement interval
N_Deltat = length(Deltat_all);

%% initialize list of all values by Deltat
fQ_bin_avg_Deltat = cell(N_Deltat,1); %transport activity - average
fQ_bin_SE_Deltat = cell(N_Deltat,1); %transport activity - standard error
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
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,sigma_z0,...
        Deltat_all,deltat_all,...
        Deltat_all(m), deltat_analysis,...
        Sites,Site_sensitivityanalysis,...
        fQ_subwindow_all, uth_subwindow_all, zU_subwindow_all,...
        theta_subwindow_all,theta_max,zL_subwindow_all,zL_max,...
        StartTime_subwindow_all,timeofday_subwindow_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);

    %% assign values by Deltat
    fQ_bin_avg_Deltat{m} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_Deltat{m} = fQ_bin_SE; %transport activity - standard error
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
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,sigma_z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_all(s),...
        Sites,Site_sensitivityanalysis,...
        fQ_subwindow_all, uth_subwindow_all, zU_subwindow_all,...
        theta_subwindow_all,theta_max,zL_subwindow_all,zL_max,...
        StartTime_subwindow_all,timeofday_subwindow_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);

    %% assign values by sampling interval
    fQ_bin_avg_deltat{s} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_deltat{s} = fQ_bin_SE; %transport activity - standard error
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
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,sigma_z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Site_sensitivityanalysis,...
        fQ_subwindow_all, uth_subwindow_all, zU_subwindow_all,...
        theta_subwindow_all,theta_max,zL_subwindow_all,zL_max,...
        StartTime_subwindow_all,timeofday_subwindow_all,...
        diurnalrange_starthour(d),diurnalrange_endhour(d));%,...
        %startdate_analysis,enddate_analysis);

    %% assign values by sampling interval
    fQ_bin_avg_diurnalrange{d} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_diurnalrange{d} = fQ_bin_SE; %transport activity - standard error
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
        uth_bin_avg,uth_bin_SE,...
        ustth_bin_avg,ustth_bin_SE,...
        tauth_bin_avg,tauth_bin_SE,...
        tauft, sigma_tauft,...
        tauit, sigma_tauit,...
        ustitftratio, sigma_ustitftratio] = ...
    ThresholdBinning(rho,z0,sigma_z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Site_sensitivityanalysis,...
        fQ_subwindow_all, uth_subwindow_all, zU_subwindow_all,...
        theta_subwindow_all,theta_max,zL_subwindow_all,zL_max,...
        StartTime_subwindow_all,timeofday_subwindow_all,...
        0,24,...
        daterange_startdate(d),daterange_enddate(d));

    %% assign values by sampling interval
    fQ_bin_avg_daterange{d} = fQ_bin_avg; %transport activity - average
    fQ_bin_SE_daterange{d} = fQ_bin_SE; %transport activity - standard error
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
figure(2); clf; hold on; %initialize plot

%plot average values
for i = 1:N_Sites
    plot(fQ_bin_avg_all{i},tauth_bin_avg_all{i},PlotMarkers_Site{i},'Color',PlotColors_Site{i});
end

%plot SE values
for i = 1:N_Sites
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],tauth_bin_avg_all{i}(k)*[1 1],'Color',PlotColors_Site{i});
        plot(fQ_bin_avg_all{i}(k)*[1 1],tauth_bin_avg_all{i}(k)+tauth_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors_Site{i});
    end
end

%plot fit lines
for i = 1:N_Sites
    plot([0 1],[tauft_all(i) tauit_all(i)],'Color',PlotColors_Site{i});
end

%plot fit intercept values
for i = 1:N_Sites
    plot([0 0],tauft_all(i)+sigma_tauft_all(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',2);
    plot([0.999 0.999],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',2);
end

%annotate plot
legend(SiteNames,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$ (Pa)','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity.png'],'-dpng');



%% plot tau_th versus fQ for measurement interval sensitivity analysis
figure(12); clf; hold on; %initialize plot

%plot average values
for m = 1:N_Deltat
    plot(fQ_bin_avg_Deltat{m},tauth_bin_avg_Deltat{m},PlotMarkers_Sensitivity{m},'Color',PlotColors_Sensitivity{m});
end

%plot SE values
for m = 1:N_Deltat
    N_fQ_bins = length(fQ_bin_avg_Deltat{m});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_Deltat{m}(k)+fQ_bin_SE_Deltat{m}(k)*[-1 1],tauth_bin_avg_Deltat{m}(k)*[1 1],'Color',PlotColors_Sensitivity{m});
        plot(fQ_bin_avg_Deltat{m}(k)*[1 1],tauth_bin_avg_Deltat{m}(k)+tauth_bin_SE_Deltat{m}(k)*[-1 1],'Color',PlotColors_Sensitivity{m});
    end
end

%plot fit lines
for m = 1:N_Deltat
    plot([0 1],[tauft_Deltat(m) tauit_Deltat(m)],'Color',PlotColors_Sensitivity{m});
end

%plot fit intercept values
for m = 1:N_Deltat
    plot([0 0],tauft_Deltat(m)+sigma_tauft_Deltat(m)*[-1 1],'Color',PlotColors_Sensitivity{m},'LineWidth',2);
    plot([0.999 0.999],tauit_Deltat(m)+sigma_tauit_Deltat(m)*[-1 1],'Color',PlotColors_Sensitivity{m},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_Deltat,1);
for m = 1:N_Deltat
    legend_items{m} = ['\Delta t = ',num2str(minutes(Deltat_all(m))),' min.'];
end
legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$ (Pa)','interpreter','latex');
title(Site_sensitivityanalysis);
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
    plot(fQ_bin_avg_deltat{s},tauth_bin_avg_deltat{s},PlotMarkers_Sensitivity{s},'Color',PlotColors_Sensitivity{s});
end

%plot SE values
for s = 1:N_deltat
    N_fQ_bins = length(fQ_bin_avg_deltat{s});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_deltat{s}(k)+fQ_bin_SE_deltat{s}(k)*[-1 1],tauth_bin_avg_deltat{s}(k)*[1 1],'Color',PlotColors_Sensitivity{s});
        plot(fQ_bin_avg_deltat{s}(k)*[1 1],tauth_bin_avg_deltat{s}(k)+tauth_bin_SE_deltat{s}(k)*[-1 1],'Color',PlotColors_Sensitivity{s});
    end
end

%plot fit lines
for s = 1:N_deltat
    plot([0 1],[tauft_deltat(s) tauit_deltat(s)],'Color',PlotColors_Sensitivity{s});
end

%plot fit intercept values
for s = 1:N_deltat
    plot([0 0],tauft_deltat(s)+sigma_tauft_deltat(s)*[-1 1],'Color',PlotColors_Sensitivity{s},'LineWidth',2);
    plot([0.999 0.999],tauit_deltat(s)+sigma_tauit_deltat(s)*[-1 1],'Color',PlotColors_Sensitivity{s},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_deltat,1);
for s = 1:N_deltat
    legend_items{s} = ['\delta{t} = ',num2str(seconds(deltat_all(s)),2),' s'];
end
legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$ (Pa)','interpreter','latex');
title(Site_sensitivityanalysis);
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
    plot(fQ_bin_avg_diurnalrange{d},tauth_bin_avg_diurnalrange{d},PlotMarkers_Sensitivity{d},'Color',PlotColors_Sensitivity{d});
end

%plot SE values
for d = 1:N_diurnalrange
    N_fQ_bins = length(fQ_bin_avg_diurnalrange{d});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_diurnalrange{d}(k)+fQ_bin_SE_diurnalrange{d}(k)*[-1 1],tauth_bin_avg_diurnalrange{d}(k)*[1 1],'Color',PlotColors_Sensitivity{d});
        plot(fQ_bin_avg_diurnalrange{d}(k)*[1 1],tauth_bin_avg_diurnalrange{d}(k)+tauth_bin_SE_diurnalrange{d}(k)*[-1 1],'Color',PlotColors_Sensitivity{d});
    end
end

%plot fit lines
for d = 1:N_diurnalrange
    plot([0 1],[tauft_diurnalrange(d) tauit_diurnalrange(d)],'Color',PlotColors_Sensitivity{d});
end

%plot fit intercept values
for d = 1:N_diurnalrange
    plot([0 0],tauft_diurnalrange(d)+sigma_tauft_diurnalrange(d)*[-1 1],'Color',PlotColors_Sensitivity{d},'LineWidth',2);
    plot([0.999 0.999],tauit_diurnalrange(d)+sigma_tauit_diurnalrange(d)*[-1 1],'Color',PlotColors_Sensitivity{d},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_diurnalrange,1);
for d = 1:N_diurnalrange
    legend_items{d} = [num2str(diurnalrange_starthour(d)),'-',num2str(diurnalrange_endhour(d)),'h'];
end

legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$ (Pa)','interpreter','latex');
title(Site_sensitivityanalysis);
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
    plot(fQ_bin_avg_daterange{d},tauth_bin_avg_daterange{d},PlotMarkers_Sensitivity{d},'Color',PlotColors_Sensitivity{d});
end

%plot SE values
for d = 1:N_daterange
    N_fQ_bins = length(fQ_bin_avg_daterange{d});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_daterange{d}(k)+fQ_bin_SE_daterange{d}(k)*[-1 1],tauth_bin_avg_daterange{d}(k)*[1 1],'Color',PlotColors_Sensitivity{d});
        plot(fQ_bin_avg_daterange{d}(k)*[1 1],tauth_bin_avg_daterange{d}(k)+tauth_bin_SE_daterange{d}(k)*[-1 1],'Color',PlotColors_Sensitivity{d});
    end
end

%plot fit lines
for d = 1:N_daterange
    plot([0 1],[tauft_daterange(d) tauit_daterange(d)],'Color',PlotColors_Sensitivity{d});
end

%plot fit intercept values
for d = 1:N_daterange
    plot([0 0],tauft_daterange(d)+sigma_tauft_daterange(d)*[-1 1],'Color',PlotColors_Sensitivity{d},'LineWidth',2);
    plot([0.999 0.999],tauit_daterange(d)+sigma_tauit_daterange(d)*[-1 1],'Color',PlotColors_Sensitivity{d},'LineWidth',2);
end

%annotate plot
legend_items = cell(N_daterange,1);
for d = 1:N_daterange
    legend_items{d} = [datestr(daterange_startdate(d),'mmm dd'),'-',datestr(daterange_enddate(d),'dd')];
end

legend(legend_items,'Location','NorthEast');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold, $$\tau_{th}$$ (Pa)','interpreter','latex');
title(Site_sensitivityanalysis);
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'threshold_activity_daterange.png'],'-dpng');