%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% Parameters
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,2); %sampling interval for analysis

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
PlotFont = 12;
PlotMarkers_Site = {'s','d','o'};
PlotColors_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
PlotMarkers_Sensitivity = {'s','d','o','<','>','^'};
PlotColors_Sensitivity = {[0.6473 0.7456 0.4188],[0.2116 0.1898 0.5777],[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330]};

%% folders for loading and saving data
folder_LoadFluxLawData = '../../AnalysisData/FluxLaw/'; %folder for loading 30 minute data
folder_LoadRoughnessData = '../../AnalysisData/Misc/'; %folder for loading roughness data
folder_LoadSubwindowData = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
folder_LoadWindowData = '../../AnalysisData/Windowing/'; %folder for loading 30 minute data
folder_Functions = '../Functions/'; %folder with functions
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% paths for loading and saving data
LoadFluxLawData_Path = strcat(folder_LoadFluxLawData,'FluxLawCalcs_30min_Restricted'); %path for 30 minute data
LoadRoughnessData_Path = strcat(folder_LoadRoughnessData,'RoughnessCalcs_30min_Restricted'); %path for 30 minute data
LoadSubwindowData_Path = strcat(folder_LoadSubwindowData,'DataIntervalAveragedSubwindowCalcs_30min_Unrestricted'); %path for loading time window data
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
sigma_tauit_intercept = sigma_tauit_linearfit_all;

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
        for s = 1:N_deltat
            zL_subwindow_all{i}{m,s} = zL_all{i}(ind_window_subwindow_all{i}{m,s});
        end
    end
end

%% adjust subwindow winds to theta_Site
theta_subwindow_adjusted_all = cell(N_Sites,1);
for i = 1:N_Sites
    theta_subwindow_adjusted_all{i} = cell(N_Deltat,N_deltat);
    for m = 1:N_Deltat
        for s = 1:N_deltat
            theta_subwindow_adjusted_all{i}{m,s} = theta_subwindow_all{i}{m,s}-theta_Site(i);
        end
    end
end

%% for Jeri, convert fQ NaN values to zeros, assuming that times with no flux data are fQ = 0 (not valid for RG and not necessary for Oceano)
ind_Site = find(strcmp(SiteNames,'Jericoacoara'));
for m = 1:N_Deltat
    for s = 1:N_deltat
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
tauitftratio_all = zeros(N_Sites,1); %ratio of impact and fluid threshold stresses - best fit
sigma_tauitftratio_all = zeros(N_Sites,1); %ratio of impact and fluid threshold stresses - uncertainty
ustft_all = zeros(N_Sites,1); %inferred fluid threshold shear velocity - best fit
sigma_ustft_all = zeros(N_Sites,1); %inferred fluid threshold shear velocity - uncertainty
ustit_all = zeros(N_Sites,1); %inferred impact threshold shear velociy - best fit
sigma_ustit_all = zeros(N_Sites,1); %inferred impact threshold shear velocity - uncertainty
ustitftratio_all = zeros(N_Sites,1); %ratio of impact and fluid threshold shear velocities - best fit
sigma_ustitftratio_all = zeros(N_Sites,1); %ratio of impact and fluid threshold shear velocities - uncertainty
Chi2nu_all = zeros(N_Sites,1); %normalized Chi2 for threshold versus fQ fit
b_all = zeros(N_Sites,1); %slope for fit
sigma_b_all = zeros(N_Sites,1); %uncertainty in slope for fit

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
        Chi2nu, sigma_b] = ...
    ThresholdBinning(rho,z0,sigma_z0,...
        Deltat_all,deltat_all,...
        Deltat_analysis, deltat_analysis,...
        Sites,Sites{i},...
        fQ_subwindow_all, uth_subwindow_all, zU_subwindow_all,...
        theta_subwindow_all,theta_max,zL_subwindow_all,zL_max,...
        StartTime_subwindow_all,timeofday_subwindow_all);%,...
        %starthour_analysis,endhour_analysis,...
        %startdate_analysis,enddate_analysis);
    
    %get equivalent shear velocities and their uncertainties
    ustft = sqrt(tauft/rho);
    sigma_ustft = (sigma_tauft/2)*sqrt(rho/tauft);
    ustit = sqrt(tauit/rho);
    sigma_ustit = (sigma_tauit/2)*sqrt(rho/tauit);

    
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
    tauitftratio_all(i) = tauit/tauft; %ratio of impact and fluid threshold shear stresses - best fit
    sigma_tauitftratio_all(i) = sqrt((sigma_tauit/tauft)^2+(sigma_tauft*tauit/tauft^2)^2); %ratio of impact and fluid threshold shear stresses - uncertainty
    ustft_all(i) = ustft; %inferred fluid threshold shear velocity - best fit
    sigma_ustft_all(i) = sigma_ustft; %inferred fluid threshold shear velocity - uncertainty
    ustit_all(i) = ustit; %inferred impact threshold shear velocity - best fit
    sigma_ustit_all(i) = sigma_ustit; %inferred impact threshold shear velocity - uncertainty
    ustitftratio_all(i) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_all(i) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
    Chi2nu_all(i) = Chi2nu; %normalized Chi2 for threshold versus fQ fit
    b_all(i) = tauit-tauft; %slope for fit
    sigma_b_all(i) = sigma_b; %uncertainty in slope for fit
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
ustft_daterange = zeros(N_daterange,1); %inferred fluid threshold shear velocity - best fit
tauft_daterange = zeros(N_daterange,1); %inferred fluid threshold stress - best fit
sigma_tauft_daterange = zeros(N_daterange,1); %inferred fluid threshold stress - uncertainty
ustit_daterange = zeros(N_daterange,1); %inferred impact threshold shear velocity - best fit
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
    ustft_daterange(d) = sqrt(tauft/rho); %inferred fluid threshold shear velocity - best fit
    tauft_daterange(d) = tauft; %inferred fluid threshold stress - best fit
    sigma_tauft_daterange(d) = sigma_tauft; %inferred fluid threshold stress - uncertainty
    ustit_daterange(d) = sqrt(tauit/rho); %inferred impact threshold shear velocity - best fit
    tauit_daterange(d) = tauit; %inferred impact threshold stress - best fit
    sigma_tauit_daterange(d) = sigma_tauit; %inferred impact threshold stress - uncertainty
    ustitftratio_daterange(d) = ustitftratio; %ratio of impact and fluid threshold shear velocities - best fit
    sigma_ustitftratio_daterange(d) = sigma_ustitftratio; %ratio of impact and fluid threshold shear velocities - uncertainty
end

%% get wind distribution by date range - Oceano
ind_Site = find(strcmp(Site_sensitivityanalysis,Sites)); %% get index of site for analysis
ind_Deltat = find(Deltat_all==Deltat_analysis); %get measurement interval
ind_deltat = find(deltat_all==deltat_analysis); %get sampling interval
starttime = StartTime_subwindow_all{ind_Site}{ind_Deltat,ind_deltat};
ubar = ubar_subwindow_all{ind_Site}{ind_Deltat,ind_deltat}; %get mean wind speeds
zU = zU_subwindow_all{ind_Site}{ind_Deltat,ind_deltat}; %get anemometer heights
z0 = z0_Site(ind_Site); %get z0
kappa = 0.4; %assume von Karman constant
ubar_daterange = cell(N_daterange,1); %initialize list of ubar for each date range
p_ubar_daterange = cell(N_daterange,1); %initialize pdf of ubar values in bin
ubar_bins_daterange = cell(N_daterange,1); %initialize ubars for each bin
uft_daterange = zeros(N_daterange,1); %initialize list of uft for each date range
uit_daterange = zeros(N_daterange,1); %initialize list of uit for each date range
for d = 1:N_daterange
    ind_date = intersect(find(days(starttime-daterange_startdate(d))>0),find(days(starttime-daterange_enddate(d))<0)); %get indices of date range
    ubar_daterange{d} = ubar(ind_date); %get ubars for this date range
    zU_bar = mean(zU(ind_date)); %get mean zU for this date range
    [N_ubar, ubar_bins] = hist(ubar_daterange{d},20); %get histogram of ubar values
    p_ubar_daterange{d} = N_ubar/(mean(diff(ubar_bins))*length(ubar_daterange{d})); %get pdf of ubar values in bin
    ubar_bins_daterange{d} = ubar_bins; %get ubars for each bin
    uft_daterange(d) = ustft_daterange(d)/kappa*log(zU_bar/z0); %get fluid threshold
    uit_daterange(d) = ustit_daterange(d)/kappa*log(zU_bar/z0); %get impact threshold
end

%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%
save(SaveData_Path,...
    'deltat_all','Deltat_all',...
    'diurnalrange_starthour','diurnalrange_endhour',...
    'daterange_startdate','daterange_enddate',...
    'fQ_bin_avg_*','fQ_bin_SE_*',...
    'tauft_*','sigma_tauft_*',...
    'tauit_*','sigma_tauit_*',...
    'ustitftratio_*','sigma_ustitftratio_*',...
    'Chi2nu_*');

%%%%%%%%%%%%
% PLOTTING %
%%%%%%%%%%%%

%% plot tau_th versus fQ for sites
figure(2); clf; 

%subplot(1,14,1:12);
subplot('Position',[0.09 0.1 0.76 0.88]); hold on;

hold on; %initialize plot

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
    %choose x-location for plotting impact threshold
    if i==2
        x_it = 0.995;
    else
        x_it = 0.999;
    end
    plot([x_it x_it],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',2);        
end

%organize plot
legend(SiteNames,'Location','SouthWest');
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold stress, $$\tau_{th}$$ (Pa)','interpreter','latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
ylim([0.08 0.18]);
ylims = ylim;
text(0.02,ylims(2)-0.005,'(a)','FontSize',PlotFont);
set(gca,'FontSize',PlotFont);

%additional plot for values from flux law fit
%subplot(1,14,14); hold on;
subplot('Position',[0.91 0.1 0.08 0.88]); hold on;
for i = 1:N_Sites
    plot([0 0]+i/2,tauit_intercept(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i}); %plot average values
    plot([0 0]+i/2,tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',2); %plot SE values
end

%organize plot
xlim([0 2]);
ylim(ylims);
ylabel('threshold from stress-flux fit, $$\tau_{th,flux}$$ (Pa)','interpreter','latex');
set(gca,'XTick',[],'XTickLabel',{''},'YTickLabel',{''},'YMinorTick','On','Box','On');
text(0.3,ylims(2)-0.005,'(b)','FontSize',PlotFont);
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
print([folder_Plots,'threshold_activity.png'],'-dpng');
print([folder_Plots,'threshold_activity.tif'],'-dtiff','-r600');

%% plot fQ versus fD

%get indices for comparison
ind_Deltat = find(Deltat_analysis == Deltat_all);
ind_deltat = find(deltat_analysis == deltat_all);

%get values for each site
fD_plot = cell(N_Sites,1);
fQ_plot = cell(N_Sites,1);
for i = 1:N_Sites
    ind_intermittent = find(fD_subwindow_all{i}{ind_Deltat,ind_deltat}>0 & fD_subwindow_all{i}{ind_Deltat,ind_deltat}<1);
    fD_plot{i} = fD_subwindow_all{i}{ind_Deltat,ind_deltat}(ind_intermittent);
    fQ_plot{i} = fQ_subwindow_all{i}{ind_Deltat,ind_deltat}(ind_intermittent);
end

%initialize plot
figure(11); clf;

%compare fD and fQ
subplot(2,1,1); hold on;
for i = 1:N_Sites
    plot(fD_plot{i},fQ_plot{i},PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',4);
end

%organize plot
legend(SiteNames,'Location','SouthEast');
xlabel('transport detection rate, $$f_D$$','interpreter','latex');
ylabel('transport activity, $$f_Q$$','interpreter','latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlim([0 1]);
ylim([0 1]);
text(0.02,0.95,'(a)','FontSize',PlotFont);
set(gca,'FontSize',PlotFont);
set(gca, 'LooseInset', get(gca,'TightInset'));

%look at fQ/fD ratio
subplot(2,1,2); hold on;
for i = 1:N_Sites
    plot(fD_plot{i},fQ_plot{i}./fD_plot{i},PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',4);
end

%organize plot
xlabel('transport detection rate, $$f_D$$','interpreter','latex');
ylabel('ratio of activity to detection rate, $$f_Q/f_D$$','interpreter','latex');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlim([0 1]);
text(0.9,1.55,'(b)','FontSize',PlotFont);
set(gca,'FontSize',PlotFont);
set(gca, 'LooseInset', get(gca,'TightInset'));

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[5 7],'PaperPosition',[0 0 5 7],'PaperPositionMode','Manual');
print([folder_Plots,'detection_activity.png'],'-dpng');
print([folder_Plots,'detection_activity.tif'],'-dtiff','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot sensitivity analyses %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12); clf;

%% plot tau_th versus fQ for sampling interval sensitivity analysis
%subplot(2,2,1); hold on; %initialize plot
subplot('Position',[0.09 0.54 0.43 0.45]); hold on;

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

%annotate and format plot
xlim([0 1]);
ylim([0.08 0.135]);
legend_items = cell(N_deltat,1);
for s = 1:N_deltat
    legend_items{s} = ['\delta{t} = ',num2str(seconds(deltat_all(s)),2),' s'];
end
legend(legend_items,'Location','NorthEast');
text(0.03,0.133,'(a)','FontSize',PlotFont);
% xlabel('transport activity, $$f_Q$$','interpreter','latex');
set(gca,'XTickLabel',[]);
ylabel('effective threshold stress, $$\tau_{th}$$ (Pa)','interpreter','latex');
% title(Site_sensitivityanalysis);
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
% print([folder_Plots,'threshold_activity_samplinginterval.png'],'-dpng');
% print([folder_Plots,'threshold_activity_samplinginterval.tif'],'-dtiff','-r600');

%% plot tau_th versus fQ for measurement interval sensitivity analysis
%subplot(2,2,2); hold on; %initialize plot
subplot('Position',[0.56 0.54 0.43 0.45]); hold on;

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

%annotate and format plot
xlim([0 1]);
ylim([0.08 0.135]);
legend_items = cell(N_Deltat,1);
for m = 1:N_Deltat
    legend_items{m} = ['\Delta t = ',num2str(minutes(Deltat_all(m))),' min.'];
end
legend(legend_items,'Location','SouthWest');
text(0.03,0.133,'(b)','FontSize',PlotFont);
% xlabel('transport activity, $$f_Q$$','interpreter','latex');
set(gca,'XTickLabel',[]);
% ylabel('effective threshold stress, $$\tau_{th}$$ (Pa)','interpreter','latex');
set(gca,'YTickLabel',[]);
% title(Site_sensitivityanalysis);
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
% print([folder_Plots,'threshold_activity_measurementinterval.png'],'-dpng');
% print([folder_Plots,'threshold_activity_measurementinterval.tif'],'-dtiff','-r600');


%% plot tau_th versus fQ for diurnal period sensitivity analysis
%subplot(2,2,3); hold on; %initialize plot
subplot('Position',[0.09 0.07 0.43 0.45]); hold on;

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

%annotate and format plot
xlim([0 1]);
ylim([0.08 0.135]);
legend_items = cell(N_diurnalrange,1);
for d = 1:N_diurnalrange
    legend_items{d} = [num2str(diurnalrange_starthour(d)),'-',num2str(diurnalrange_endhour(d)),'h'];
end
legend(legend_items,'Location','NorthEast');
text(0.03,0.133,'(c)','FontSize',PlotFont);
xlabel('transport activity, $$f_Q$$','interpreter','latex');
ylabel('effective threshold stress, $$\tau_{th}$$ (Pa)','interpreter','latex');
% title(Site_sensitivityanalysis);
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
% print([folder_Plots,'threshold_activity_diurnalrange.png'],'-dpng');
% print([folder_Plots,'threshold_activity_diurnalrange.tif'],'-dtiff','-r600');

%% plot tau_th versus fQ for date period sensitivity analysis
%subplot(2,2,4); hold on; %initialize plot
subplot('Position',[0.56 0.07 0.43 0.45]); hold on;

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

%annotate and format plot
xlim([0 1]);
ylim([0.08 0.135]);
legend_items = cell(N_daterange,1);
for d = 1:N_daterange
    legend_items{d} = [datestr(daterange_startdate(d),'mmm dd'),'-',datestr(daterange_enddate(d),'dd')];
end
legend(legend_items,'Location','NorthEast');
text(0.03,0.133,'(d)','FontSize',PlotFont);
xlabel('transport activity, $$f_Q$$','interpreter','latex');
% ylabel('effective threshold stress, $$\tau_{th}$$ (Pa)','interpreter','latex');
set(gca,'YTickLabel',[]);
% title(Site_sensitivityanalysis);
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
% print([folder_Plots,'threshold_activity_daterange.png'],'-dpng');
% print([folder_Plots,'threshold_activity_daterange.tif'],'-dtiff','-r600');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[7 7],'PaperPosition',[0 0 7 7],'PaperPositionMode','Manual');
print([folder_Plots,'threshold_activity_sensitivityanalyses.png'],'-dpng');
print([folder_Plots,'threshold_activity_sensitivityanalyses.tif'],'-dtiff','-r600');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot thresholds versus known relationships %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_min_plot = 0.1; %mm
d_max_plot = 1.0; %mm
d_plot = logspace(log10(d_min_plot),log10(d_max_plot),100);
rho_p = 2650; %kg/m^3
g = 9.8; %m/s^2
P = 1e5; %Pa
T = 293.15; %K (=20 C)
rho_a = 1.2; %kg/m^3
mu = 0.00751*(T/291.15)^(1.5)/(T+120); %the viscosity in kg/(m s) from the Sutherland relation

%fluid threshold predictions
tauft_Bagnold1941 = 0.1^2*rho_p*g*(d_plot/1000); %Pa (Kok et al. 2012, Eq. 2.5)
tauft_ShaoLu2000 = 0.111^2*(rho_p*g*(d_plot/1000) + 0.29./(d_plot*rho_a)); %Pa (Kok et al., 2012, Eq. 2.8)

%fluid threshold prediction - Iversen and White (1982) (also Greeley and Iversen (1985) - from Jasper's code)
ustft_IversenWhite1982 = zeros(size(d_plot)); %initialize ust values
for k = 1:length(d_plot)
    i = 1; %initial value for loop
    error = 1; %initial error for iterative calculation
    ustft_calc = 2; %initial guess at fluid threshold, Greeley and Iversen
    d = d_plot(k)*1e-3; %grain size for calculation, m
    while (error > 0.000001) %calculating the fluid threshold, following Eqs. (3.5) and (3.12) in Iversen and White (1982)        
        i = i+1; %iterate i
        Re_IW = rho_a*d*ustft_calc/mu; %Reynolds number for Greeley and Iversen calculation
        if (Re_IW > 10)
            A = 0.120*sqrt(1+0.0000006/(rho_p*g*d^(2.5)))*(1-0.0858*exp(-0.0617*(Re_IW-10))); %calculation of parameter A for Iversen and White (1982) equation
        elseif (Re_IW > 0.3)
            A = 0.129*sqrt(1+0.0000006/(rho_p*g*d^(2.5)))/sqrt(1.928*Re_IW^(0.092)-1); %calculation of parameter A for Iversen and White (1982) equation
        elseif (Re_IW > 0.03)
            A = 0.200*sqrt(1+0.0000006/(rho_p*g*d^(2.5)))/sqrt(1+2.5*Re_IW); %calculation of parameter A for Iversen and White (1982) equation
        else
            error('Reynolds number too small for Greeley & Iversen scheme');
        end
        ustft_last = ustft_calc; %record last ust_ft
        ustft_calc = A*sqrt((rho_p/rho_a)*g*d); %compute new ust_ft
        error = abs(ustft_calc/ustft_last-1); %compute difference
    end
    ustft_IversenWhite1982(k) = ustft_calc; %record final ust_ft value from iterative calculation
end
tauft_IversenWhite1982 = rho_a*ustft_IversenWhite1982.^2; %compute fluid threshold shear stress for Greeley and Iversen

%impact threshold predictions
tauit_Bagnold1941 = 0.082^2*rho_p*g*(d_plot/1000); %Pa (Kok et al. 2012, Eq. 2.5)
d_KokEtAl2012 = 1e-3*[50;57;65;80;100;125;150;200;250;325;400;500;650;800;1000;1150;1350;1650;2000]; %mm (Kok et al. 2012, Fig 21a)
ustit_KokEtAl2012 = [0.1914325;0.1587225;0.139529375;0.13898875;0.14131125;0.1490075;0.156135;0.178005;0.19810625;0.22675;0.25335375;0.28567;0.326585625;0.36428125;0.40729875;0.436010625;0.47370875;0.521250625;0.5729775]; %m/s (Kok et al. 2012, Fig 21a)
tauit_KokEtAl2012 = rho_a.*ustit_KokEtAl2012.^2; %Pa (Kok et al. 2012, Fig 21a)

% initialize plot 
figure(13); clf; hold on;

% observations - fluid threshold
for i = 1:N_Sites
    plot(d50_Site(i),tauft_all(i),['b',PlotMarkers_Site{i}],'MarkerSize',10,'LineWidth',1); %fluid threshold
    %plot(d50_Site(i),tauft_all(i),['k',PlotMarkers_Site{i}],'MarkerSize',10,'MarkerFaceColor','b'); %fluid threshold
end

% predictions - fluid threshold
plot(d_plot,tauft_Bagnold1941,'b-','LineWidth',1);
plot(d_plot,tauft_IversenWhite1982,'b--','LineWidth',1);
plot(d_plot,tauft_ShaoLu2000,'b-.','LineWidth',1);

% observations - impact threshold
for i = 1:N_Sites
    plot(d50_Site(i),tauit_all(i),['g',PlotMarkers_Site{i}],'MarkerSize',10,'LineWidth',1); %impact threshold
    %plot(d50_Site(i),tauit_all(i),['k',PlotMarkers_Site{i}],'MarkerSize',10,'MarkerFaceColor','g'); %impact threshold
end

% predictions - impact threshold
plot(d_plot,tauit_Bagnold1941,'g-','LineWidth',1);
plot(d_KokEtAl2012,tauit_KokEtAl2012,'g--','LineWidth',1);

% observations - error bars
for i = 1:N_Sites
    plot(d50_Site(i)*[1 1],tauft_all(i)+sigma_tauft_all(i)*[-1 1],'b','LineWidth',1); %sigma ft
    plot(d50_Site(i)*[1 1],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'g','LineWidth',1); %sigma it
end

%create legend
h_legend = legend('\tau_{ft,Rancho Guadalupe}','\tau_{ft,Jericoacoara}','\tau_{ft,Oceano}',...
    '\tau_{ft,Bagnold (1941)}','\tau_{ft,Iversen & White (1982)}','\tau_{ft,Shao & Lu (2000)}',...
    '\tau_{it,Jericoacoara}','\tau_{it,Rancho Guadalupe}','\tau_{it,Oceano}',...
    '\tau_{it,Bagnold (1941)}','\tau_{it,Kok et al. (2012)}',...
    'Location','NorthWest');
set(h_legend,'FontSize',12);

% annotate plot
xlabel('particle diameter, $$d$$ (mm)','interpreter','latex');
ylabel('shear stress, $$\tau$$ (Pa)','interpreter','latex');
set(gca,'FontSize',PlotFont);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'XTick',0.1:0.1:1,'XTickLabel',{'0.1','0.2','','','0.5','','','','','1.0'});
xlim([0.1 1.0]);
ylim([0 0.3]);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[7.5 6],'PaperPosition',[0 0 7.5 6],'PaperPositionMode','Manual');
print([folder_Plots,'threshold_prediction.png'],'-dpng');
print([folder_Plots,'threshold_prediction.tif'],'-dtiff','-r600');

%%%%%%%%%%%%%%%%%%%%%%%
%% wind distributions %
%%%%%%%%%%%%%%%%%%%%%%%
figure(14); clf; hold on;
plot(ubar_bins_daterange{2},p_ubar_daterange{2},'-','Color','k','LineWidth',1); %2nd date range wind distribution
plot(uit_daterange(2)*[1 1],[0 0.5],'-','Color','g','LineWidth',1); %2nd date range impact threshold
plot(uft_daterange(2)*[1 1],[0 0.5],'-','Color','b','LineWidth',1); %2nd date range fluid threshold
plot(ubar_bins_daterange{3},p_ubar_daterange{3},'--','Color','k','LineWidth',1); %3rd date range wind distribution
plot(uit_daterange(3)*[1 1],[0 0.5],'--','Color','g','LineWidth',1); %3rd date range impact threshold
plot(uft_daterange(3)*[1 1],[0 0.5],'--','Color','b','LineWidth',1); %3rd date range fluid threshold

%create legend
legend_items = cell(6,1);
legend_items{1} = [datestr(daterange_startdate(2),'mmm dd'),'-',datestr(daterange_enddate(2),'dd'),': \Phi_{u,1 minute}'];
legend_items{2} = [datestr(daterange_startdate(2),'mmm dd'),'-',datestr(daterange_enddate(2),'dd'),': u_{it}'];
legend_items{3} = [datestr(daterange_startdate(2),'mmm dd'),'-',datestr(daterange_enddate(2),'dd'),': u_{ft}'];
legend_items{4} = [datestr(daterange_startdate(3),'mmm dd'),'-',datestr(daterange_enddate(3),'dd'),': \Phi_{u,1 minute}'];
legend_items{5} = [datestr(daterange_startdate(3),'mmm dd'),'-',datestr(daterange_enddate(3),'dd'),': u_{it}'];
legend_items{6} = [datestr(daterange_startdate(3),'mmm dd'),'-',datestr(daterange_enddate(3),'dd'),': u_{ft}'];
h_legend = legend(legend_items,'Location','NorthWest');
set(h_legend,'FontSize',PlotFont);

%annotate plot
xlim([0 10]);
xlabel('$$u$$, streamwise wind speed (m/s)','Interpreter','Latex','FontSize',PlotFont);
ylabel('probability density (m$$^{-1}$$ s)','Interpreter','Latex','FontSize',PlotFont);
set(gca,'FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperSize',[6.5 4],'PaperPosition',[0 0 6.5 4],'PaperPositionMode','Manual');
print([folder_Plots,'wind_distribution_daterange.png'],'-dpng');
print([folder_Plots,'wind_distribution_daterange.tif'],'-dtiff','-r600');