%% initialize
clearvars;

%%
%%%%%%%%%%%%%%
% PARAMETERS %
%%%%%%%%%%%%%%

%% physical parameters
rho_s = 2650*1e3; %sediment density, g/m^3
A_w = 30*0.6*(1e-3)^2; %area of Wenglor, m^2

%% set reference grain size
dref_type = 'd50';
%dref_type = 'dmodal';

%% set range of trustworthy grain sizes
d_min = 0.13; %minimum trustworthy grain size
d_max = 1; %maximum trustworthy grain size

%% Information about sites, clusters, and cluster names
SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};
% rho_a = [1.16, 1.22, 1.22]; %air density by Site
% tauth_Site = [0.135, 0.110, 0.094]; %threshold stress by Site
% sigma_tauth_Site = [0.015, 0.021, 0.006]; %uncertainty in threshold stress by Site
% ustth_Site = sqrt(tauth_Site./rho_a); %threshold shear velocity by Site

%% set threshold values for cluster, based on thresholds paper
tauth_Cluster = [0.111, 0.110, 0.083, 0.091]; %shear stress threshold by Cluster
sigma_tauth_Cluster = [0.002, 0.002, 0.001, 0.001]; %shear stress uncertainty by Cluster
rho_a = [1.16, 1.22, 1.22, 1.22]; %air density by Cluster
ustth_Cluster = sqrt(tauth_Cluster./rho_a); %threshold shear velocity by Cluster
sigma_ustth_Cluster = sigma_tauth_Cluster./(2*rho_a.*ustth_Cluster); %uncertainty in threshold shear velocity by Cluster

N_Sites = length(SiteNames);
Cluster_StartDate = {...
    datetime(2014,11,13),...
    datetime(2015,3,23),...
    [datetime(2015,5,15); datetime(2015,5,23)]
    %[datetime(2015,5,15); datetime(2015,5,23); datetime(2015,6,1)]
    };
Cluster_EndDate = {...
    datetime(2014,11,20),...
    datetime(2015,3,24),...
    [datetime(2015,5,22); datetime(2015,6,4)]
    %[datetime(2015,5,22); datetime(2015,5,31); datetime(2015,6,4)]
    };
N_Cluster_Site = cellfun(@length,Cluster_StartDate); %number of surface grain size clusters for each site (Jeri, RG, Oceano);
N_Cluster = sum(N_Cluster_Site); %total number of clusters

%generate names of clusters
ClusterNames = cell(N_Cluster,1); %names of clusters
for i = 1:N_Sites
    for j = 1:N_Cluster_Site(i)
        if i == 1
            ind_Cluster = j;
        else
            ind_Cluster = sum(N_Cluster_Site(1:i-1))+j;
        end
        if N_Cluster_Site(i) == 1
            ClusterNames{ind_Cluster} = SiteNames{i};
        else
            ClusterNames{ind_Cluster} = [SiteNames{i},' ',int2str(j)];
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR BINNING BY GRAIN SIZE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% partial flux grain-size fixed bin values
binning_type = 'fixed';
N_bins = 6; %number of bins
d_bin_min = d_min; %minimum bin (mm)
d_bin_max = d_max; %maximum bin (mm)
d_bin_edges = logspace(log10(d_bin_min),log10(d_bin_max),N_bins+1);
d_bin_lower = d_bin_edges(1:N_bins); %lower edge of bin (mm)
d_bin_upper = d_bin_edges(2:N_bins+1); %upper edge of bin (mm)
d_bin_mid = geomean([d_bin_lower; d_bin_upper],1); %midpoint of bin (mm)

% %% define bins in terms of values with respect to dref (i.e., d/dref)
% binning_type = 'dhat';
% N_bins = 6; %number of bins
% dhat_bin_min = 0.37; %bottom of lower bin -- associated with 0.13 mm minimum grain diameter for BSNE collection efficiency
% dhat_bin_edges = logspace(log10(dhat_bin_min),log10(1/dhat_bin_min),N_bins+2); %edges of bins -- log-spaced
% dhat_bin_lower = dhat_bin_edges(1:N_bins); %d/dref values for size bins - lower limit
% dhat_bin_upper = dhat_bin_edges(2:N_bins+1); %d/dref values for size bins - upper limit
% dhat_bin_mid = geomean([dhat_bin_lower; dhat_bin_upper],1); %d/dref values for size bins - midpoint
% % 
% %% define bins in terms of CDF of airborne values
% binning_type = 'CDFa';
% N_bins = 9; %number of bins
% CDFa_bin_min = 0.05; %minimum bin
% CDFa_bin_max = 0.95; %maximum bin
% CDFa_bin_edges = linspace(CDFa_bin_min,CDFa_bin_max,N_bins+1);
% CDFa_bin_lower = CDFa_bin_edges(1:N_bins); %CDF values for size bins - lower limit
% CDFa_bin_upper = CDFa_bin_edges(2:N_bins+1); %CDF values for size bins - upper limit
% CDFa_bin_mid = mean([CDFa_bin_lower; CDFa_bin_upper],1); %CDF values for size bins - midpoint
% f_d_surface_CDFa = CDFa_bin_upper - CDFa_bin_lower; %fractions of distribution for CDF values

% %% define bins in terms of CDF of surface values
% binning_type = 'CDFs';
% N_bins = 9; %number of bins
% CDFs_bin_min = 0.05; %minimum bin
% CDFs_bin_max = 0.95; %maximum bin
% CDFs_bin_edges = linspace(CDFs_bin_min,CDFs_bin_max,N_bins+1);
% CDFs_bin_lower = CDFs_bin_edges(1:N_bins); %CDF values for size bins - lower limit
% CDFs_bin_upper = CDFs_bin_edges(2:N_bins+1); %CDF values for size bins - upper limit
% CDFs_bin_mid = mean([CDFs_bin_lower; CDFs_bin_upper],1); %CDF values for size bins - midpoint
% f_d_surface_CDFs = CDFs_bin_upper - CDFs_bin_lower; %fractions of distribution for CDF values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR BINNING BY SHEAR VELOCITY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ust_bin_min = 0.2:0.05:0.55; %u* for bottom of shear velocity bins
ust_bin_max = 0.25:0.05:0.6; %u* for top of shear velocity bins
ust_bin_mid = mean([ust_bin_min; ust_bin_max],1); %midpoint u* for shear velocity bins
N_ust_bins = length(ust_bin_min);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR BINNING BY NORMALIZED SHEAR VELOCITY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ustnorm_min = 1.0;
ustnorm_max = 2.0;
N_ustnorm_bins = 4;
% ustnorm_bin_edges = [0,logspace(log10(ustnorm_min),log10(ustnorm_max),N_ustnorm_bins)];
% ustnorm_bin_min = ustnorm_bin_edges(1:N_ustnorm_bins); %u*norm for bottom of normalized shear velocity bins
% ustnorm_bin_max = ustnorm_bin_edges(2:N_ustnorm_bins+1); %u*norm for top of normalized shear velocity bins
ustnorm_bin_min = [0,1,1.25,1.6]; %set limits manually
ustnorm_bin_max = [1,1.25,1.6,2]; %set limits manually

%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR FITTING %
%%%%%%%%%%%%%%%%%%%%%%%%%%
taunorm_max_Qhat_tau_fit = 2; %maximum normalized shear stress for fitting Qhat versus tau

%%
%%%%%%%%%%%%%%%%
% DATA SOURCES %
%%%%%%%%%%%%%%%%

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
folder_SaltationData = '../../AnalysisData/Windowing/'; %folder for saltation flux data
GrainSizeData_Path = strcat(folder_GrainSizeData,'GrainSizeData'); %path for loading mean grain size data
SaltationFluxData_Path = strcat(folder_SaltationData,'DataWindowCalcs_30min_Restricted'); %path for loading saltation data
folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%
% PLOTTING INFO %
%%%%%%%%%%%%%%%%%

%% plotting information
PlotFont = 10; %font for labels
LineWidth_Plot = 1; %width of lines
Marker_Cluster = {'s','d','o','p','h','^','v','>'}; %markers for Clusters
Color_Cluster = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]}; %colors for sites / clusters
Label_Cluster = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'}; %markers for Clusters
Marker_bin = {'o','s','d','^','v','>','<','p','h','x','*','+'}; %markers for bins
Color_bin = cell(N_bins,1);
for i = 1:N_bins
    Color_bin{i} = [0, (i-1)/(N_bins-1), 1-(i-1)/(N_bins-1)];
end
Color_ust_bin = cell(N_ust_bins,1);
for i = 1:N_ust_bins
    Color_ust_bin{i} = [(i-1)/(N_ust_bins-1), 0, 1-(i-1)/(N_ust_bins-1)];
end
Color_ustnorm_bin = cell(N_ustnorm_bins,1);
for i = 1:N_ustnorm_bins
    Color_ustnorm_bin{i} = [(i-1)/(N_ustnorm_bins-1), 0, 1-(i-1)/(N_ustnorm_bins-1)];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(GrainSizeData_Path); %grain size data
load(SaltationFluxData_Path); %saltation data
addpath(folder_Functions); %point MATLAB to location of functions


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AGGREGATE SURFACE AND AIRBORNE SAMPLES INTO CLUSTERS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize information for surface sample clusters
d_surface_mid_Cluster = cell(N_Cluster,1); %diameter of midpoints of grain size bins
d_surface_lower_Cluster = cell(N_Cluster,1); %diameter of lower edges of grain size bins
d_surface_upper_Cluster = cell(N_Cluster,1); %diameter of upper edges of grain size bins
dlogd_surface_Cluster = cell(N_Cluster,1); %normalized diameter range of grain size bins
dV_bar_surface_Cluster = cell(N_Cluster,1); %mean volume fraction of surface samples
dVdlogd_bar_surface_Cluster = cell(N_Cluster,1); %mean differential volume fraction of surface samples
d10_bar_surface_Cluster = zeros(N_Cluster,1); %mean d10 of surface samples
d50_bar_surface_Cluster = zeros(N_Cluster,1); %mean d50 of surface samples
d90_bar_surface_Cluster = zeros(N_Cluster,1); %mean d90 of surface samples
dbar_bar_surface_Cluster = zeros(N_Cluster,1); %mean dbar of surface samples
dmodal_bar_surface_Cluster = zeros(N_Cluster,1); %mean dmodal of surface samples
d10_sigma_surface_Cluster = zeros(N_Cluster,1); %d10 uncertainty of surface samples
d50_sigma_surface_Cluster = zeros(N_Cluster,1); %d50 uncertainty of surface samples
d90_sigma_surface_Cluster = zeros(N_Cluster,1); %d90 uncertainty of surface samples
dbar_sigma_surface_Cluster = zeros(N_Cluster,1); %dbar uncertainty of surface samples
dmodal_sigma_surface_Cluster = zeros(N_Cluster,1); %dmodal uncertainty of surface samples

%% initialize information for airborne sample clusters

%grain size bins
d_airborne_mid_Cluster = cell(N_Cluster,1); %diameter of midpoints of grain size bins
d_airborne_lower_Cluster = cell(N_Cluster,1); %diameter of lower edges of grain size bins
d_airborne_upper_Cluster = cell(N_Cluster,1); %diameter of upper edges of grain size bins
dlogd_airborne_Cluster = cell(N_Cluster,1); %normalized diameter range of grain size bins

%information for individual samples (not separated into profiles)
dV_airborne_Cluster = cell(N_Cluster,1); %volume fraction of airborne individual sample
dVdlogd_airborne_Cluster = cell(N_Cluster,1); %differential volume fraction of airborne individual sample
d10_airborne_Cluster = cell(N_Cluster,1); %d10 of airborne samples
d50_airborne_Cluster = cell(N_Cluster,1); %d50 of airborne samples
d90_airborne_Cluster = cell(N_Cluster,1); %d90 of airborne samples
dbar_airborne_Cluster = cell(N_Cluster,1); %mean d of airborne samples
z_airborne_Cluster = cell(N_Cluster,1); %height associated with individual sample

%information for profiles
dV_profile_airborne_Cluster = cell(N_Cluster,1); %volume fraction of airborne individual sample
dVdlogd_profile_airborne_Cluster = cell(N_Cluster,1); %differential volume fraction of airborne individual sample
d10_profile_airborne_Cluster = cell(N_Cluster,1); %d10 of airborne sample profiles
d50_profile_airborne_Cluster = cell(N_Cluster,1); %d50 of airborne sample profiles
d90_profile_airborne_Cluster = cell(N_Cluster,1); %d90 of airborne sample profiles
dbar_profile_airborne_Cluster = cell(N_Cluster,1); %mean d of airborne sample profiles
q_profile_Cluster = cell(N_Cluster,1); %partial flux associated with individual sample
sigma_q_profile_Cluster = cell(N_Cluster,1); %partial flux uncertainty associated with individual sample
z_profile_Cluster = cell(N_Cluster,1); %height associated with individual sample
sigma_z_profile_Cluster = cell(N_Cluster,1); %height uncertainty associated with individual sample
Q_profile_Cluster = cell(N_Cluster,1); %total flux associated with sample profile
sigma_Q_profile_Cluster = cell(N_Cluster,1); %total flux associated with sample profile
zq_profile_Cluster = cell(N_Cluster,1); %e-folding height associated with sample profile
tau_profile_Cluster = cell(N_Cluster,1); %shear stress associated with sample profile
ust_profile_Cluster = cell(N_Cluster,1); %shear velocity associated with sample profile
taunorm_profile_Cluster = cell(N_Cluster,1); %normalized shear stress (tau/tauth) associated with sample profile
ustnorm_profile_Cluster = cell(N_Cluster,1); %normalized shear velocity (ust/ustth) associated with sample profile

dV_profilebar_airborne_Cluster = cell(N_Cluster,1); %profile mean volume fraction of airborne samples
dVdlogd_profilebar_airborne_Cluster = cell(N_Cluster,1); %profile mean differential volume fraction of airborne samples
d10_profilebar_airborne_Cluster = cell(N_Cluster,1); %profile mean d10 of airborne samples
d50_profilebar_airborne_Cluster = cell(N_Cluster,1); %profile mean d50 of airborne samples
d90_profilebar_airborne_Cluster = cell(N_Cluster,1); %profile mean d90 of airborne samples
dbar_profilebar_airborne_Cluster = cell(N_Cluster,1); %profile mean d of airborne sample profiles
d10_profilesigma_airborne_Cluster = cell(N_Cluster,1); %profile std d10 of airborne samples
d50_profilesigma_airborne_Cluster = cell(N_Cluster,1); %profile std d50 of airborne samples
d90_profilesigma_airborne_Cluster = cell(N_Cluster,1); %profile std d90 of airborne samples
dbar_profilesigma_airborne_Cluster = cell(N_Cluster,1); %profile std mean d of airborne sample profiles

dV_bar_airborne_Cluster = cell(N_Cluster,1); %cluster mean volume fraction of airborne samples in cluster
dVdlogd_bar_airborne_Cluster = cell(N_Cluster,1); %cluster mean differential volume fraction of airborne samples in cluster
d10_bar_airborne_Cluster = zeros(N_Cluster,1); %mean d10 of airborne samples
d50_bar_airborne_Cluster = zeros(N_Cluster,1); %mean d50 of airborne samples
d90_bar_airborne_Cluster = zeros(N_Cluster,1); %mean d90 of airborne samples
dbar_bar_airborne_Cluster = zeros(N_Cluster,1); %mean dbar of airborne samples
d10_sigma_airborne_Cluster = zeros(N_Cluster,1); %d10 uncertainty of airborne samples
d50_sigma_airborne_Cluster = zeros(N_Cluster,1); %d50 uncertainty of airborne samples
d90_sigma_airborne_Cluster = zeros(N_Cluster,1); %d90 uncertainty of airborne samples
dbar_sigma_airborne_Cluster = zeros(N_Cluster,1); %dbar uncertainty of airborne samples

%% go through sites
for i = 1:N_Sites
    
    %% surface samples
    % load relevant information about surface samples
    GrainSize_surface = GrainSizeData_all{i}.GrainSize_Surface; %surface sample array
    Date_surface = [GrainSize_surface.Date]; %dates of surface samples
    N_surface = length(GrainSize_surface); %number of surface samples
    
    %get surface size bins from first surface sample
    d_surface_mid = [GrainSize_surface(1).gsd(2:end-1).Sizeclass_mid_mm];
    d_surface_lower = [GrainSize_surface(1).gsd(2:end-1).Sizeclass_lower_mm];
    d_surface_upper = [GrainSize_surface(1).gsd(2:end-1).Sizeclass_upper_mm];
    dlogd_surface = log(d_surface_upper) - log(d_surface_lower);
    N_d_surface = length(d_surface_mid); %number of grain size bins
    
    %get each surface size distribution
    dV_surface = zeros(N_surface,N_d_surface); % initialize volume fraction in bin
    dVdlogd_surface = zeros(N_surface,N_d_surface); %initialize dV/dlogd in bin
    for j = 1:N_surface %go through each surface sample
        dV_surface(j,:) = [GrainSize_surface(j).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
        dVdlogd_surface(j,:) = dV_surface(j,:)./dlogd_surface; %by normalized volume
    end
    
    %% airborne samples
    % load relevant information about airborne samples
    GrainSize_airborne_samples = GrainSizeData_all{i}.GrainSize_BSNE; %airborne sample array
    Name_airborne_samples = {GrainSize_airborne_samples.NameBSNE}; %BSNE names of airborne samples
    d10_airborne_samples = [GrainSize_airborne_samples.d_10_mm]; %airborne sample d10
    d50_airborne_samples = [GrainSize_airborne_samples.d_50_mm]; %airborne sample d50
    d90_airborne_samples = [GrainSize_airborne_samples.d_90_mm]; %airborne sample d90
    Date_airborne_samples = [GrainSize_airborne_samples.Date]'; %initial time of airborne sample
    StartTime_airborne_samples = [GrainSize_airborne_samples.StartTime]'; %initial time of airborne sample
    EndTime_airborne_samples = [GrainSize_airborne_samples.EndTime]'; %end time of airborne sample
    N_airborne_samples = length(GrainSize_airborne_samples); %number of airborne samples
    
    % get saltation flux data for site
    Flux_BSNE = FluxBSNE_all{i}; %saltation flux array for BSNEs 
    N_Flux = length(Flux_BSNE); %number of BSNE flux intervals
    
    %get airborne size bins from first airborne sample
    d_airborne_mid = [GrainSize_airborne_samples(1).gsd(2:end-1).Sizeclass_mid_mm];
    d_airborne_lower = [GrainSize_airborne_samples(1).gsd(2:end-1).Sizeclass_lower_mm];
    d_airborne_upper = [GrainSize_airborne_samples(1).gsd(2:end-1).Sizeclass_upper_mm];
    dlogd_airborne = log(d_airborne_upper) - log(d_airborne_lower);
    N_d_airborne = length(d_airborne_mid); %number of grain size bins
    
    %get each airborne size distribution
    dV_airborne_samples = zeros(N_airborne_samples,N_d_airborne); % initialize volume fraction in bin
    dVdlogd_airborne_samples = zeros(N_airborne_samples,N_d_airborne); %initialize dV/dlogd in bin
    for j = 1:N_airborne_samples %go through each airborne sample
        dV_airborne_samples(j,:) = [GrainSize_airborne_samples(j).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
        dVdlogd_airborne_samples(j,:) = dV_airborne_samples(j,:)./dlogd_airborne; %by normalized volume
    end
    
    %get mean d of each airborne size distribution
    dbar_airborne_samples = zeros(N_airborne_samples,1); %initialize list of dbar
    for j = 1:N_airborne_samples %go through each airborne sample
        dbar_airborne_samples(j,:) = exp(sum(dV_airborne_samples(j,:).*log(d_airborne_mid))); %compute weighted geometric mean of particle size bins
    end
    
    %get height of each airborne sample
    z_airborne_samples = zeros(N_airborne_samples,1)*NaN; %height of each airborne sample
    for j = 1:N_airborne_samples %go through each airborne sample
        ind_BSNE = find([Flux_BSNE.StartTime] >= StartTime_airborne_samples(j) & [Flux_BSNE.StartTime] <= EndTime_airborne_samples(j), 1);
        ind_z = find(strcmp(Flux_BSNE(ind_BSNE).name,Name_airborne_samples{j}));
        if ~isempty(ind_z)
            z_airborne_samples(j) = Flux_BSNE(ind_BSNE).z.z(ind_z);
        end
    end
        
    %% combine airborne grain size distributions for each BSNE collection interval  
    
    %get dates for profiles
    Date_profile_airborne = [Flux_BSNE.Date];
    
    %initialize profile values
    dV_profile_airborne = cell(N_Flux,1); %volume fraction of airborne individual sample
    dVdlogd_profile_airborne = cell(N_Flux,1); %differential volume fraction of airborne individual sample
    d10_profile_airborne = cell(N_Flux,1); %d10 of airborne individual sample
    d50_profile_airborne = cell(N_Flux,1); %d50 of airborne individual sample
    d90_profile_airborne = cell(N_Flux,1); %d90 of airborne individual sample
    dbar_profile_airborne = cell(N_Flux,1); %dbar of airborne individual sample
    q_profile_airborne = cell(N_Flux,1); %partial flux associated with individual sample
    sigma_q_profile_airborne = cell(N_Flux,1); %partial flux uncertainty associated with individual sample
    z_profile_airborne = cell(N_Flux,1); %height associated with individual sample
    sigma_z_profile_airborne = cell(N_Flux,1); %height uncertainties associated with individual sample
    Q_profile_airborne = zeros(N_Flux,1); %total flux associated with sample profile
    sigma_Q_profile_airborne = zeros(N_Flux,1); %total flux uncertainty associated with sample profile
    zq_profile_airborne = zeros(N_Flux,1); %e-folding height associated with profile
    tau_profile_airborne = zeros(N_Flux,1); %shear stress associated with sample profile
    ust_profile_airborne = zeros(N_Flux,1); %shear velocity associated with sample profile
    
    %initialize profile mean values
    dV_profilebar_airborne = zeros(N_Flux,N_d_airborne); %profile mean volume fraction of airborne samples
    dVdlogd_profilebar_airborne = zeros(N_Flux,N_d_airborne); %profile mean differential volume fraction of airborne samples
    d10_profilebar_airborne = zeros(N_Flux,1)*NaN; %profile mean d10 of airborne samples
    d50_profilebar_airborne = zeros(N_Flux,1)*NaN; %profile mean d50 of airborne samples
    d90_profilebar_airborne = zeros(N_Flux,1)*NaN; %profile mean d90 of airborne samples
    dbar_profilebar_airborne = zeros(N_Flux,1)*NaN; %profile mean dbar of airborne samples
    d10_profilesigma_airborne = zeros(N_Flux,1)*NaN; %profile std d10 of airborne samples
    d50_profilesigma_airborne = zeros(N_Flux,1)*NaN; %profile std d50 of airborne samples
    d90_profilesigma_airborne = zeros(N_Flux,1)*NaN; %profile std d90 of airborne samples
    dbar_profilesigma_airborne = zeros(N_Flux,1)*NaN; %profile std dbar of airborne samples
    
    %go through each BSNE collection interval
    for j = 1:N_Flux
        
        %get total flux for BSNE
        Q_profile_airborne(j) = Flux_BSNE(j).Q.Q; %total flux associated with sample profile
        sigma_Q_profile_airborne(j) = Flux_BSNE(j).Q.sigma_Q; %total flux associated with sample profile
        zq_profile_airborne(j) = Flux_BSNE(j).z.zq; %e-folding height associated with sample profile
       
        %get values for BSNE
        q_BSNE = Flux_BSNE(j).qz.qz; %get partial flux for BSNE
        sigma_q_BSNE = Flux_BSNE(j).qz.sigma; %get partial flux uncertainty for BSNE
        z_BSNE = Flux_BSNE(j).z.z; %get height for BSNE
        sigma_z_BSNE = Flux_BSNE(j).z.sigma_z; %get height uncertainty for BSNE
        Name_BSNE = Flux_BSNE(j).name; %get names of BSNEs in profile
        StartTime_BSNE = Flux_BSNE(j).StartTime; %get start time for BSNE interval
        EndTime_BSNE = Flux_BSNE(j).EndTime; %get end time for BSNE interval
        
        %get wind associated with BSNE flux  
        ind_wind = find(StartTimes_all{i}>=StartTime_BSNE & EndTimes_all{i}<=EndTime_BSNE); %indices of wind for BSNE interval
        tau_profile_airborne(j) = mean(tauRe_all{i}(ind_wind)); %mean shear stress for interval
        ust_profile_airborne(j) = mean(ustRe_all{i}(ind_wind)); %mean shear velocity for interval
        
        %get indices of grain size samples associated with BSNE flux
        ind_samples_BSNE = find(StartTime_airborne_samples<=StartTime_BSNE & EndTime_airborne_samples>StartTime_BSNE); %get indices within GrainSize_airborne corresponding to time interval
        N_samples_profile = length(ind_samples_BSNE); %number of airborne samples associated with flux
        
        %get each airborne size distribution
        dV_profile_airborne{j} = dV_airborne_samples(ind_samples_BSNE,:); %by volume
        dVdlogd_profile_airborne{j} = dVdlogd_airborne_samples(ind_samples_BSNE,:); %by normalized volume
        d10_profile_airborne{j} = d10_airborne_samples(ind_samples_BSNE); %d10 of airborne individual sample
        d50_profile_airborne{j} = d50_airborne_samples(ind_samples_BSNE); %d50 of airborne individual sample
        d90_profile_airborne{j} = d90_airborne_samples(ind_samples_BSNE); %d90 of airborne individual sample
        dbar_profile_airborne{j} = dbar_airborne_samples(ind_samples_BSNE); %dbar of airborne individual sample
        
        %initialize matrices of airborne samples
        q_profile_airborne{j} = zeros(N_samples_profile,1); %partial fluxes for grain size samples
        sigma_q_profile_airborne{j} = zeros(N_samples_profile,1); %partial flux uncertainties for grain size samples
        z_profile_airborne{j} = zeros(N_samples_profile,1); %heights for grain size samples
        sigma_z_profile_airborne{j} = zeros(N_samples_profile,1); %height uncertainties for grain size samples
        
        %get each airborne size distribution
        for k = 1:N_samples_profile
            try
                ind_k = find(strcmp(Name_BSNE,Name_airborne_samples(ind_samples_BSNE(k)))); %get index for BSNE
                q_profile_airborne{j}(k) = q_BSNE(ind_k); %get flux associated with BSNE
                sigma_q_profile_airborne{j}(k) = sigma_q_BSNE(ind_k); %get flux uncertainty associated with BSNE
                z_profile_airborne{j}(k) = z_BSNE(ind_k); %get height associated with BSNE
                sigma_z_profile_airborne{j}(k) = sigma_z_BSNE(ind_k); %get height uncertainty associated with BSNE
            catch
                q_profile_airborne{j}(k) = NaN;
                sigma_q_profile_airborne{j}(k) = NaN;
                z_profile_airborne{j}(k) = NaN;
                sigma_z_profile_airborne{j}(k) = NaN;
            end           
        end
        
        %eliminate NaN smaples
        ind_hasdata = find(~isnan(q_profile_airborne{j}));
        N_samples_profile = length(ind_hasdata);
        q_profile_airborne{j} = q_profile_airborne{j}(ind_hasdata);
        sigma_q_profile_airborne{j} = sigma_q_profile_airborne{j}(ind_hasdata);
        z_profile_airborne{j} = z_profile_airborne{j}(ind_hasdata);
        sigma_z_profile_airborne{j} = sigma_z_profile_airborne{j}(ind_hasdata);
        dV_profile_airborne{j} = dV_profile_airborne{j}(ind_hasdata,:);
        dVdlogd_profile_airborne{j} = dVdlogd_profile_airborne{j}(ind_hasdata,:);
        
        %create flux weighting matrix to compute flux-weighted profile mean grain-size distribution
        dV_weighted_matrix = zeros(N_samples_profile,N_d_airborne); %create matrix of dV's multiplied by weights
        for k = 1:N_samples_profile
            weight_k = q_profile_airborne{j}(k)/sum(q_profile_airborne{j}); %get weight for sample
            dV_weighted_matrix(k,:) = dV_profile_airborne{j}(k,:).*weight_k; %compute weighted values for row
        end
        dV_profilebar_airborne(j,:) = sum(dV_weighted_matrix); %add weighted rows to get flux-weighted grain-size distribution           
        dVdlogd_profilebar_airborne(j,:) = dV_profilebar_airborne(j,:)./dlogd_airborne; %convert to dVdlogd         
        
        %get profile mean airborne size distribution reference sizes
        if sum(q_profile_airborne{j})>0 && ~isempty(q_profile_airborne{j}) %compute only if there is flux associated with grain size
            [d10, d50, d90] = ReferenceGrainSizes(dV_profilebar_airborne(j,:), d_airborne_lower, d_airborne_upper);
            d10_profilebar_airborne(j) = d10; %d10 of airborne profile
            d50_profilebar_airborne(j) = d50; %d50 of airborne profile
            d90_profilebar_airborne(j) = d90; %d90 of airborne profile
            dbar_profilebar_airborne(j) = exp(sum(dV_profilebar_airborne(j,:).*log(d_airborne_mid))); %compute weighted geometric mean of particle size bins
            d10_profilesigma_airborne(j) = std(d10_profile_airborne{j}); %std of d10 in profile
            d50_profilesigma_airborne(j) = std(d50_profile_airborne{j}); %std of d50 in profile
            d90_profilesigma_airborne(j) = std(d90_profile_airborne{j}); %std of d90 in profile
            dbar_profilesigma_airborne(j) = std(dbar_profile_airborne{j}); %std of dbar in profile
        end
    end
    
    %% keep only intervals with data
    ind_hasdata = find(~isnan(d10_profilebar_airborne)); %indices of profiles with data
    Date_profile_airborne = Date_profile_airborne(ind_hasdata); %dates for profiles  
    dV_profile_airborne = dV_profile_airborne(ind_hasdata); %volume fraction of airborne individual sample
    dVdlogd_profile_airborne = dVdlogd_profile_airborne(ind_hasdata); %differential volume fraction of airborne individual sample
    d10_profile_airborne = d10_profile_airborne(ind_hasdata); %d10 of airborne individual sample
    d50_profile_airborne = d50_profile_airborne(ind_hasdata); %d50 of airborne individual sample
    d90_profile_airborne = d90_profile_airborne(ind_hasdata); %d90 of airborne individual sample
    dbar_profile_airborne = dbar_profile_airborne(ind_hasdata); %dbar of airborne individual sample
    q_profile_airborne = q_profile_airborne(ind_hasdata); %partial flux associated with individual sample
    sigma_q_profile_airborne = sigma_q_profile_airborne(ind_hasdata); %partial flux uncertainty associated with individual sample
    z_profile_airborne = z_profile_airborne(ind_hasdata); %height associated with individual sample
    sigma_z_profile_airborne = sigma_z_profile_airborne(ind_hasdata); %height uncertainty associated with individual sample
    Q_profile_airborne = Q_profile_airborne(ind_hasdata); %total flux associated with sample profile
    sigma_Q_profile_airborne = sigma_Q_profile_airborne(ind_hasdata); %total flux uncertainty associated with sample profile
    zq_profile_airborne = zq_profile_airborne(ind_hasdata); %e-folding height associated with sample profile
    tau_profile_airborne = tau_profile_airborne(ind_hasdata); %mean shear stress for interval
    ust_profile_airborne = ust_profile_airborne(ind_hasdata); %mean shear velocity for interval
    dV_profilebar_airborne = dV_profilebar_airborne(ind_hasdata,:); %profile mean volume fraction of airborne samples
    dVdlogd_profilebar_airborne = dVdlogd_profilebar_airborne(ind_hasdata,:); %profile mean differential volume fraction of airborne samples
    d10_profilebar_airborne = d10_profilebar_airborne(ind_hasdata); %profile mean d10 of airborne samples
    d50_profilebar_airborne = d50_profilebar_airborne(ind_hasdata); %profile mean d50 of airborne samples
    d90_profilebar_airborne = d90_profilebar_airborne(ind_hasdata); %profile mean d90 of airborne samples
    dbar_profilebar_airborne = dbar_profilebar_airborne(ind_hasdata); %profile mean dbar of airborne samples
    d10_profilesigma_airborne = d10_profilesigma_airborne(ind_hasdata); %profile std d10 of airborne samples
    d50_profilesigma_airborne = d50_profilesigma_airborne(ind_hasdata); %profile std d50 of airborne samples
    d90_profilesigma_airborne = d90_profilesigma_airborne(ind_hasdata); %profile std d90 of airborne samples
    dbar_profilesigma_airborne = dbar_profilesigma_airborne(ind_hasdata); %profile std dbar of airborne samples
    
    %% aggregate data into clusters   
    for j = 1:N_Cluster_Site(i)
               
        %% get index of cluster
        if i == 1
            ind_Cluster = j;
        else
            ind_Cluster = sum(N_Cluster_Site(1:i-1))+j;
        end
               
        %% calculations for surface samples
        
        %add surface grain size bins into cluster
        d_surface_mid_Cluster{ind_Cluster} = d_surface_mid;
        d_surface_lower_Cluster{ind_Cluster} = d_surface_lower;
        d_surface_upper_Cluster{ind_Cluster} = d_surface_upper;
        dlogd_surface_Cluster{ind_Cluster} = dlogd_surface;
        
        %get indices of surface samples associated with cluster
        ind_surface_Cluster = intersect(find(Date_surface>=Cluster_StartDate{i}(j)),find(Date_surface<=Cluster_EndDate{i}(j)));
        
        %compute mean dV and dVdlogd of surface samples for cluster        
        dV_bar_surface_Cluster{ind_Cluster} = mean(dV_surface(ind_surface_Cluster,:),1);
        dVdlogd_bar_surface_Cluster{ind_Cluster} = dV_bar_surface_Cluster{ind_Cluster}./dlogd_surface;      

        %get mean d10, d50, d90 of surface samples
        [d10, d50, d90] = ReferenceGrainSizes(dV_bar_surface_Cluster{ind_Cluster}, d_surface_lower, d_surface_upper);
        d10_bar_surface_Cluster(ind_Cluster) = d10;
        d50_bar_surface_Cluster(ind_Cluster) = d50;
        d90_bar_surface_Cluster(ind_Cluster) = d90;
        dbar_bar_surface_Cluster(ind_Cluster) =  exp(sum(dV_bar_surface_Cluster{ind_Cluster}.*log(d_surface_mid))); %compute weighted geometric mean of particle size bins
        d10_sigma_surface_Cluster(ind_Cluster) = std([GrainSize_surface(ind_surface_Cluster).d_10_mm])./sqrt(length(ind_surface_Cluster));
        d50_sigma_surface_Cluster(ind_Cluster) = std([GrainSize_surface(ind_surface_Cluster).d_50_mm])./sqrt(length(ind_surface_Cluster));
        d90_sigma_surface_Cluster(ind_Cluster) = std([GrainSize_surface(ind_surface_Cluster).d_90_mm])./sqrt(length(ind_surface_Cluster));
        %dbar_sigma_surface_Cluster(ind_Cluster) = ???
   
        %get mean dbar of surface samples and its uncertainty
        dbar_bar = exp(sum(dV_bar_surface_Cluster{ind_Cluster}.*log(d_surface_mid))); %dbar of average distribution
        dbar = zeros(length(ind_surface_Cluster),1);
        for k = 1:length(ind_surface_Cluster)
            dbar(k) = mean(exp(sum(dV_surface(ind_surface_Cluster(k),:).*log(d_surface_mid)))); %dbar of average distribution
        end
        sigma_dbar = std(dbar);
        dbar_bar_surface_Cluster(ind_Cluster) = dbar_bar; %mean dbar of surface samples
        dbar_sigma_surface_Cluster(ind_Cluster) = sigma_dbar; %dbar uncertainty of surface samples
        
        %get mean dmodal of surface samples and its uncertainty
        dmodal_bar = d_surface_mid(dV_bar_surface_Cluster{ind_Cluster}==max(dV_bar_surface_Cluster{ind_Cluster})); %dmodal of average distribution
        dmodal = zeros(length(ind_surface_Cluster),1);
        for k = 1:length(ind_surface_Cluster)
            dmodal(k) = mean(d_surface_mid(dV_surface(ind_surface_Cluster(k),:)==max(dV_surface(ind_surface_Cluster(k),:)))); %dmodal of average distribution
        end
        sigma_dmodal = std(dmodal);
        dmodal_bar_surface_Cluster(ind_Cluster) = dmodal_bar; %mean dmodal of surface samples
        dmodal_sigma_surface_Cluster(ind_Cluster) = sigma_dmodal; %dmodal uncertainty of surface samples
        
        %% airborne samples
        
        %add airborne grain size bins into cluster
        d_airborne_mid_Cluster{ind_Cluster} = d_airborne_mid;
        d_airborne_lower_Cluster{ind_Cluster} = d_airborne_lower;
        d_airborne_upper_Cluster{ind_Cluster} = d_airborne_upper;
        dlogd_airborne_Cluster{ind_Cluster} = dlogd_airborne;
                
        %get indices of individual airborne samples associated with cluster
        ind_airborne_samples_Cluster = intersect(find(Date_airborne_samples>=Cluster_StartDate{i}(j)),find(Date_airborne_samples<=Cluster_EndDate{i}(j)));
        N_airborne_samples_Cluster = length(ind_airborne_samples_Cluster); %number of airborne samples in cluster
        
        %group information about individual airborne samples
        dV_airborne_Cluster{ind_Cluster} = dV_airborne_samples(ind_airborne_samples_Cluster,:); %volume fraction of airborne individual sample
        dVdlogd_airborne_Cluster{ind_Cluster} = dVdlogd_airborne_samples(ind_airborne_samples_Cluster,:); %differential volume fraction of airborne individual sample
        d10_airborne_Cluster{ind_Cluster} = d10_airborne_samples(ind_airborne_samples_Cluster); %d10 of airborne sample profiles
        d50_airborne_Cluster{ind_Cluster} = d50_airborne_samples(ind_airborne_samples_Cluster); %d50 of airborne sample profiles
        d90_airborne_Cluster{ind_Cluster} = d90_airborne_samples(ind_airborne_samples_Cluster); %d90 of airborne sample profiles
        dbar_airborne_Cluster{ind_Cluster} = dbar_airborne_samples(ind_airborne_samples_Cluster); %dbar of airborne sample profiles

        z_airborne_Cluster{ind_Cluster} = z_airborne_samples(ind_airborne_samples_Cluster); %height associated with individual sample
                
        %get indices of airborne profiles associated with cluster
        ind_airborne_Cluster = intersect(find(Date_profile_airborne>=Cluster_StartDate{i}(j)),find(Date_profile_airborne<=Cluster_EndDate{i}(j)));
        N_airborne_Cluster = length(ind_airborne_Cluster); %number of airborne samples in cluster
               
        %assign profile values within cluster
        dV_profile_airborne_Cluster{ind_Cluster} = dV_profile_airborne(ind_airborne_Cluster); %volume fraction of airborne individual sample
        dVdlogd_profile_airborne_Cluster{ind_Cluster} = dVdlogd_profile_airborne(ind_airborne_Cluster); %differential volume fraction of airborne individual sample
        d10_profile_airborne_Cluster{ind_Cluster} = d10_profile_airborne(ind_airborne_Cluster); %profile mean d10 of airborne samples
        d50_profile_airborne_Cluster{ind_Cluster} = d50_profile_airborne(ind_airborne_Cluster); %profile mean d50 of airborne samples
        d90_profile_airborne_Cluster{ind_Cluster} = d90_profile_airborne(ind_airborne_Cluster); %profile mean d90 of airborne samples
        dbar_profile_airborne_Cluster{ind_Cluster} = dbar_profile_airborne(ind_airborne_Cluster); %profile mean dbar of airborne samples
        q_profile_Cluster{ind_Cluster} = q_profile_airborne(ind_airborne_Cluster); %partial fluxes associated with individual sample
        sigma_q_profile_Cluster{ind_Cluster} = sigma_q_profile_airborne(ind_airborne_Cluster); %partial flux uncertainties associated with individual sample
        z_profile_Cluster{ind_Cluster} = z_profile_airborne(ind_airborne_Cluster); %heights associated with individual sample
        sigma_z_profile_Cluster{ind_Cluster} = sigma_z_profile_airborne(ind_airborne_Cluster); %height uncertainties associated with individual sample
        Q_profile_Cluster{ind_Cluster} = Q_profile_airborne(ind_airborne_Cluster); %total flux associated with sample profile
        sigma_Q_profile_Cluster{ind_Cluster} = sigma_Q_profile_airborne(ind_airborne_Cluster); %total flux uncertainty associated with sample profile
        zq_profile_Cluster{ind_Cluster} = zq_profile_airborne(ind_airborne_Cluster); %e-folding height associated with sample profile
        tau_profile_Cluster{ind_Cluster} = tau_profile_airborne(ind_airborne_Cluster); %shear stress associated with sample profile
        ust_profile_Cluster{ind_Cluster} = ust_profile_airborne(ind_airborne_Cluster); %shear velocity associated with sample profile
        taunorm_profile_Cluster{ind_Cluster} = tau_profile_airborne(ind_airborne_Cluster)/tauth_Cluster(ind_Cluster); %normalized shear stress associated with sample profile
        ustnorm_profile_Cluster{ind_Cluster} = ust_profile_airborne(ind_airborne_Cluster)/ustth_Cluster(ind_Cluster); %normalized shear velocity associated with sample profile
        
        dV_profilebar_airborne_Cluster{ind_Cluster} = dV_profilebar_airborne(ind_airborne_Cluster,:); %profile mean volume fraction of airborne samples
        dVdlogd_profilebar_airborne_Cluster{ind_Cluster} = dV_profilebar_airborne_Cluster{ind_Cluster}./dlogd_airborne; %profile mean differential volume fraction of airborne samples
        d10_profilebar_airborne_Cluster{ind_Cluster} = d10_profilebar_airborne(ind_airborne_Cluster); %profile mean d10 of airborne samples
        d50_profilebar_airborne_Cluster{ind_Cluster} = d50_profilebar_airborne(ind_airborne_Cluster); %profile mean d50 of airborne samples
        d90_profilebar_airborne_Cluster{ind_Cluster} = d90_profilebar_airborne(ind_airborne_Cluster); %profile mean d90 of airborne samples
        dbar_profilebar_airborne_Cluster{ind_Cluster} = dbar_profilebar_airborne(ind_airborne_Cluster); %profile mean dbar of airborne samples
        d10_profilesigma_airborne_Cluster{ind_Cluster} = d10_profilesigma_airborne(ind_airborne_Cluster); %profile std d10 of airborne samples
        d50_profilesigma_airborne_Cluster{ind_Cluster} = d50_profilesigma_airborne(ind_airborne_Cluster); %profile std d50 of airborne samples
        d90_profilesigma_airborne_Cluster{ind_Cluster} = d90_profilesigma_airborne(ind_airborne_Cluster); %profile std d90 of airborne samples
        dbar_profilesigma_airborne_Cluster{ind_Cluster} = dbar_profilesigma_airborne(ind_airborne_Cluster); %profile std dbar of airborne samples
        
        %compute mean dV and dVdlogd of airborne samples for cluster        
        dV_bar_airborne_Cluster{ind_Cluster} = mean(dV_profilebar_airborne(ind_airborne_Cluster,:),1);
        dVdlogd_bar_airborne_Cluster{ind_Cluster} = dV_bar_airborne_Cluster{ind_Cluster}./dlogd_airborne; 
%        dVdlogd_bar_airborne_Cluster{ind_Cluster} = mean(dVdlogd_profilebar_airborne(ind_airborne_Cluster,:),1); 
                
        %get mean d10, d50, d90 of mean dV
        [d10, d50, d90] = ReferenceGrainSizes(dV_bar_airborne_Cluster{ind_Cluster}, d_airborne_lower, d_airborne_upper);
        d10_bar_airborne_Cluster(ind_Cluster) = d10;
        d50_bar_airborne_Cluster(ind_Cluster) = d50;
        d90_bar_airborne_Cluster(ind_Cluster) = d90;
        dbar_bar_airborne_Cluster(ind_Cluster) =  exp(sum(dV_bar_airborne_Cluster{ind_Cluster}.*log(d_airborne_mid))); %compute weighted geometric mean of particle size bins
        d10_sigma_airborne_Cluster(ind_Cluster) = std(d10_profilebar_airborne_Cluster{ind_Cluster});
        d50_sigma_airborne_Cluster(ind_Cluster) = std(d50_profilebar_airborne_Cluster{ind_Cluster});
        d90_sigma_airborne_Cluster(ind_Cluster) = std(d90_profilebar_airborne_Cluster{ind_Cluster});
        dbar_sigma_airborne_Cluster(ind_Cluster) = std(dbar_profilebar_airborne_Cluster{ind_Cluster});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE BOUNDS OF BINS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_bin_lower_Cluster = zeros(N_Cluster,N_bins); %lower edge diameter of bin (mm)
d_bin_upper_Cluster = zeros(N_Cluster,N_bins); %upper edge diameter of bin (mm)
d_bin_mid_Cluster = zeros(N_Cluster,N_bins); %midpoint diameter of bin (mm)
dhat_bin_lower_Cluster = zeros(N_Cluster,N_bins); %lower edge d/dref of bin
dhat_bin_upper_Cluster = zeros(N_Cluster,N_bins); %upper edge d/dref of bin
dhat_bin_mid_Cluster = zeros(N_Cluster,N_bins); %midpoint d/dref of bin
CDFa_bin_lower_Cluster = zeros(N_Cluster,N_bins); %airborne CDF values for size bins - lower limit
CDFa_bin_upper_Cluster = zeros(N_Cluster,N_bins); %airborne CDF values for size bins - upper limit
CDFa_bin_mid_Cluster = zeros(N_Cluster,N_bins); %airborne CDF values for size bins - midpoint
CDFs_bin_lower_Cluster = zeros(N_Cluster,N_bins); %surface CDF values for size bins - lower limit
CDFs_bin_upper_Cluster = zeros(N_Cluster,N_bins); %surface CDF values for size bins - upper limit
CDFs_bin_mid_Cluster = zeros(N_Cluster,N_bins); %surface CDF values for size bins - midpoint

if strcmp(binning_type,'fixed')==1
    for i = 1:N_Cluster
        d_bin_lower_Cluster(i,:) = d_bin_lower; %lower edge diameter of bin (mm)
        d_bin_upper_Cluster(i,:) = d_bin_upper; %upper edge diameter of bin (mm)
        d_bin_mid_Cluster(i,:) = d_bin_mid; %midpoint diameter of bin (mm)
        if strcmp(dref_type,'d50')==1
            dhat_bin_lower_Cluster(i,:) = d_bin_lower/d50_bar_surface_Cluster(i); %lower edge d/d50 of bin
            dhat_bin_upper_Cluster(i,:) = d_bin_upper/d50_bar_surface_Cluster(i); %upper edge d/d50 of bin
            dhat_bin_mid_Cluster(i,:) = d_bin_mid/d50_bar_surface_Cluster(i); %midpoint d/d50 of bin
        elseif strcmp(dref_type,'dmodal')==1
            dhat_bin_lower_Cluster(i,:) = d_bin_lower/dmodal_bar_surface_Cluster(i); %lower edge d/dmodal of bin
            dhat_bin_upper_Cluster(i,:) = d_bin_upper/dmodal_bar_surface_Cluster(i); %upper edge d/dmodal of bin
            dhat_bin_mid_Cluster(i,:) = d_bin_mid/dmodal_bar_surface_Cluster(i); %midpoint d/dmodal of bin
        end
        for j = 1:N_bins
            CDFa_bin_lower_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_lower(j)); %airborne CDF values for size bins - lower limit
            CDFa_bin_upper_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_upper(j)); %airborne CDF values for size bins - upper limit
            CDFa_bin_mid_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_mid(j)); %airborne CDF values for size bins - midpoint
            CDFs_bin_lower_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_lower(j)); %surface CDF values for size bins - lower limit
            CDFs_bin_upper_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_upper(j)); %surface CDF values for size bins - upper limit
            CDFs_bin_mid_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_mid(j)); %surface CDF values for size bins - midpoint
        end
    end
elseif strcmp(binning_type,'dhat')==1
    for i = 1:N_Cluster
        if strcmp(dref_type,'d50')==1
            d_bin_lower_Cluster(i,:) = dhat_bin_lower*d50_bar_surface_Cluster(i); %lower edge diameter of bin (mm)
            d_bin_upper_Cluster(i,:) = dhat_bin_upper*d50_bar_surface_Cluster(i); %upper edge diameter of bin (mm)
            d_bin_mid_Cluster(i,:) = dhat_bin_mid*d50_bar_surface_Cluster(i); %midpoint diameter of bin (mm)
        elseif strcmp(dref_type,'dmodal')==1
            d_bin_lower_Cluster(i,:) = dhat_bin_lower*dmodal_bar_surface_Cluster(i); %lower edge diameter of bin (mm)
            d_bin_upper_Cluster(i,:) = dhat_bin_upper*dmodal_bar_surface_Cluster(i); %upper edge diameter of bin (mm)
            d_bin_mid_Cluster(i,:) = dhat_bin_mid*dmodal_bar_surface_Cluster(i); %midpoint diameter of bin (mm)
        end
        dhat_bin_lower_Cluster(i,:) = dhat_bin_lower; %lower edge d/d50 of bin
        dhat_bin_upper_Cluster(i,:) = dhat_bin_upper; %upper edge d/d50 of bin
        dhat_bin_mid_Cluster(i,:) = dhat_bin_mid; %midpoint d/d50 of bin
        for j = 1:N_bins
            CDFa_bin_lower_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_lower_Cluster(i,j)); %airborne CDF values for size bins - lower limit
            CDFa_bin_upper_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_upper_Cluster(i,j)); %airborne CDF values for size bins - upper limit
            CDFa_bin_mid_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_mid_Cluster(i,j)); %airborne CDF values for size bins - midpoint
            CDFs_bin_lower_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_lower_Cluster(i,j)); %surface CDF values for size bins - lower limit
            CDFs_bin_upper_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_upper_Cluster(i,j)); %surface CDF values for size bins - upper limit
            CDFs_bin_mid_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_mid_Cluster(i,j)); %surface CDF values for size bins - midpoint
        end
    end    
elseif strcmp(binning_type,'CDFa')==1
    for i = 1:N_Cluster
        for j = 1:N_bins
            d_bin_lower_Cluster(i,j) = ReferenceGrainSizes_arbitrary(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, CDFa_bin_lower(j)); %lower edge diameter of bin (mm)
            d_bin_upper_Cluster(i,j) = ReferenceGrainSizes_arbitrary(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, CDFa_bin_upper(j)); %upper edge diameter of bin (mm)
            d_bin_mid_Cluster(i,j) = ReferenceGrainSizes_arbitrary(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, CDFa_bin_mid(j)); %midpoint diameter of bin (mm)
        end
        if strcmp(dref_type,'d50')==1
            dhat_bin_lower_Cluster(i,:) = d_bin_lower_Cluster(i,:)/d50_bar_surface_Cluster(i); %lower edge d/d50 of bin
            dhat_bin_upper_Cluster(i,:) = d_bin_upper_Cluster(i,:)/d50_bar_surface_Cluster(i); %upper edge d/d50 of bin
            dhat_bin_mid_Cluster(i,:) = d_bin_mid_Cluster(i,:)/d50_bar_surface_Cluster(i); %midpoint d/d50 of bin
        elseif strcmp(dref_type,'dmodal')==1
            dhat_bin_lower_Cluster(i,:) = d_bin_lower_Cluster(i,:)/dmodal_bar_surface_Cluster(i); %lower edge d/dmodal of bin
            dhat_bin_upper_Cluster(i,:) = d_bin_upper_Cluster(i,:)/dmodal_bar_surface_Cluster(i); %upper edge d/dmodal of bin
            dhat_bin_mid_Cluster(i,:) = d_bin_mid_Cluster(i,:)/dmodal_bar_surface_Cluster(i); %midpoint d/dmodal of bin
        end
        CDFa_bin_lower_Cluster(i,:) = CDFa_bin_lower; %airborne CDF values for size bins - lower limit
        CDFa_bin_upper_Cluster(i,:) = CDFa_bin_upper; %airborne CDF values for size bins - upper limit
        CDFa_bin_mid_Cluster(i,:) = CDFa_bin_mid; %airborne CDF values for size bins - midpoint
        for j = 1:N_bins
            CDFs_bin_lower_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_lower_Cluster(i,j)); %surface CDF values for size bins - lower limit
            CDFs_bin_upper_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_upper_Cluster(i,j)); %surface CDF values for size bins - upper limit
            CDFs_bin_mid_Cluster(i,j) = CDF_GrainSize(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_mid_Cluster(i,j)); %surface CDF values for size bins - midpoint
        end
    end
elseif strcmp(binning_type,'CDFs')==1
    for i = 1:N_Cluster
        for j = 1:N_bins
            d_bin_lower_Cluster(i,j) = ReferenceGrainSizes_arbitrary(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, CDFs_bin_lower(j)); %lower edge diameter of bin (mm)
            d_bin_upper_Cluster(i,j) = ReferenceGrainSizes_arbitrary(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, CDFs_bin_upper(j)); %upper edge diameter of bin (mm)
            d_bin_mid_Cluster(i,j) = ReferenceGrainSizes_arbitrary(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, CDFs_bin_mid(j)); %midpoint diameter of bin (mm)
        end
        if strcmp(dref_type,'d50')==1
            dhat_bin_lower_Cluster(i,:) = d_bin_lower_Cluster(i,:)/d50_bar_surface_Cluster(i); %lower edge d/d50 of bin
            dhat_bin_upper_Cluster(i,:) = d_bin_upper_Cluster(i,:)/d50_bar_surface_Cluster(i); %upper edge d/d50 of bin
            dhat_bin_mid_Cluster(i,:) = d_bin_mid_Cluster(i,:)/d50_bar_surface_Cluster(i); %midpoint d/d50 of bin
        elseif strcmp(dref_type,'dmodal')==1
            dhat_bin_lower_Cluster(i,:) = d_bin_lower_Cluster(i,:)/dmodal_bar_surface_Cluster(i); %lower edge d/dmodal of bin
            dhat_bin_upper_Cluster(i,:) = d_bin_upper_Cluster(i,:)/dmodal_bar_surface_Cluster(i); %upper edge d/dmodal of bin
            dhat_bin_mid_Cluster(i,:) = d_bin_mid_Cluster(i,:)/dmodal_bar_surface_Cluster(i); %midpoint d/dmodal of bin
        end
        for j = 1:N_bins
            CDFa_bin_lower_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_lower_Cluster(i,j)); %airborne CDF values for size bins - lower limit
            CDFa_bin_upper_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_upper_Cluster(i,j)); %airborne CDF values for size bins - upper limit
            CDFa_bin_mid_Cluster(i,j) = CDF_GrainSize(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_mid_Cluster(i,j)); %airborne CDF values for size bins - midpoint
        end
        CDFs_bin_lower_Cluster(i,:) = CDFs_bin_lower; %surface CDF values for size bins - lower limit
        CDFs_bin_upper_Cluster(i,:) = CDFs_bin_upper; %surface CDF values for size bins - upper limit
        CDFs_bin_mid_Cluster(i,:) = CDFs_bin_mid; %surface CDF values for size bins - midpoint
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FRACTION OF SAMPLE IN EACH BIN AND SIZE CONDITIONED FLUXES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_surface_Cluster = cell(N_Cluster,1); %fraction of surface distribution in bin
f_bar_airborne_Cluster = cell(N_Cluster,1); %fraction of mean airborne distribution in bin
f_profilebar_airborne_Cluster = cell(N_Cluster,1); %fraction of profile airborne distribution in bin
f_airborne_Cluster = cell(N_Cluster,1); %fraction of airborne sample distribution in bin
qi_Cluster = cell(N_Cluster,1); %size-selective partial flux in cluster
sigma_qi_Cluster = cell(N_Cluster,1); %size-selective partial flux uncertainty in cluster
Qi_Cluster = cell(N_Cluster,1); %size-selective flux in cluster
sigma_Qi_Cluster = cell(N_Cluster,1); %size-selective flux uncertainty in cluster
zqi_Cluster = cell(N_Cluster,1); %size-selective zq in cluster
sigma_zqi_Cluster = cell(N_Cluster,1); %size-selective zq uncertainty in cluster
Qhat_Cluster = cell(N_Cluster,1); %normalized size-selective flux in cluster
sigma_Qhat_Cluster = cell(N_Cluster,1); %normalized size-selective flux uncertainty in cluster

for i = 1:N_Cluster
    %fraction of surface distribution per bin
    f_surface_Cluster{i} = GrainSizeBinFraction(dV_bar_surface_Cluster{i}, d_surface_lower_Cluster{i}, d_surface_upper_Cluster{i}, d_bin_lower_Cluster(i,:), d_bin_upper_Cluster(i,:));
    
    %fraction of mean airborne distribution per bin
    f_bar_airborne_Cluster{i} = GrainSizeBinFraction(dV_bar_airborne_Cluster{i}, d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_lower_Cluster(i,:), d_bin_upper_Cluster(i,:));

    %fraction of profile airborne distribution per bin
    N_profile_Cluster = size(dV_profilebar_airborne_Cluster{i},1);
    f_profilebar_airborne_Cluster{i} = zeros(N_profile_Cluster,N_bins);
    for j = 1:N_profile_Cluster
        f_profilebar_airborne_Cluster{i}(j,:) = GrainSizeBinFraction(dV_profilebar_airborne_Cluster{i}(j,:), d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_lower_Cluster(i,:), d_bin_upper_Cluster(i,:));
    end
    
    %go through each sample in each profile to get size-selective partial fluxes
    f_airborne_Cluster{i} = cell(N_profile_Cluster,1); %fraction of sample per bin
    qi_Cluster{i} = cell(N_profile_Cluster,1); %size-selective partial flux in cluster
    sigma_qi_Cluster{i} = cell(N_profile_Cluster,1); %size-selective partial flux uncertainty in cluster
    for j = 1:N_profile_Cluster
    	N_z = length(z_profile_Cluster{i}{j});
        f_airborne_Cluster{i}{j} = zeros(N_z,N_bins);
        qi_Cluster{i}{j} = zeros(N_z,N_bins);
        sigma_qi_Cluster{i}{j} = zeros(N_z,N_bins);
        for k = 1:N_z
            f_airborne_Cluster{i}{j}(k,:) = GrainSizeBinFraction(dV_profile_airborne_Cluster{i}{j}(k,:), d_airborne_lower_Cluster{i}, d_airborne_upper_Cluster{i}, d_bin_lower_Cluster(i,:), d_bin_upper_Cluster(i,:)); %fraction of sample in each bin
            qi_Cluster{i}{j}(k,:) = q_profile_Cluster{i}{j}(k)*f_airborne_Cluster{i}{j}(k,:); %size-selective partial flux in each bin
            sigma_qi_Cluster{i}{j}(k,:) = sigma_q_profile_Cluster{i}{j}(k)*f_airborne_Cluster{i}{j}(k,:); %size-selective partial flux uncertainty in each bin
        end
    end
    
    %compute total fluxes
    Qi_Cluster{i} = zeros(N_profile_Cluster,N_bins); %size-selective flux in cluster
    sigma_Qi_Cluster{i} = zeros(N_profile_Cluster,N_bins); %size-selective flux uncertainty in cluster
    zqi_Cluster{i} = zeros(N_profile_Cluster,N_bins); %size-selective zq in cluster
    sigma_zqi_Cluster{i} = zeros(N_profile_Cluster,N_bins); %size-selective zq uncertainty in cluster
    Qhat_Cluster{i} = zeros(N_profile_Cluster,N_bins); %normalized size-selective flux in cluster
    sigma_Qhat_Cluster{i} = zeros(N_profile_Cluster,N_bins); %normalized size-selective flux uncertainty in cluster
    
    for j = 1:N_profile_Cluster
        for k = 1:N_bins
            [~,zqi,Qi,~,sigma_zqi,sigma_Qi] = qz_profilefit_exponential(qi_Cluster{i}{j}(:,k), z_profile_Cluster{i}{j}, sigma_qi_Cluster{i}{j}(:,k), sigma_z_profile_Cluster{i}{j}, zq_profile_Cluster{i}(j)); %perform profile fitting
            Qi_Cluster{i}(j,k) = Qi;
            sigma_Qi_Cluster{i}(j,k) = sigma_Qi;
            zqi_Cluster{i}(j,k) = zqi;
            sigma_zqi_Cluster{i}(j,k) = sigma_zqi;
        end
        Qhat_Cluster{i}(j,:) = Qi_Cluster{i}(j,:)./f_surface_Cluster{i}; %compute Qhat
        sigma_Qhat_Cluster{i}(j,:) = sigma_Qi_Cluster{i}(j,:)./f_surface_Cluster{i}; %compute sigma_Qhat
    end
end

%%%%%%%%%%%%%%%%%%%%%%
% FIT TO QHAT VS TAU %
%%%%%%%%%%%%%%%%%%%%%%

% initialize values
C_Qhat_tau = cell(N_Cluster,1);
tauth_Qhat_tau = cell(N_Cluster,1);
sigma_C_Qhat_tau = cell(N_Cluster,1);
sigma_tauth_Qhat_tau = cell(N_Cluster,1);
tau_fit_Qhat = cell(N_Cluster,1);
Qhat_fit = cell(N_Cluster,1);
sigma_Qhat_fit = cell(N_Cluster,1);

for i = 1:N_Cluster
    
    %inialize fit values
    C_Qhat_tau{i} = zeros(N_bins,1);
    tauth_Qhat_tau{i} = zeros(N_bins,1);
    sigma_C_Qhat_tau{i} = zeros(N_bins,1);
    sigma_tauth_Qhat_tau{i} = zeros(N_bins,1);
    tau_fit_Qhat{i} = cell(N_bins,1);
    Qhat_fit{i} = cell(N_bins,1);  
    sigma_Qhat_fit{i} = cell(N_bins,1);  
    
    %go through bins
    for k = 1:N_bins
        ind_fit = intersect(find(~isnan(Qhat_Cluster{i}(:,k))),find(taunorm_profile_Cluster{i}<=taunorm_max_Qhat_tau_fit)); %get indices of non-NaN values and taunorm < taunorm_max for fitting
        taunorm_max_Qhat_tau_fit = 2; %maximum normalized shear stress for fitting Qhat versus tau

        tau_fit_Qhat{i}{k} = tau_profile_Cluster{i}(ind_fit);
        Qhat = Qhat_Cluster{i}(ind_fit,k);
        sigma_Qhat = sigma_Qhat_Cluster{i}(ind_fit,k);
        
        [a, C, sigma_a, sigma_C, Qhat_fit_values, sigma_Qhat_fit_values] = linearfit(tau_fit_Qhat{i}{k},Qhat,sigma_Qhat);

        tauth = -a/C;
        sigma_tauth = -sigma_a/C;
        C_Qhat_tau{i}(k) = C;
        tauth_Qhat_tau{i}(k) = tauth;
        sigma_C_Qhat_tau{i}(k) = sigma_C;
        sigma_tauth_Qhat_tau{i}(k) = sigma_tauth;
        Qhat_fit{i}{k} = Qhat_fit_values;
        sigma_Qhat_fit{i}{k} = sigma_Qhat_fit_values;
    end
end

%%%%%%%%%%%%%%%%%%%%%
% FIT TO ZQI VS TAU %
%%%%%%%%%%%%%%%%%%%%%

% initialize values
a_zqi_tau = cell(N_Cluster,1); %fitting intercept
b_zqi_tau = cell(N_Cluster,1); %fitting slope
sigma_a_zqi_tau = cell(N_Cluster,1);
sigma_b_zqi_tau = cell(N_Cluster,1);
tau_fit_zqi = cell(N_Cluster,1);
zqi_fit = cell(N_Cluster,1);
sigma_zqi_fit = cell(N_Cluster,1);
zqi_bar_Cluster = cell(N_Cluster,1); %mean zqi
zqi_sigma_Cluster = cell(N_Cluster,1); %std dev of zqi

for i = 1:N_Cluster
    
    %inialize fit values
    a_zqi_tau{i} = zeros(N_bins,1);
    b_zqi_tau{i} = zeros(N_bins,1);
    sigma_a_zqi_tau{i} = zeros(N_bins,1);
    sigma_b_zqi_tau{i} = zeros(N_bins,1);
    tau_fit_zqi{i} = cell(N_bins,1);
    zqi_fit{i} = cell(N_bins,1);  
    sigma_zqi_fit{i} = cell(N_bins,1);
    zqi_bar_Cluster{i} = zeros(N_bins,1); %mean zqi
    zqi_sigma_Cluster{i} = zeros(N_bins,1); %std dev of zqi
    
    %go through bins
    for k = 1:N_bins
        ind_fit = find(~isnan(zqi_Cluster{i}(:,k))); %get indices of non-NaN values for fitting
        tau_fit_zqi{i}{k} = tau_profile_Cluster{i}(ind_fit);
        zqi = zqi_Cluster{i}(ind_fit,k);
        sigma_zqi = sigma_zqi_Cluster{i}(ind_fit,k);
        
        [a, b, sigma_a, sigma_b, zqi_fit_values, sigma_zqi_fit_values] = linearfit(tau_fit_zqi{i}{k},zqi,sigma_zqi);

        a_zqi_tau{i}(k) = a;
        b_zqi_tau{i}(k) = b;
        sigma_a_zqi_tau{i}(k) = sigma_a;
        sigma_b_zqi_tau{i}(k) = sigma_b;
        zqi_fit{i}{k} = zqi_fit_values;
        sigma_zqi_fit{i}{k} = sigma_zqi_fit_values;
        zqi_bar_Cluster{i}(k) = mean(zqi);
        zqi_sigma_Cluster{i}(k) = std(zqi);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT TO ZQI/D VS TAU/TAUTH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize values
a_zqinorm_taunorm = cell(N_Cluster,1); %fitting intercept
b_zqinorm_taunorm = cell(N_Cluster,1); %fitting slope
sigma_a_zqinorm_taunorm = cell(N_Cluster,1);
sigma_b_zqinorm_taunorm = cell(N_Cluster,1);
taunorm_fit_zqinorm = cell(N_Cluster,1);
zqinorm_fit = cell(N_Cluster,1);
sigma_zqinorm_fit = cell(N_Cluster,1);
zqinorm_bar_Cluster = cell(N_Cluster,1); %mean zqi/d
zqinorm_sigma_Cluster = cell(N_Cluster,1); %std dev of zqi/d

for i = 1:N_Cluster
    
    %inialize fit values
    a_zqinorm_taunorm{i} = zeros(N_bins,1);
    b_zqinorm_taunorm{i} = zeros(N_bins,1);
    sigma_a_zqinorm_taunorm{i} = zeros(N_bins,1);
    sigma_b_zqinorm_taunorm{i} = zeros(N_bins,1);
    taunorm_fit_zqinorm{i} = cell(N_bins,1);
    zqinorm_fit{i} = cell(N_bins,1);  
    sigma_zqinorm_fit{i} = cell(N_bins,1);
    zqinorm_bar_Cluster{i} = zeros(N_bins,1); %mean zqi
    zqinorm_sigma_Cluster{i} = zeros(N_bins,1); %std dev of zqi
    
    %go through bins
    for k = 1:N_bins
        ind_fit = find(~isnan(zqi_Cluster{i}(:,k))); %get indices of non-NaN values for fitting
        taunorm_fit_zqinorm{i}{k} = tau_profile_Cluster{i}(ind_fit)/tauth_Cluster(i);
        zqinorm = 1000*zqi_Cluster{i}(ind_fit,k)/d_bin_mid_Cluster(i,k);
        sigma_zqinorm = 1000*sigma_zqi_Cluster{i}(ind_fit,k)/d_bin_mid_Cluster(i,k);
        
        [a, b, sigma_a, sigma_b, zqinorm_fit_values, sigma_zqinorm_fit_values] = linearfit(taunorm_fit_zqinorm{i}{k},zqinorm,sigma_zqinorm);

        a_zqinorm_taunorm{i}(k) = a;
        b_zqinorm_taunorm{i}(k) = b;
        sigma_a_zqinorm_taunorm{i}(k) = sigma_a;
        sigma_b_zqinorm_taunorm{i}(k) = sigma_b;
        zqinorm_fit{i}{k} = zqinorm_fit_values;
        sigma_zqinorm_fit{i}{k} = sigma_zqinorm_fit_values;
        zqinorm_bar_Cluster{i}(k) = mean(zqinorm);
        zqinorm_sigma_Cluster{i}(k) = std(zqinorm);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT TO F_AIR/F_BED VS TAU %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize values
a_fratio_tau = cell(N_Cluster,1);
b_fratio_tau = cell(N_Cluster,1);
sigma_a_fratio_tau = cell(N_Cluster,1);
sigma_b_fratio_tau = cell(N_Cluster,1);
tau_fit_fratio = cell(N_Cluster,1);
fratio_fit = cell(N_Cluster,1);
fratio_bar_Cluster = cell(N_Cluster,1);
fratio_sigma_Cluster = cell(N_Cluster,1);

for i = 1:N_Cluster
    
    %inialize fit values
    a_fratio_tau{i} = zeros(N_bins,1);
    b_fratio_tau{i} = zeros(N_bins,1);
    sigma_a_fratio_tau{i} = zeros(N_bins,1);
    sigma_b_fratio_tau{i} = zeros(N_bins,1);
    tau_fit_fratio{i} = cell(N_bins,1);
    fratio_fit{i} = cell(N_bins,1);  
    fratio_bar_Cluster{i} = zeros(N_bins,1);
    fratio_sigma_Cluster{i} = zeros(N_bins,1);
    
    %go through bins
    for k = 1:N_bins
        fratio = f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k);
        ind_fit = find(~isnan(fratio)); %get indices of non-NaN values for fitting
        tau_fit = tau_profile_Cluster{i}(ind_fit);
        fratio = fratio(ind_fit);
        
        [a, b, sigma_a, sigma_b, fratio_fit_values, sigma_fratio_fit_values] = linearfit(tau_fit,fratio);

        a_fratio_tau{i}(k) = a;
        b_fratio_tau{i}(k) = b;
        sigma_a_fratio_tau{i}(k) = sigma_a;
        sigma_b_fratio_tau{i}(k) = sigma_b;
        
        [tau_fit, ind_sort] = sort(tau_fit);
        tau_fit_fratio{i}{k} = tau_fit;
        fratio_fit{i}{k} = fratio_fit_values(ind_sort);
        
        fratio_bar_Cluster{i}(k) = mean(fratio);
        fratio_sigma_Cluster{i}(k) = std(fratio);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT TO F_AIR/F_BED VS TAU/TAUTH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize values
a_fratio_taunorm = cell(N_Cluster,1);
b_fratio_taunorm = cell(N_Cluster,1);
sigma_a_fratio_taunorm = cell(N_Cluster,1);
sigma_b_fratio_taunorm = cell(N_Cluster,1);
taunorm_fit_fratio = cell(N_Cluster,1);
fratio_fit = cell(N_Cluster,1);

for i = 1:N_Cluster
    
    %inialize fit values
    a_fratio_taunorm{i} = zeros(N_bins,1);
    b_fratio_taunorm{i} = zeros(N_bins,1);
    sigma_a_fratio_taunorm{i} = zeros(N_bins,1);
    sigma_b_fratio_taunorm{i} = zeros(N_bins,1);
    taunorm_fit_fratio{i} = cell(N_bins,1);
    fratio_fit{i} = cell(N_bins,1);  
    
    %go through bins
    for k = 1:N_bins
        fratio = f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k);
        ind_fit = find(~isnan(fratio)); %get indices of non-NaN values for fitting
        taunorm_fit = tau_profile_Cluster{i}(ind_fit)/tauth_Cluster(i);
        fratio = fratio(ind_fit);
        
        [a, b, sigma_a, sigma_b, fratio_fit_values, sigma_fratio_fit_values] = linearfit(taunorm_fit,fratio);

        a_fratio_taunorm{i}(k) = a;
        b_fratio_taunorm{i}(k) = b;
        sigma_a_fratio_taunorm{i}(k) = sigma_a;
        sigma_b_fratio_taunorm{i}(k) = sigma_b;
        
        [taunorm_fit, ind_sort] = sort(taunorm_fit);
        taunorm_fit_fratio{i}{k} = taunorm_fit;
        fratio_fit{i}{k} = fratio_fit_values(ind_sort); 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP SIZE DISTRIBUTIONS BY U* BINS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_ustbar_airborne_mid_Cluster = cell(N_Cluster,1); %grain sizes for this analysis, include only grain size bins with d>d_min and d<2
dV_ustbar_airborne_Cluster = cell(N_Cluster,1);
dVdlogd_ustbar_airborne_Cluster = cell(N_Cluster,1);
for i = 1:N_Cluster
    ind_d = find(d_airborne_lower_Cluster{i}>=d_min & d_airborne_mid_Cluster{i}<=2); %include only grain size bins with d>d_min and d<2
    N_d = length(ind_d); %number of grain size bins
    d_ustbar_airborne_mid_Cluster{i} = d_airborne_mid_Cluster{i}(ind_d); %grain sizes for this analysis, include only grain size bins with d>d_min
    dV_ustbar_airborne_Cluster{i} = zeros(N_ust_bins,N_d);
    dVdlogd_ustbar_airborne_Cluster{i} = zeros(N_ust_bins,N_d);
    for j = 1:N_ust_bins
        ind_ust = find(ust_profile_Cluster{i}>=ust_bin_min(j) & ust_profile_Cluster{i}<=ust_bin_max(j));
        dV_ustbar_airborne_Cluster{i}(j,:) = mean(dV_profilebar_airborne_Cluster{i}(ind_ust,ind_d),1)/...
            sum(mean(dV_profilebar_airborne_Cluster{i}(ind_ust,ind_d),1)); %normalize by volume fraction > d_min
        dVdlogd_ustbar_airborne_Cluster{i}(j,:) = dV_ustbar_airborne_Cluster{i}(j,:)./dlogd_airborne_Cluster{i}(ind_d);
%         dVdlogd_ustbar_airborne_Cluster{i}(j,:) = mean(dVdlogd_profilebar_airborne_Cluster{i}(ind_ust,:),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP SIZE DISTRIBUTIONS BY U*/U*TH BINS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_ustnormbar_airborne_mid_Cluster = cell(N_Cluster,1); %grain sizes for this analysis, include only grain size bins with d>d_min and d<2
dV_ustnormbar_airborne_Cluster = cell(N_Cluster,1);
dVdlogd_ustnormbar_airborne_Cluster = cell(N_Cluster,1);
for i = 1:N_Cluster
    ind_d = find(d_airborne_lower_Cluster{i}>=d_min & d_airborne_mid_Cluster{i}<=2); %include only grain size bins with d>d_min and d<2
    N_d = length(ind_d); %number of grain size bins
    d_ustnormbar_airborne_mid_Cluster{i} = d_airborne_mid_Cluster{i}(ind_d); %grain sizes for this analysis, include only grain size bins with d>d_min
    dV_ustnormbar_airborne_Cluster{i} = zeros(N_ustnorm_bins,N_d);
    dVdlogd_ustnormbar_airborne_Cluster{i} = zeros(N_ustnorm_bins,N_d);
    for j = 1:N_ustnorm_bins
        ind_ustnorm = find(ustnorm_profile_Cluster{i}>=ustnorm_bin_min(j) & ustnorm_profile_Cluster{i}<=ustnorm_bin_max(j));
        dV_ustnormbar_airborne_Cluster{i}(j,:) = mean(dV_profilebar_airborne_Cluster{i}(ind_ustnorm,ind_d),1)...
            /sum(mean(dV_profilebar_airborne_Cluster{i}(ind_ustnorm,ind_d),1)); %normalize by volume fraction > d_min 
        dVdlogd_ustnormbar_airborne_Cluster{i}(j,:) = dV_ustnormbar_airborne_Cluster{i}(j,:)./dlogd_airborne_Cluster{i}(ind_d);
%         dVdlogd_ustnormbar_airborne_Cluster{i}(j,:) = mean(dVdlogd_profilebar_airborne_Cluster{i}(ind_ustnorm,ind_d),1);
    end
end

%%%%%%%%%
% PLOTS %
%%%%%%%%%

% %% Plot ust-conditioned airborne versus surface size distributions
% figure(1); clf;
% 
% %initialize subplots
% h_subplot = gobjects(N_Cluster,1);
% 
% for i = 1:N_Cluster
%         
%     %initialize subplot
%     if i < round(N_Cluster/2)
%         h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i); hold on;
%     elseif i == round(N_Cluster/2)
%         h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i:(i+1)); hold on;
%     else
%         h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i+1); hold on;
%     end
%         
%     %plot surface distribution
%     plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
%     
%     %plot airborne distribution
%     for k = 1:N_ust_bins
%         plot(d_ustbar_airborne_mid_Cluster{i},dVdlogd_ustbar_airborne_Cluster{i}(k,:),['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_ust_bin{k})
%     end
%     
%     %format plot
%     set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
%     xlim([0.06, 3]);
%     ylim([1e-4 1e1]);
%     text(0.08, 6, Label_Cluster{i},'FontSize',12);
% 
%     %label plot
%     title(ClusterNames{i});
%     if i>round(N_Cluster/2)
%         xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('Normalized volume, $$\frac{dV}{d\textrm{log}d}$$','Interpreter','Latex');
%     end
%         
%     %create legend
%     legend_items = cell(N_ust_bins+N_bins+1,1);
%     legend_items{1} = 'Surface';
%     for j = 1:N_ust_bins
%         legend_items{j+1} = ['u_{*} = ',num2str(ust_bin_min(j)),'-',num2str(ust_bin_max(j)),' m/s'];
%     end
%     for j = 1:N_bins
%         legend_items{j+N_ust_bins+1} = ['bin ',int2str(j)];
%     end
%     if i == round(N_Cluster/2)
%         h_legend = legend(legend_items,'Location','EastOutside');
%         set(h_legend,'FontSize',12);
%     end    
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[10 6],'PaperPosition',[0 0 10 6],'PaperPositionMode','Manual');
% print([folder_Plots,'GSD_ust.png'],'-dpng');
% 
% % add bins to subplots
% for i = 1:N_Cluster
%     subplot(h_subplot(i))
%     ylims = ylim;
%     for j = 1:N_bins
%         d_plot = [[1 1]*d_bin_lower_Cluster(i,j),[1 1]*d_bin_upper_Cluster(i,j)];
%         y_plot = [ylims([2,1]),ylims];
%         plot(d_plot,y_plot,'-','Color',Color_bin{j});
%     end
%     if i == round(N_Cluster/2)
%         h_legend = legend(legend_items,'Location','EastOutside');
%         set(h_legend,'FontSize',12);
%     end
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[10 6],'PaperPosition',[0 0 10 6],'PaperPositionMode','Manual');
% print([folder_Plots,'GSD_ust_bins_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% Plot ustnorm-conditioned airborne versus surface size distributions
figure(2); clf;

%initialize subplots
h_subplot = gobjects(N_Cluster,1);

for i = 1:N_Cluster

    %initialize subplot
    if N_Cluster == 4 %defined subplot sizes for four clusters
        if i == 1
            h_subplot(1) = subplot('position',[0.10 0.56 0.4 0.4]); hold on;
        elseif i == 2
            h_subplot(2) = subplot('position',[0.58 0.56 0.4 0.4]); hold on;
        elseif i == 3
            h_subplot(3) = subplot('position',[0.10 0.08 0.4 0.4]); hold on;
        else
            h_subplot(4) = subplot('position',[0.58 0.08 0.4 0.4]); hold on;
        end
    else %otherwise, automated subplot sizes
        if i < round(N_Cluster/2)
            h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i); hold on;
        elseif i == round(N_Cluster/2)
            h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i:(i+1)); hold on;
        else
            h_subplot(i) = subplot(2,round((N_Cluster+1)/2),i+1); hold on;
        end
    end
        
    %plot surface distribution
    plot(d_surface_mid_Cluster{i},dVdlogd_bar_surface_Cluster{i},'k','LineWidth',2);
    
    %plot airborne distribution
    for k = 1:N_ustnorm_bins
        plot(d_ustnormbar_airborne_mid_Cluster{i},dVdlogd_ustnormbar_airborne_Cluster{i}(k,:),['-',Marker_bin{k}],'MarkerSize',3,'Color',Color_ustnorm_bin{k})
    end
    
    %format plot
    set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlim([0.06, 2]);
    set(gca,'xtick',[0.06:0.01:0.1, 0.2:0.1:1, 2]);
    set(gca,'xticklabels',{'0.06','','','','0.1','0.2','','0.4','','0.6','','','','1','2'});
    ylim([1e-4 1e1]);
    text(0.07, 6, Label_Cluster{i},'FontSize',12);

    %label plot
    title(ClusterNames{i});
    if i>round(N_Cluster/2)
        xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex')
    end
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('Normalized volume, $$\frac{dV}{d\textrm{ln}(d)}$$','Interpreter','Latex');
    end
        
    %create legend
    legend_items = cell(N_ustnorm_bins+N_bins+1,1);
    legend_items{1} = 'Surface';
    for j = 1:N_ustnorm_bins
        if j == 1
            legend_items{j+1} = ['u_{*}/u_{*,th} \leq ',num2str(ustnorm_bin_max(j),'%10.2f')];       
        else
            legend_items{j+1} = [num2str(ustnorm_bin_min(j),'%10.2f'),' < u_{*}/u_{*,th} \leq ',num2str(ustnorm_bin_max(j),'%10.2f')];
        end
    end
    for j = 1:N_bins
        legend_items{j+N_ustnorm_bins+1} = ['bin ',int2str(j)];
    end
    if i == 1
        h_legend = legend(legend_items,'Location','SouthWest');
        set(h_legend,'FontSize',10);
    end    
%     if i == round(N_Cluster/2)
%         h_legend = legend(legend_items,'Location','EastOutside');
%         set(h_legend,'FontSize',12);
%     end    
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
print([folder_Plots,'GSD_ustnorm.png'],'-dpng');

% % add bins to subplots
% for i = 1:N_Cluster
%     subplot(h_subplot(i))
%     ylims = ylim;
%     for j = 1:N_bins
%         d_plot = [[1 1]*d_bin_lower_Cluster(i,j),[1 1]*d_bin_upper_Cluster(i,j)];
%         y_plot = [ylims([2,1]),ylims];
%         plot(d_plot,y_plot,'-','Color',Color_bin{j});
%     end
%     if i == round(N_Cluster/2)
%         h_legend = legend(legend_items,'Location','EastOutside');
%         set(h_legend,'FontSize',12);
%     end
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[10 6],'PaperPosition',[0 0 10 6],'PaperPositionMode','Manual');
% print([folder_Plots,'GSD_ustnorm_bins_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


% %% PLOT size-conditioned normalized Q VS u* - dhat bins
% figure(3); clf;
% 
% for i = 1:N_Cluster
%     
%     %initialize subplot
%     subplot(2,round(N_Cluster/2),i); hold on;
%     
%     %plot data
%     for k = 1:N_bins
%         plot(tau_profile_Cluster{i},Qhat_Cluster{i}(:,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
%     end
%     
%     %plot error bars
%     for k = 1:N_bins
%         N_tau = length(tau_profile_Cluster{i});
%         for j = 1:N_tau
%             plot(tau_profile_Cluster{i}(j)*[1 1],Qhat_Cluster{i}(j,k)*[1 1]+sigma_Qhat_Cluster{i}(j,k)*[-1 1],'Color',Color_bin{k});
%         end
%     end
%     
%     %plot fit
%     for k = 1:N_bins
%         plot([tauth_Qhat_tau{i}(k); tau_fit_Qhat{i}{k}],[0; Qhat_fit{i}{k}],'Color',Color_bin{k}); %plot fit
%     end   
%     
%     %format plot
%     xlim([0 0.5]);
%     ylims = ylim;
%     ylim([0 ylims(2)]);
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     if i>round(N_Cluster/2)
%         xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
%     end
%     ylabel('Bed fraction norm. size-spec. flux, $$\hat{Q_{i}}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
%     text(0.05, ylims(2)*0.95, Label_Cluster{i},'FontSize',9);
%     
%     %create legend
%     legend_items = cell(N_bins,1);
%     for k = 1:N_bins
%         legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
%     end
%     h_legend = legend(legend_items,'Location','EastOutside');
%     set(h_legend,'FontSize',6);
% 
%     title(ClusterNames{i});
%     set(gca,'FontSize',PlotFont);
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
% print([folder_Plots,'Qhat_tau_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');
% 
% 
% %% PLOT Qhat,i versus tau fit values
% figure(4); clf;
% 
% %plot C values and error bars
% subplot(1,2,1); hold on;
% for i = 1:N_Cluster
%     plot(dhat_bin_mid_Cluster(i,:),C_Qhat_tau{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
% end
% for i = 1:N_Cluster
%     for k = 1:N_bins
%         plot(dhat_bin_mid_Cluster(i,k)*[1 1],C_Qhat_tau{i}(k)+[-1 1]*sigma_C_Qhat_tau{i}(k),'Color',Color_Cluster{i}); %error bars
%     end
% end
% 
% %format plot
% xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
% ylabel('flux coefficient, $$C$$ (s)','Interpreter','Latex');
% legend(ClusterNames,'Location','SouthWest');
% set(gca,'Yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
% xlims = xlim;
% ylims = ylim;
% text(xlims(1)+0.1*range(xlims), ylims(1)+10^(0.95*log10(range(ylims))), '(a)','FontSize',9);
% 
% %plot tauth values
% subplot(1,2,2); hold on;
% for i = 1:N_Cluster
%     plot(dhat_bin_mid_Cluster(i,:),tauth_Qhat_tau{i},Marker_Cluster{i},'Color',Color_Cluster{i});
% end
% for i = 1:N_Cluster
%     for k = 1:N_bins
%         plot(dhat_bin_mid_Cluster(i,k)*[1 1],tauth_Qhat_tau{i}(k)+[-1 1]*sigma_tauth_Qhat_tau{i}(k),'Color',Color_Cluster{i}); %error bars
%     end
% end
% 
% %format plot
% xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
% ylabel('threshold, $$\tau_{th}$$ (Pa)','Interpreter','Latex');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% xlims = xlim;
% ylims = ylim;
% text(xlims(1)+0.1*range(xlims), ylims(1)+0.99*range(ylims), '(b)','FontSize',9);
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'C_tauth_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');
% 
% %plot only positive tauth values
% ylim([0 ylims(2)]);
% print([folder_Plots,'C_tauth_dhat_positive_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');
% 
% 
% %% plot f_airborne / f_surface versus d - binned
% figure(5); clf;
% 
% for i = 1:N_Cluster
%     
%     %initialize subplot
%     subplot(2,round(N_Cluster/2),i); hold on;
%     
%     %plot data
%     for k = 1:N_bins
%         plot(tau_profile_Cluster{i},f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
%     end
%     
%     %plot fit
%     for k = 1:N_bins
%         plot(tau_fit_fratio{i}{k},fratio_fit{i}{k},'Color',Color_bin{k}); %plot fit
%     end   
%     
%     %format plot
%     set(gca,'YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
%     xlim([0 0.5]);
%     ylims = ylim;
%     if i>round(N_Cluster/2)
%         xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
%     end
%     ylabel('Airborne / surface fraction, $$f_{air}/f_{sfc}$$','Interpreter','Latex');
%     text(0.05, ylims(2)*0.95, Label_Cluster{i},'FontSize',9);
%     
%     %create legend
%     legend_items = cell(N_bins,1);
%     for k = 1:N_bins
%         legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
%     end
%     h_legend = legend(legend_items,'Location','EastOutside');
%     set(h_legend,'FontSize',6);
% 
%     title(ClusterNames{i});
%     set(gca,'FontSize',PlotFont);
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
% print([folder_Plots,'fratio_tau_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');
% 
% 
% %% plot f_airborne / f_surface versus d
% figure(6); clf;
% 
% for i = 1:N_Cluster
%     
%     %initialize subplot
%     subplot(2,round(N_Cluster/2),i); hold on;
%     
%     %plot data
%     for k = 1:N_bins
%         plot(tau_profile_Cluster{i},f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
%     end
%     
%     %plot fit
%     for k = 1:N_bins
%         plot(tau_fit_fratio{i}{k},fratio_fit{i}{k},'Color',Color_bin{k}); %plot fit
%     end   
%     
%     %format plot
%     set(gca,'YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
%     xlim([0 0.5]);
%     ylims = ylim;
%     xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
%     if i==1
%         ylabel('Airborne / surface fraction, $$f_{air}/f_{bed}$$','Interpreter','Latex');
%     end
%     text(0.05, ylims(2)*0.95, Label_Cluster{i},'FontSize',9);
%     
%     %create legend
%     legend_items = cell(N_bins,1);
%     for k = 1:N_bins
%         legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
%     end
%     h_legend = legend(legend_items,'Location','EastOutside');
%     set(h_legend,'FontSize',6);
% 
%     title(ClusterNames{i});
%     set(gca,'FontSize',PlotFont);
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
% print([folder_Plots,'fratio_tau_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');
% 
% 
% %% PLOT mean fratio / fsurface - unbinned
% figure(7); clf; hold on;
% 
% %plot fratio values
% for i = 1:N_Cluster
%     if strcmp(dref_type,'d50') 
%         plot(d_surface_mid_Cluster{i}./d50_bar_surface_Cluster(i),dV_bar_airborne_Cluster{i}./dV_bar_surface_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
%     elseif strcmp(dref_type,'dmodal')
%         plot(d_surface_mid_Cluster{i}./dmodal_bar_surface_Cluster(i),dV_bar_airborne_Cluster{i}./dV_bar_surface_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
%     end
% end
% 
% %format plot
% if strcmp(dref_type,'d50')==1
%     xlabel('Normalized grain size, $$d / d_{50}$$','Interpreter','Latex');
% elseif strcmp(dref_type,'dmodal')==1
%     xlabel('Normalized grain size, $$d / d_{modal}$$','Interpreter','Latex');
% end
% ylabel('Mean airborne / bed fraction, $$\langle f_{air}/f_{bed} \rangle$$','Interpreter','Latex');
% legend(ClusterNames,'Location','SouthWest');
% set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
% print([folder_Plots,'fratiobar_tau_unbinned_',dref_type,'.png'],'-dpng');
% 
% 
% %% PLOT mean fratio / fsurface - unbinned - limited d 
% figure(8); clf; hold on;
% 
% %plot fratio values
% for i = 1:N_Cluster
%     ind_plot = find(d_surface_lower_Cluster{i} >= d_min & d_surface_upper_Cluster{i} <= d_max); %get only values in range of trustworthy d
%     plot(d_surface_lower_Cluster{i}(ind_plot),dV_bar_airborne_Cluster{i}(ind_plot)./dV_bar_surface_Cluster{i}(ind_plot),Marker_Cluster{i},'Color',Color_Cluster{i}); %values
% end
% 
% %format plot
% xlabel('Grain size, $$d$$ (mm)','Interpreter','Latex');
% ylabel('Mean airborne / bed fraction, $$\langle f_{air}/f_{bed} \rangle$$','Interpreter','Latex');
% legend(ClusterNames,'Location','SouthWest');
% set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
% xlim([d_min d_max]);
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
% print([folder_Plots,'fratiobar_d_unbinned_limited_',dref_type,'.png'],'-dpng');


%% PLOT mean fratio / fsurface - unbinned - limited dhat
figure(9); clf; hold on;

%plot fratio values
for i = 1:N_Cluster
    ind_plot = find(d_surface_mid_Cluster{i} >= d_min & d_surface_mid_Cluster{i} <= d_max); %get only values in range of trustworthy d
    
    if strcmp(dref_type,'d50') 
        plot(d_surface_mid_Cluster{i}(ind_plot)./d50_bar_surface_Cluster(i),dV_bar_airborne_Cluster{i}(ind_plot)./dV_bar_surface_Cluster{i}(ind_plot),Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    elseif strcmp(dref_type,'dmodal')
        plot(d_surface_mid_Cluster{i}(ind_plot)./dmodal_bar_surface_Cluster(i),dV_bar_airborne_Cluster{i}(ind_plot)./dV_bar_surface_Cluster{i}(ind_plot),Marker_Cluster{i},'Color',Color_Cluster{i}); %values
    end
end

%format plot
if strcmp(dref_type,'d50')==1
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
elseif strcmp(dref_type,'dmodal')==1
    xlabel('Normalized grain size, $$d / d_{modal,bed}$$','Interpreter','Latex');
end
ylabel('Mean airborne / bed fraction, $$\langle f_{air} \rangle/f_{bed}$$','Interpreter','Latex');
legend(ClusterNames,'Location','NorthEast');
set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
ylim([1e-4 1e1]);

%inset plot with dimensional values
axes('Position',[.25 .22 .35 .38]); hold on;
%plot fratio values
for i = 1:N_Cluster
    ind_plot = find(d_surface_mid_Cluster{i} >= d_min & d_surface_mid_Cluster{i} <= d_max); %get only values in range of trustworthy d
    plot(d_surface_mid_Cluster{i}(ind_plot),dV_bar_airborne_Cluster{i}(ind_plot)./dV_bar_surface_Cluster{i}(ind_plot),Marker_Cluster{i},'Color',Color_Cluster{i}); %values
end
xlabel('Grain size, $$d$$ (mm)','Interpreter','Latex');
ylabel('$$\langle f_{air} \rangle /f_{bed}$$','Interpreter','Latex');
set(gca,'XScale','log','YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
xlim([d_min d_max]);
ylim([1e-4 1e1]);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 5],'PaperPosition',[0 0 7 5],'PaperPositionMode','Manual');
print([folder_Plots,'fratiobar_dhat_unbinned_limited_',dref_type,'.png'],'-dpng');


% %% PLOT mean fratio / fsurface - binned
% figure(10); clf; hold on;
% 
% %plot fratio values
% for i = 1:N_Cluster
%     plot(d_bin_mid_Cluster(i,:),f_bar_airborne_Cluster{i}./f_surface_Cluster{i},Marker_Cluster{i},'Color',Color_Cluster{i}); %values
% end
% 
% %format plot
% xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
% ylabel('Mean airborne / surface fraction, $$\langle f_{air}/f_{sfc} \rangle$$','Interpreter','Latex');
% legend(ClusterNames,'Location','SouthWest');
% set(gca,'YScale','log','XMinorTick','On','YMinorTick','On','Box','On');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 4],'PaperPosition',[0 0 6 4],'PaperPositionMode','Manual');
% print([folder_Plots,'fratiobar_tau_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');
% 
% 
%% PLOT normalized size-conditioned zq VS u*
figure(11); clf;

for i = 1:N_Cluster
    
    %initialize subplot
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot data
    for k = 1:N_bins
        plot(tau_profile_Cluster{i},1000*zqi_Cluster{i}(:,k)./d_bin_mid_Cluster(i,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
    end
    
    %plot error bars
    for k = 1:N_bins
        N_tau = length(tau_profile_Cluster{i});
        for j = 1:N_tau
            plot(tau_profile_Cluster{i}(j)*[1 1],1000*(zqi_Cluster{i}(j,k)*[1 1]+sigma_zqi_Cluster{i}(j,k)*[-1 1])./d_bin_mid_Cluster(i,k),'Color',Color_bin{k});
        end
    end
    
    %plot fit
    for k = 1:N_bins
        if ~isempty(tau_fit_zqi{i}{k})
            plot(tau_fit_zqi{i}{k},1000*zqi_fit{i}{k}./d_bin_mid_Cluster(i,k),'Color',Color_bin{k}); %plot fit
        end
    end   
    
    %format plot
    xlim([0 0.5]);
    ylim([0 400]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    if i>round(N_Cluster/2)
        xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
    end
    ylabel('Norm. size-sel. salt. ht., $$z_{q,i}/d_{i}$$','Interpreter','Latex');
    text(0.02, 380, Label_Cluster{i},'FontSize',9);
    
    %create legend
    legend_items = cell(N_bins,1);
    for k = 1:N_bins
        legend_items{k} = ['d/d_{ref} = ',num2str(dhat_bin_mid_Cluster(i,k),2)];
    end
    h_legend = legend(legend_items,'Location','EastOutside');
    set(h_legend,'FontSize',6);

    title(ClusterNames{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[12 6],'PaperPosition',[0 0 12 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_tau_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT mean size-conditioned zq VS dhat
figure(12); clf; hold on;

%plot data
for i = 1:N_Cluster
    plot(dhat_bin_mid_Cluster(i,:),zqinorm_bar_Cluster{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],zqinorm_bar_Cluster{i}(k)*[1 1]+zqinorm_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end
   
%format plot
xlim([0.28 2]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
set(gca,'XTick',[0.3:0.1:1,2],'XTickLabel',{'0.3','0.4','0.5','','0.7','','','1','2'});
if strcmp(dref_type,'dmodal')
    xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
else
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
end
ylabel('Mean normalized size-selective saltation height, $$\langle z_{q,i} \rangle/d_{i}$$','Interpreter','Latex');

%create legend
h_legend = legend(ClusterNames,'Location','NorthEast');
set(h_legend,'FontSize',10);

%inset plot with dimensional values
axes('Position',[.23 .2 .35 .3]); hold on;

%plot data
for i = 1:N_Cluster
    plot(d_bin_mid_Cluster(i,:),zqi_bar_Cluster{i}',[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(d_bin_mid_Cluster(i,k)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end

%format plot
%ylim([-0.1 0.15]);
xlim([0.12 1.0]);
xlabel('grain size, $$d$$ (mm)','Interpreter','Latex');
ylabel('Mean size-sel. salt. ht., $$\langle z_{q,i} \rangle$$ (m)','Interpreter','Latex');
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_dhat_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


%% PLOT mean size-conditioned zq VS d
figure(13); clf; hold on;

%plot data
for i = 1:N_Cluster
    plot(d_bin_mid_Cluster(i,:),zqi_bar_Cluster{i}',[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
        plot(d_bin_mid_Cluster(i,k)*[1 1],zqi_bar_Cluster{i}(k)*[1 1]+zqi_sigma_Cluster{i}(k)*[-1 1],'Color',Color_Cluster{i});
    end
end

%format plot
%ylim([-0.1 0.15]);
xlim([0.12 1.0]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
xlabel('Grain size, $$d$$ (mm)','Interpreter','Latex');
ylabel('Mean size-selective saltation height, $$\langle z_{q,i} \rangle$$ (m)','Interpreter','Latex');

%create legend
h_legend = legend(ClusterNames,'Location','SouthWest');
set(h_legend,'FontSize',10);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
print([folder_Plots,'zq_d_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');



%% PLOT slope of fair/fbed VS tau/tauth and slope of zq/d VS tau/tauth
figure(14); clf;

% PLOT slope of fair/fbed
subplot(2,1,1); hold on;

%plot data
for i = 1:N_Cluster
%    plot(dhat_bin_mid_Cluster(i,:),b_fratio_taunorm{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
    plot(dhat_bin_mid_Cluster(i,:),b_fratio_taunorm{i}./fratio_bar_Cluster{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
%        plot(dhat_bin_mid_Cluster(i,k)*[1 1],b_fratio_taunorm{i}(k)*[1 1]+sigma_b_fratio_taunorm{i}(k)*[-1 1],'Color',Color_Cluster{i});
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],(b_fratio_taunorm{i}(k)*[1 1]+sigma_b_fratio_taunorm{i}(k)*[-1 1])./fratio_bar_Cluster{i}(k),'Color',Color_Cluster{i});
    end
end
   
%format plot
xlim([0.4 2]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
set(gca,'XTick',[0.4:0.1:1,2],'XTickLabel',{'0.4','0.5','','0.7','','','1','2'});
if strcmp(dref_type,'dmodal')
    xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
else
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
end
%ylabel('Airborne fraction vs. shear stress trend, $$\frac{\textrm{d}f_{air}/f_{bed}}{\textrm{d}\tau/\tau_{th}}$$','Interpreter','Latex');
ylabel('Norm air/bed frac vs stress trend, $$(\frac{\textrm{d}f_{air}/f_{bed}}{\textrm{d}\tau/\tau_{th}})/\langle f_{air}/f_{bed} \rangle$$','Interpreter','Latex');

% PLOT slope of zq/d
subplot(2,1,2); hold on;

%plot data
for i = 1:N_Cluster
%    plot(dhat_bin_mid_Cluster(i,:),b_zqinorm_taunorm{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
    plot(dhat_bin_mid_Cluster(i,:),b_zqinorm_taunorm{i}./zqinorm_bar_Cluster{i},[Marker_Cluster{i},'-'],'Color',Color_Cluster{i});
end

%plot error bars
for i = 1:N_Cluster
    for k = 1:N_bins
%        plot(dhat_bin_mid_Cluster(i,k)*[1 1],b_zqinorm_taunorm{i}(k)*[1 1]+sigma_b_zqinorm_taunorm{i}(k)*[-1 1],'Color',Color_Cluster{i});
        plot(dhat_bin_mid_Cluster(i,k)*[1 1],(b_zqinorm_taunorm{i}(k)*[1 1]+sigma_b_zqinorm_taunorm{i}(k)*[-1 1])./zqinorm_bar_Cluster{i}(k),'Color',Color_Cluster{i});
    end
end
   
%format plot
xlim([0.4 2]);
set(gca,'XScale','log','XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
set(gca,'XTick',[0.4:0.1:1,2],'XTickLabel',{'0.4','0.5','','0.7','','','1','2'});
if strcmp(dref_type,'dmodal')
    xlabel('Normalized grain size, $$d / d_{ref}$$','Interpreter','Latex');
else
    xlabel('Normalized grain size, $$d / d_{50,bed}$$','Interpreter','Latex');
end
%ylabel('Saltation height vs. shear stress trend, $$\frac{\textrm{d}z_{q,i}/d_{i}}{\textrm{d}\tau/\tau_{th}}$$','Interpreter','Latex');
ylabel('Norm salt ht vs stress trend, $$(\frac{\textrm{d}z_{q,i}/d_{i}}{\textrm{d}\tau/\tau_{th}})/\langle z_{q,i} \rangle$$','Interpreter','Latex');

%create legend
h_legend = legend(ClusterNames,'Location','NorthWest');
set(h_legend,'FontSize',10);

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 8],'PaperPosition',[0 0 6 8],'PaperPositionMode','Manual');
print([folder_Plots,'fratio_zqinorm_taunorm_d_',binning_type,'_',int2str(N_bins),'_',dref_type,'.png'],'-dpng');


% %% plot variation in reference grain sizes with saltation flux
% figure(15); clf;
% for i = 1:N_Cluster
%     subplot(2,round(N_Cluster/2),i); hold on;
%     h90_air = plot(Q_profile_Cluster{i},d90_profilebar_airborne_Cluster{i},'^');
%     h50_air = plot(Q_profile_Cluster{i},d50_profilebar_airborne_Cluster{i},'o');
%     h10_air = plot(Q_profile_Cluster{i},d10_profilebar_airborne_Cluster{i},'v');
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     h90_sfc = plot([0 60],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_sfc = plot([0 60],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_sfc = plot([0 60],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     xlim([0 60]);
%     ylim([0 0.9]);
%     if i>round(N_Cluster/2)
%         xlabel('saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air], 'd_{90,air}','d_{50,air}','d_{10,air}');
%     elseif i == N_Cluster
%         legend([h90_sfc, h50_sfc, h10_sfc], 'd_{90,surface}','d_{50,surface}','d_{10,surface}');
%     end
%     title(ClusterNames{i});
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'IndexGrainSize_SaltationFlux.png'],'-dpng');
% 
% 
%% plot variation in reference grain sizes with shear velocity
figure(16); clf;
for i = 1:N_Cluster
    subplot(2,round(N_Cluster/2),i); hold on;
    h90_air = plot(ust_profile_Cluster{i},d90_profilebar_airborne_Cluster{i},'^');
    h50_air = plot(ust_profile_Cluster{i},d50_profilebar_airborne_Cluster{i},'o');
    h10_air = plot(ust_profile_Cluster{i},d10_profilebar_airborne_Cluster{i},'v');
    hbar_air = plot(ust_profile_Cluster{i},dbar_profilebar_airborne_Cluster{i},'s');
    c90 = get(h90_air,'Color');
    c50 = get(h50_air,'Color');
    c10 = get(h10_air,'Color');
    cbar = get(hbar_air,'Color');
    h90_bed = plot([0.2 0.6],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
    h50_bed = plot([0.2 0.6],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
    h10_bed = plot([0.2 0.6],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
    hbar_bed = plot([0.2 0.6],dbar_bar_surface_Cluster(i)*[1 1],'Color',cbar);
    xlim([0.2 0.6]);
    ylim([0 0.9]);
    if i>round(N_Cluster/2)
        xlabel('shear velocity, $$u_{*}$$ (m s$$^{-1}$$)','Interpreter','Latex')
    end
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
    end
    if i == N_Cluster - 1
        legend([h90_air, h50_air, h10_air,hbar_air], 'd_{90,air}','d_{50,air}','d_{10,air}','d_{bar,air}');
    elseif i == N_Cluster
        legend([h90_bed, h50_bed, h10_bed,hbar_bed], 'd_{90,bed}','d_{50,bed}','d_{10,bed}','d_{bar,bed}');
    end
    title(ClusterNames{i});
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 7],'PaperPosition',[0 0 6 7],'PaperPositionMode','Manual');
print([folder_Plots,'ReferenceGrainSize_ShearVelocity.png'],'-dpng');

%% plot variation in reference grain sizes with shear stress
figure(17); clf;
for i = 1:N_Cluster
    subplot(2,round(N_Cluster/2),i); hold on;
    
    %plot airborne sizes
    h90_air = plot(tau_profile_Cluster{i},d90_profilebar_airborne_Cluster{i},'^');
    h50_air = plot(tau_profile_Cluster{i},d50_profilebar_airborne_Cluster{i},'o');
    h10_air = plot(tau_profile_Cluster{i},d10_profilebar_airborne_Cluster{i},'v');
    hbar_air = plot(tau_profile_Cluster{i},dbar_profilebar_airborne_Cluster{i},'s');
    
    %get info about airborne plots
    c90 = get(h90_air,'Color');
    c50 = get(h50_air,'Color');
    c10 = get(h10_air,'Color');
    cbar = get(hbar_air,'Color');
    
    %plot airborne uncertainties
    for j = 1:length(tau_profile_Cluster{i})
        plot(tau_profile_Cluster{i}(j)*[1 1],d90_profilebar_airborne_Cluster{i}(j)*[1 1]+d90_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',c90);
        plot(tau_profile_Cluster{i}(j)*[1 1],d50_profilebar_airborne_Cluster{i}(j)*[1 1]+d50_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',c50);
        plot(tau_profile_Cluster{i}(j)*[1 1],d10_profilebar_airborne_Cluster{i}(j)*[1 1]+d10_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',c10);
        plot(tau_profile_Cluster{i}(j)*[1 1],dbar_profilebar_airborne_Cluster{i}(j)*[1 1]+dbar_profilesigma_airborne_Cluster{i}(j)*[-1 1],'Color',cbar);
    end
        
    %make surface plots
    h90_bed = plot([0 0.45],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
    h50_bed = plot([0 0.45],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
    h10_bed = plot([0 0.45],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
    hbar_bed = plot([0 0.45],dbar_bar_surface_Cluster(i)*[1 1],'Color',cbar);       
    xlim([0 0.45]);
    ylim([0 0.9]);
    if i>round(N_Cluster/2)
        xlabel('shear stress, $$\tau$$ (Pa)','Interpreter','Latex')
    end
    if mod(i,round(N_Cluster/2)) == 1
        ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
    end
    if i == N_Cluster - 1
        legend([h90_air, h50_air, h10_air,hbar_air], 'd_{90,air}','d_{50,air}','d_{10,air}','d_{bar,air}');
    elseif i == N_Cluster
        legend([h90_bed, h50_bed, h10_bed,hbar_bed], 'd_{90,bed}','d_{50,bed}','d_{10,bed}','d_{bar,bed}');
    end
    title(ClusterNames{i});
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gcf,'PaperUnits','inches','PaperSize',[6 7],'PaperPosition',[0 0 6 7],'PaperPositionMode','Manual');
print([folder_Plots,'ReferenceGrainSize_ShearStress.png'],'-dpng');

% 
% %% plot variation in reference grain sizes with height
% figure(17); clf;
% for i = 1:N_Cluster
%     subplot(2,round(N_Cluster/2),i); hold on;
%     h90_air = plot(z_airborne_Cluster{i},d90_airborne_Cluster{i},'^');
%     h50_air = plot(z_airborne_Cluster{i},d50_airborne_Cluster{i},'o');
%     h10_air = plot(z_airborne_Cluster{i},d10_airborne_Cluster{i},'v');
%     c90 = get(h90_air,'Color');
%     c50 = get(h50_air,'Color');
%     c10 = get(h10_air,'Color');
%     h90_sfc = plot([0 0.6],d90_bar_surface_Cluster(i)*[1 1],'Color',c90);
%     h50_sfc = plot([0 0.6],d50_bar_surface_Cluster(i)*[1 1],'Color',c50);
%     h10_sfc = plot([0 0.6],d10_bar_surface_Cluster(i)*[1 1],'Color',c10);
%     xlim([0 0.6]);
%     ylim([0 0.9]);
%     if i>round(N_Cluster/2)
%         xlabel('height above surface, $$z$$ (m)','Interpreter','Latex')
%     end
%     if mod(i,round(N_Cluster/2)) == 1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     end
%     if i == N_Cluster - 1
%         legend([h90_air, h50_air, h10_air], 'd_{90,air}','d_{50,air}','d_{10,air}');
%     elseif i == N_Cluster
%         legend([h90_sfc, h50_sfc, h10_sfc], 'd_{90,surface}','d_{50,surface}','d_{10,surface}');
%     end
%     title(ClusterNames{i});
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[6 5],'PaperPosition',[0 0 6 5],'PaperPositionMode','Manual');
% print([folder_Plots,'IndexGrainSize_Height.png'],'-dpng');