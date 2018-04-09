%% initialize
clearvars;

%%
%%%%%%%%%%%%%%%%
% DATA SOURCES %
%%%%%%%%%%%%%%%%

%% information about where to load and save data
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
folder_SaltationData = '../../AnalysisData/Windowing/'; %folder for saltation flux data
folder_AnalysisData = '../../AnalysisData/SizeSelective/'; %folder for saving analysis data
GrainSizeData_Path = strcat(folder_GrainSizeData,'GrainSizeData'); %path for loading mean grain size data
SaltationFluxData_Path = strcat(folder_SaltationData,'DataWindowCalcs_30min_Restricted'); %path for loading saltation data
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(GrainSizeData_Path); %grain size data
load(SaltationFluxData_Path); %saltation data
addpath(folder_Functions); %point MATLAB to location of functions

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
tauth_Cluster = [0.111, 0.110, 0.083, 0.083, 0.091, 0.091]'; %shear stress threshold by Cluster based on thresholds paper
sigma_tauth_Cluster = [0.002, 0.002, 0.001, 0.001, 0.001, 0.001]'; %shear stress uncertainty by Cluster based on thresholds paper
rho_a = [1.16, 1.22, 1.22, 1.22, 1.22, 1.22]; %air density by Cluster
ustth_Cluster = sqrt(tauth_Cluster./rho_a); %threshold shear velocity by Cluster
sigma_ustth_Cluster = sigma_tauth_Cluster./(2*rho_a.*ustth_Cluster); %uncertainty in threshold shear velocity by Cluster

N_Sites = length(SiteNames);
Cluster_StartDate = {...
    datetime(2014,11,13),...
    datetime(2015,3,23),...
    [datetime(2015,5,15); datetime(2015,5,15); datetime(2015,5,23); datetime(2015,5,23)]
    };
Cluster_EndDate = {...
    datetime(2014,11,20),...
    datetime(2015,3,24),...
    [datetime(2015,5,22); datetime(2015,5,22); datetime(2015,6,4); datetime(2015,6,4)]
    };
Cluster_ymin = {...
    -Inf,...
    -Inf,...
    [0; -Inf; 0; -Inf]
    };
Cluster_ymax = {...
    Inf,...
    Inf,...
    [Inf; 0; Inf; 0]
    };
Cluster_Location = {...
    {{'BSNE_A','BSNE_B','Tower_A','Tower_B','U1','Upwind_X0m','Upwind_X5m','Upwind_X10m','Upwind_X15m','Upwind_X20m'}},...
    {{'A','B','C','Wenglor','A','B','C','Wenglor'}},...
    {{'A','B','A1','A2','A3','A4'},{'C','D','B1','B2_B3','B4'},{'A','B','A1','A2','A3','A4'},{'C','D','B1','B2_B3','B4'}},...
    };
N_Cluster_Site = cellfun(@length,Cluster_StartDate); %number of surface grain size clusters for each site (Jeri, RG, Oceano);
N_Cluster = sum(N_Cluster_Site); %total number of clusters

%generate names of clusters
ClusterNames = {'Jericoacoara','Rancho Guadalupe','Oceano - N1','Oceano - S1','Oceano - N2','Oceano - S2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR PROFILE AVERAGING %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

znorm_maxlower = 2.5; %max lower z/zq
znorm_minupper = 4.5; %min upper z/zq
N_min_profile = 3; %minimum number of values in profile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR BINNING BY GRAIN SIZE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% partial flux grain-size fixed bin values
% binning_type = 'fixed';
% N_bins = 6; %number of bins
% d_bin_min = d_min; %minimum bin (mm)
% d_bin_max = d_max; %maximum bin (mm)
% d_bin_edges = logspace(log10(d_bin_min),log10(d_bin_max),N_bins+1);
% d_bin_lower = d_bin_edges(1:N_bins); %lower edge of bin (mm)
% d_bin_upper = d_bin_edges(2:N_bins+1); %upper edge of bin (mm)
% d_bin_mid = geomean([d_bin_lower; d_bin_upper],1); %midpoint of bin (mm)

%% define bins in terms of values with respect to dref (i.e., d/dref)
binning_type = 'dhat';
N_bins = 6; %number of bins
dhat_bin_min = 0.37; %bottom of lower bin -- associated with 0.13 mm minimum grain diameter for BSNE collection efficiency
dhat_bin_edges = logspace(log10(dhat_bin_min),log10(1/dhat_bin_min),N_bins+2); %edges of bins -- log-spaced
dhat_bin_lower = dhat_bin_edges(1:N_bins); %d/dref values for size bins - lower limit
dhat_bin_upper = dhat_bin_edges(2:N_bins+1); %d/dref values for size bins - upper limit
dhat_bin_mid = geomean([dhat_bin_lower; dhat_bin_upper],1); %d/dref values for size bins - midpoint
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
% ustnorm_min = 1.0;
% ustnorm_max = 2.0;
% N_ustnorm_bins = 4;
% ustnorm_bin_edges = [0,logspace(log10(ustnorm_min),log10(ustnorm_max),N_ustnorm_bins)];
% ustnorm_bin_min = ustnorm_bin_edges(1:N_ustnorm_bins); %u*norm for bottom of normalized shear velocity bins
% ustnorm_bin_max = ustnorm_bin_edges(2:N_ustnorm_bins+1); %u*norm for top of normalized shear velocity bins
ustnorm_bin_min = [0,1,1.25,1.6]; %set limits manually
ustnorm_bin_max = [1,1.25,1.6,2]; %set limits manually
N_ustnorm_bins = length(ustnorm_bin_min); %get number of bins

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR BINNING BY NORMALIZED SHEAR STRESS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
taunorm_min_bins = [0,1,1.5,2.2,3]; %set limits manually
taunorm_max_bins = [1,1.5,2.2,3,4]; %set limits manually
N_taunorm_bins = length(taunorm_min_bins); %get number of bins


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
znorm_airborne_Cluster = cell(N_Cluster,1); %normalized height associated with individual sample

%information for flux profiles
dV_profile_airborne_Cluster = cell(N_Cluster,1); %volume fraction of airborne individual sample
dVdlogd_profile_airborne_Cluster = cell(N_Cluster,1); %differential volume fraction of airborne individual sample
d10_profile_airborne_Cluster = cell(N_Cluster,1); %d10 of airborne sample profiles
d50_profile_airborne_Cluster = cell(N_Cluster,1); %d50 of airborne sample profiles
d90_profile_airborne_Cluster = cell(N_Cluster,1); %d90 of airborne sample profiles
dbar_profile_airborne_Cluster = cell(N_Cluster,1); %mean d of airborne sample profiles
q_profile_Cluster = cell(N_Cluster,1); %partial flux profiles
sigma_q_profile_Cluster = cell(N_Cluster,1); %partial flux uncertainty profiles
z_profile_Cluster = cell(N_Cluster,1); %height profiles
sigma_z_profile_Cluster = cell(N_Cluster,1); %height uncertainty profiles
Q_profile_Cluster = cell(N_Cluster,1); %total flux associated with sample profile
sigma_Q_profile_Cluster = cell(N_Cluster,1); %total flux associated with sample profile
zq_profile_Cluster = cell(N_Cluster,1); %e-folding height associated with sample profile
znorm_profile_Cluster = cell(N_Cluster,1); %normalized e-folding height associated with sample profile
tau_profile_Cluster = cell(N_Cluster,1); %shear stress associated with sample profile
ust_profile_Cluster = cell(N_Cluster,1); %shear velocity associated with sample profile
taunorm_profile_Cluster = cell(N_Cluster,1); %normalized shear stress (tau/tauth) associated with sample profile
ustnorm_profile_Cluster = cell(N_Cluster,1); %normalized shear velocity (ust/ustth) associated with sample profile
Date_profile_Cluster = cell(N_Cluster,1); %Dates
Time_profile_Cluster = cell(N_Cluster,1); %Times

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

zq_bar_Cluster = zeros(N_Cluster,1); %mean zq value for each Cluster
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

ind_usable_profile_Cluster = cell(N_Cluster,1); %set to 1 if profile is usable for analysis

%% go through sites
for i = 1:N_Sites
    
    %% surface samples
    % load relevant information about surface samples
    GrainSize_surface = GrainSizeData_all{i}.GrainSize_Surface; %surface sample array
    Date_surface = [GrainSize_surface.Date]; %dates of surface samples
    Location_surface = {GrainSize_surface.Location}; %locations of surface samples
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
        dV_airborne_samples(j,:) = dV_airborne_samples(j,:)/sum(dV_airborne_samples(j,:)); %normalize so that sum equals 1
        dVdlogd_airborne_samples(j,:) = dV_airborne_samples(j,:)./dlogd_airborne; %by normalized volume
    end
    
    %get mean d of each airborne size distribution
    dbar_airborne_samples = zeros(N_airborne_samples,1); %initialize list of dbar
    for j = 1:N_airborne_samples %go through each airborne sample
        dbar_airborne_samples(j,:) = exp(sum(dV_airborne_samples(j,:).*log(d_airborne_mid))); %compute weighted geometric mean of particle size bins
    end
    
    %get spanwise location and height of each airborne sample
    y_airborne_samples = zeros(N_airborne_samples,1)*NaN; %spanwise location of each airborne sample
    z_airborne_samples = zeros(N_airborne_samples,1)*NaN; %height of each airborne sample
    for j = 1:N_airborne_samples %go through each airborne sample
        ind_BSNE = find([Flux_BSNE.StartTime] >= StartTime_airborne_samples(j) & [Flux_BSNE.StartTime] <= EndTime_airborne_samples(j), 1);
        ind_z = find(strcmp(Flux_BSNE(ind_BSNE).name,Name_airborne_samples{j}));
        if ~isempty(ind_z)
            y_airborne_samples(j) = Flux_BSNE(ind_BSNE).y(ind_z);
            z_airborne_samples(j) = Flux_BSNE(ind_BSNE).z.z(ind_z);
        end
    end
        
    %% combine airborne grain size distributions for each BSNE collection interval  
    
    %get dates for profiles
    Date_profile_airborne = [Flux_BSNE.Date];
    Time_profile_airborne = [Flux_BSNE.EndTime];
    
    %initialize profile values
    dV_profile_airborne = cell(N_Flux,1); %volume fraction of airborne individual sample
    dVdlogd_profile_airborne = cell(N_Flux,1); %differential volume fraction of airborne individual sample
    d10_profile_airborne = cell(N_Flux,1); %d10 of airborne individual sample
    d50_profile_airborne = cell(N_Flux,1); %d50 of airborne individual sample
    d90_profile_airborne = cell(N_Flux,1); %d90 of airborne individual sample
    dbar_profile_airborne = cell(N_Flux,1); %dbar of airborne individual sample
    q_profile_airborne = cell(N_Flux,1); %partial flux profiles
    sigma_q_profile_airborne = cell(N_Flux,1); %partial flux uncertainty profiles
    y_profile_airborne = cell(N_Flux,1); %spanwise location profiles
    z_profile_airborne = cell(N_Flux,1); %height profiles
    sigma_z_profile_airborne = cell(N_Flux,1); %height uncertainties profiles
    Q_profile_airborne = zeros(N_Flux,1); %total flux associated with sample profile
    Q_1_profile_airborne = zeros(N_Flux,1)*NaN; %total flux associated with +y side of sample profile
    Q_2_profile_airborne = zeros(N_Flux,1)*NaN; %total flux associated with -y side of sample profile
    sigma_Q_profile_airborne = zeros(N_Flux,1); %total flux uncertainty associated with sample profile
    sigma_Q_1_profile_airborne = zeros(N_Flux,1)*NaN; %total flux uncertainty associated with +y side of sample profile
    sigma_Q_2_profile_airborne = zeros(N_Flux,1)*NaN; %total flux uncertainty associated with -y side of sample profile
    zq_profile_airborne = zeros(N_Flux,1); %e-folding height associated with profile
    zq_1_profile_airborne = zeros(N_Flux,1)*NaN; %e-folding height associated with +y side of sample profile
    zq_2_profile_airborne = zeros(N_Flux,1)*NaN; %e-folding height associated with -y side of sample profile
    tau_profile_airborne = zeros(N_Flux,1); %shear stress associated with sample profile
    ust_profile_airborne = zeros(N_Flux,1); %shear velocity associated with sample profile
    
    %go through each BSNE collection interval
    for j = 1:N_Flux
        
        %get values for profiles
        Q_profile_airborne(j) = Flux_BSNE(j).Q.Q; %total flux associated with sample profile
        sigma_Q_profile_airborne(j) = Flux_BSNE(j).Q.sigma_Q; %total flux associated with sample profile
        zq_profile_airborne(j) = Flux_BSNE(j).z.zq; %e-folding height associated with sample profile
        
        if strcmp(Sites{i},'Oceano')
            Q_1_profile_airborne(j) = Flux_BSNE(j).Q.Q_1; %total flux associated with +y side of sample profile
            sigma_Q_1_profile_airborne(j) = Flux_BSNE(j).Q.sigma_Q_1; %total flux associated with +y side of sample profile
            zq_1_profile_airborne(j) = Flux_BSNE(j).z.zq_1; %e-folding height associated with +y side of sample profile
            Q_2_profile_airborne(j) = Flux_BSNE(j).Q.Q_2; %total flux associated with -y side of sample profile
            sigma_Q_2_profile_airborne(j) = Flux_BSNE(j).Q.sigma_Q_2; %total flux associated with -y side of sample profile
            zq_2_profile_airborne(j) = Flux_BSNE(j).z.zq_2; %e-folding height associated with y side of sample profile    
        end
            
        %get values for BSNE
        q_BSNE = Flux_BSNE(j).qz.qz; %get partial flux for BSNE
        sigma_q_BSNE = Flux_BSNE(j).qz.sigma; %get partial flux uncertainty for BSNE
        y_BSNE = Flux_BSNE(j).y; %get y-position for BSNE
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
        y_profile_airborne{j} = zeros(N_samples_profile,1); %y-positions for grain size samples
        z_profile_airborne{j} = zeros(N_samples_profile,1); %heights for grain size samples
        sigma_z_profile_airborne{j} = zeros(N_samples_profile,1); %height uncertainties for grain size samples
        
        %get each airborne size distribution
        for k = 1:N_samples_profile
            try
                ind_k = find(strcmp(Name_BSNE,Name_airborne_samples(ind_samples_BSNE(k)))); %get index for BSNE
                q_profile_airborne{j}(k) = q_BSNE(ind_k); %get flux associated with BSNE
                sigma_q_profile_airborne{j}(k) = sigma_q_BSNE(ind_k); %get flux uncertainty associated with BSNE
                y_profile_airborne{j}(k) = y_BSNE(ind_k); %get y-position associated with BSNE
                z_profile_airborne{j}(k) = z_BSNE(ind_k); %get height associated with BSNE
                sigma_z_profile_airborne{j}(k) = sigma_z_BSNE(ind_k); %get height uncertainty associated with BSNE
            catch %if there is no airborne profile, set values as NaN
                q_profile_airborne{j}(k) = NaN;
                sigma_q_profile_airborne{j}(k) = NaN;
                y_profile_airborne{j}(k) = NaN;
                z_profile_airborne{j}(k) = NaN;
                sigma_z_profile_airborne{j}(k) = NaN;
            end           
        end
        
        %eliminate NaN or 0 samples
        ind_hasdata = find((~isnan(q_profile_airborne{j})) | (q_profile_airborne{j}==0));
        N_samples_profile = length(ind_hasdata);
        q_profile_airborne{j} = q_profile_airborne{j}(ind_hasdata);
        sigma_q_profile_airborne{j} = sigma_q_profile_airborne{j}(ind_hasdata);
        y_profile_airborne{j} = y_profile_airborne{j}(ind_hasdata);
        z_profile_airborne{j} = z_profile_airborne{j}(ind_hasdata);
        sigma_z_profile_airborne{j} = sigma_z_profile_airborne{j}(ind_hasdata);
        dV_profile_airborne{j} = dV_profile_airborne{j}(ind_hasdata,:);
        dVdlogd_profile_airborne{j} = dVdlogd_profile_airborne{j}(ind_hasdata,:);
    end
        
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
        ind_surface_Date = intersect(find(Date_surface>=Cluster_StartDate{i}(j)),find(Date_surface<=Cluster_EndDate{i}(j))); %indices of dates associated with cluster
        ind_surface_Location = []; %inialize list of indices of locations associated with cluster
        for k = 1:length(Cluster_Location{i}{j})
            ind_surface_Location = [ind_surface_Location, find(strcmp(Location_surface,Cluster_Location{i}{j}{k})==1)]; %get indices of locations associated with cluster
        end     
        ind_surface_Cluster = intersect(ind_surface_Date,ind_surface_Location); %get combined indices associated with cluster
        
        %compute mean dV and dVdlogd of surface samples for cluster        
        dV_bar_surface_Cluster{ind_Cluster} = mean(dV_surface(ind_surface_Cluster,:),1)./sum(mean(dV_surface(ind_surface_Cluster,:),1));
        dVdlogd_bar_surface_Cluster{ind_Cluster} = dV_bar_surface_Cluster{ind_Cluster}./dlogd_surface;      

        %get mean d10, d50, d90 of surface samples
        [d10, d50, d90] = ReferenceGrainSizes(dV_bar_surface_Cluster{ind_Cluster}, d_surface_lower, d_surface_upper);
        d10_bar_surface_Cluster(ind_Cluster) = d10;
        d50_bar_surface_Cluster(ind_Cluster) = d50;
        d90_bar_surface_Cluster(ind_Cluster) = d90;
        d10_sigma_surface_Cluster(ind_Cluster) = std([GrainSize_surface(ind_surface_Cluster).d_10_mm])./sqrt(length(ind_surface_Cluster));
        d50_sigma_surface_Cluster(ind_Cluster) = std([GrainSize_surface(ind_surface_Cluster).d_50_mm])./sqrt(length(ind_surface_Cluster));
        d90_sigma_surface_Cluster(ind_Cluster) = std([GrainSize_surface(ind_surface_Cluster).d_90_mm])./sqrt(length(ind_surface_Cluster));
   
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
        
        %% calculations for individual airborne samples
        
        %add airborne grain size bins into cluster
        d_airborne_mid_Cluster{ind_Cluster} = d_airborne_mid;
        d_airborne_lower_Cluster{ind_Cluster} = d_airborne_lower;
        d_airborne_upper_Cluster{ind_Cluster} = d_airborne_upper;
        dlogd_airborne_Cluster{ind_Cluster} = dlogd_airborne;
                
        %get indices of individual airborne samples associated with cluster
        ind_Date = intersect(find(Date_airborne_samples>=Cluster_StartDate{i}(j)),find(Date_airborne_samples<=Cluster_EndDate{i}(j))); %get indices of samples with Date in cluster
        ind_y = intersect(find(y_airborne_samples>=Cluster_ymin{i}(j)),find(y_airborne_samples<=Cluster_ymax{i}(j))); %get indices of samples with y-position in cluster
        ind_samples_airborne_Cluster = intersect(ind_Date,ind_y); %get indices of samples in Cluster
        N_samples_airborne_Cluster = length(ind_samples_airborne_Cluster); %number of airborne samples in cluster
        
        %group information about individual airborne samples
        dV_airborne_Cluster{ind_Cluster} = dV_airborne_samples(ind_samples_airborne_Cluster,:); %volume fraction of airborne individual sample
        dVdlogd_airborne_Cluster{ind_Cluster} = dVdlogd_airborne_samples(ind_samples_airborne_Cluster,:); %differential volume fraction of airborne individual sample
        d10_airborne_Cluster{ind_Cluster} = d10_airborne_samples(ind_samples_airborne_Cluster); %d10 of airborne sample profiles
        d50_airborne_Cluster{ind_Cluster} = d50_airborne_samples(ind_samples_airborne_Cluster); %d50 of airborne sample profiles
        d90_airborne_Cluster{ind_Cluster} = d90_airborne_samples(ind_samples_airborne_Cluster); %d90 of airborne sample profiles
        dbar_airborne_Cluster{ind_Cluster} = dbar_airborne_samples(ind_samples_airborne_Cluster); %dbar of airborne sample profiles
        z_airborne_Cluster{ind_Cluster} = z_airborne_samples(ind_samples_airborne_Cluster); %height profiles
        
        %% get airborne profiles for cluster
        
        %get indices of airborne profiles associated with cluster
        ind_profile_airborne_Cluster = intersect(find(Date_profile_airborne>=Cluster_StartDate{i}(j)),find(Date_profile_airborne<=Cluster_EndDate{i}(j)));
        N_profiles_airborne_Cluster = length(ind_profile_airborne_Cluster); %number of airborne profiles in cluster

        %get wind values for these profiles
        tau_profile_Cluster{ind_Cluster} = tau_profile_airborne(ind_profile_airborne_Cluster); %shear stress associated with sample profile
        ust_profile_Cluster{ind_Cluster} = ust_profile_airborne(ind_profile_airborne_Cluster); %shear velocity associated with sample profile
        taunorm_profile_Cluster{ind_Cluster} = tau_profile_airborne(ind_profile_airborne_Cluster)/tauth_Cluster(ind_Cluster); %normalized shear stress associated with sample profile
        ustnorm_profile_Cluster{ind_Cluster} = ust_profile_airborne(ind_profile_airborne_Cluster)/ustth_Cluster(ind_Cluster); %normalized shear velocity associated with sample profile
        Date_profile_Cluster{ind_Cluster} = Date_profile_airborne(ind_profile_airborne_Cluster); %date of cluster
        Time_profile_Cluster{ind_Cluster} = Time_profile_airborne(ind_profile_airborne_Cluster); %time of cluster
        
        %initialize profile and flux values
        dV_profile_airborne_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %volume fraction of airborne profile
        dVdlogd_profile_airborne_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %differential volume fraction of airborne profile
        d10_profile_airborne_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %d10 values in profiles
        d50_profile_airborne_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %d50 values in profiles
        d90_profile_airborne_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %d90 values in profiles
        dbar_profile_airborne_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %dbar values in profiles
        q_profile_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %partial fluxes associated with profile
        sigma_q_profile_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %partial flux uncertainties associated with profile
        z_profile_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %heights associated with profile
        sigma_z_profile_Cluster{ind_Cluster} = cell(N_profiles_airborne_Cluster,1); %height uncertainties associated with profile
        Q_profile_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1); %total flux associated with sample profile
        sigma_Q_profile_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1); %total flux uncertainty associated with sample profile
        zq_profile_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1); %e-folding height associated with sample profile
        
        %go through each profile, determine which values to put into cluster based on y
        for k = 1:N_profiles_airborne_Cluster
            ind_y = intersect(find(y_profile_airborne{ind_profile_airborne_Cluster(k)}>=Cluster_ymin{i}(j)),find(y_profile_airborne{ind_profile_airborne_Cluster(k)}<=Cluster_ymax{i}(j))); %get indices of samples in profile with y-position in cluster       
            dV_profile_airborne_Cluster{ind_Cluster}{k} = dV_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y,:); %volume fraction of airborne profile
            dVdlogd_profile_airborne_Cluster{ind_Cluster}{k} = dVdlogd_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y,:); %differential volume fraction of airborne profile
            q_profile_Cluster{ind_Cluster}{k} = q_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %partial fluxes associated with profile
            sigma_q_profile_Cluster{ind_Cluster}{k} = sigma_q_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %partial flux uncertainties associated with profile
            z_profile_Cluster{ind_Cluster}{k} = z_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %heights associated with profile
            sigma_z_profile_Cluster{ind_Cluster}{k} = sigma_z_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %height uncertainties associated with profile
            d10_profile_airborne_Cluster{ind_Cluster}{k} = d10_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %profile d10 of airborne samples
            d50_profile_airborne_Cluster{ind_Cluster}{k} = d50_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %profile d50 of airborne samples
            d90_profile_airborne_Cluster{ind_Cluster}{k} = d90_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %profile d90 of airborne samples
            dbar_profile_airborne_Cluster{ind_Cluster}{k} = dbar_profile_airborne{ind_profile_airborne_Cluster(k)}(ind_y); %profile dbar of airborne samples           
            
            if Cluster_ymin{i}(j)==-Inf && Cluster_ymax{i}(j)==Inf %full profile
                Q_profile_Cluster{ind_Cluster}(k) = Q_profile_airborne(ind_profile_airborne_Cluster(k)); %total flux associated with sample profile
                sigma_Q_profile_Cluster{ind_Cluster}(k) = sigma_Q_profile_airborne(ind_profile_airborne_Cluster(k)); %total flux uncertainty associated with sample profile
                zq_profile_Cluster{ind_Cluster}(k) = zq_profile_airborne(ind_profile_airborne_Cluster(k)); %e-folding height associated with sample profile
            elseif Cluster_ymin{i}(j)==0 %y+ profile
                Q_profile_Cluster{ind_Cluster}(k) = Q_1_profile_airborne(ind_profile_airborne_Cluster(k)); %total flux associated with sample profile
                sigma_Q_profile_Cluster{ind_Cluster}(k) = sigma_Q_1_profile_airborne(ind_profile_airborne_Cluster(k)); %total flux uncertainty associated with sample profile
                zq_profile_Cluster{ind_Cluster}(k) = zq_1_profile_airborne(ind_profile_airborne_Cluster(k)); %e-folding height associated with sample profile
            elseif Cluster_ymax{i}(j)==0 %y- profile
                Q_profile_Cluster{ind_Cluster}(k) = Q_2_profile_airborne(ind_profile_airborne_Cluster(k)); %total flux associated with sample profile
                sigma_Q_profile_Cluster{ind_Cluster}(k) = sigma_Q_2_profile_airborne(ind_profile_airborne_Cluster(k)); %total flux uncertainty associated with sample profile
                zq_profile_Cluster{ind_Cluster}(k) = zq_2_profile_airborne(ind_profile_airborne_Cluster(k)); %e-folding height associated with sample profile
            end
        end

        %eliminate profiles with no data points
        ind_hasdata = find(cellfun(@length,d10_profile_airborne_Cluster{ind_Cluster})>0); %indices of profiles with values in them
        N_profiles_airborne_Cluster = length(ind_hasdata); %get new number of profiles
        dV_profile_airborne_Cluster{ind_Cluster} = dV_profile_airborne_Cluster{ind_Cluster}(ind_hasdata,:);
        dVdlogd_profile_airborne_Cluster{ind_Cluster} = dVdlogd_profile_airborne_Cluster{ind_Cluster}(ind_hasdata,:);
        d10_profile_airborne_Cluster{ind_Cluster} = d10_profile_airborne_Cluster{ind_Cluster}(ind_hasdata);
        d50_profile_airborne_Cluster{ind_Cluster} = d50_profile_airborne_Cluster{ind_Cluster}(ind_hasdata);
        d90_profile_airborne_Cluster{ind_Cluster} = d90_profile_airborne_Cluster{ind_Cluster}(ind_hasdata);
        dbar_profile_airborne_Cluster{ind_Cluster} = dbar_profile_airborne_Cluster{ind_Cluster}(ind_hasdata);  
        q_profile_Cluster{ind_Cluster} = q_profile_Cluster{ind_Cluster}(ind_hasdata);
        sigma_q_profile_Cluster{ind_Cluster} = sigma_q_profile_Cluster{ind_Cluster}(ind_hasdata);
        z_profile_Cluster{ind_Cluster} = z_profile_Cluster{ind_Cluster}(ind_hasdata);
        sigma_z_profile_Cluster{ind_Cluster} = sigma_z_profile_Cluster{ind_Cluster}(ind_hasdata);
        Q_profile_Cluster{ind_Cluster} = Q_profile_Cluster{ind_Cluster}(ind_hasdata);
        sigma_Q_profile_Cluster{ind_Cluster} = sigma_Q_profile_Cluster{ind_Cluster}(ind_hasdata);
        zq_profile_Cluster{ind_Cluster} = zq_profile_Cluster{ind_Cluster}(ind_hasdata);
        tau_profile_Cluster{ind_Cluster} = tau_profile_Cluster{ind_Cluster}(ind_hasdata);
        ust_profile_Cluster{ind_Cluster} = ust_profile_Cluster{ind_Cluster}(ind_hasdata);
        taunorm_profile_Cluster{ind_Cluster} = taunorm_profile_Cluster{ind_Cluster}(ind_hasdata);
        ustnorm_profile_Cluster{ind_Cluster} = ustnorm_profile_Cluster{ind_Cluster}(ind_hasdata);  
        Date_profile_Cluster{ind_Cluster} = Date_profile_Cluster{ind_Cluster}(ind_hasdata);
        Time_profile_Cluster{ind_Cluster} = Time_profile_Cluster{ind_Cluster}(ind_hasdata);
        
        %calculate z/zq values
        zq_bar_Cluster(ind_Cluster) = mean(zq_profile_Cluster{ind_Cluster}(cellfun(@length,z_profile_Cluster{ind_Cluster})>N_min_profile)); %compute mean zq for Cluster - needed to compute znorm         
        for k = 1:N_profiles_airborne_Cluster %go through profiles
            znorm_profile_Cluster{ind_Cluster}{k} = z_profile_Cluster{ind_Cluster}{k}./zq_bar_Cluster(ind_Cluster); %get z/zq for profile
        end
                        
        %determine which profiles are usable for profile analyses
        ind_sufficient_data = find(cellfun(@length,z_profile_Cluster{ind_Cluster})>=N_min_profile); %indices of profiles with sufficient data
        ind_sufficient_range = find(cellfun(@min,znorm_profile_Cluster{ind_Cluster})<=znorm_maxlower & cellfun(@max,znorm_profile_Cluster{ind_Cluster})>=znorm_minupper); %indices of profiles covering sufficient range
        ind_usable_profile_Cluster{ind_Cluster} = intersect(ind_sufficient_data, ind_sufficient_range); %set to 1 if profile is usable for analysis
             
        %% get normalized height profiles for all traps
        znorm_airborne_Cluster{ind_Cluster} = z_airborne_Cluster{ind_Cluster}./zq_bar_Cluster(ind_Cluster);
        
        %% compute profile mean values
                
        %initialize profile mean values
        dV_profilebar_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,N_d_airborne)*NaN; %profile mean volume fraction of airborne samples
        dVdlogd_profilebar_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,N_d_airborne)*NaN; %profile mean differential volume fraction of airborne samples
        d10_profilebar_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile mean d10 of airborne samples
        d50_profilebar_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile mean d50 of airborne samples
        d90_profilebar_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile mean d90 of airborne samples
        dbar_profilebar_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile mean dbar of airborne samples
        d10_profilesigma_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile std d10 of airborne samples
        d50_profilesigma_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile std d50 of airborne samples
        d90_profilesigma_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile std d90 of airborne samples
        dbar_profilesigma_airborne_Cluster{ind_Cluster} = zeros(N_profiles_airborne_Cluster,1)*NaN; %profile std dbar of airborne samples        
           
        %calculate profile mean values 
        for k = 1:N_profiles_airborne_Cluster
            
            %get profiles
            z_profile = z_profile_Cluster{ind_Cluster}{k}; %get z profile
            zq_profile = zq_profile_Cluster{ind_Cluster}(k); %get zq for profile
            znorm_profile = znorm_profile_Cluster{ind_Cluster}{k}; %get z/zq for profile
            q_profile = q_profile_Cluster{ind_Cluster}{k}; %get q profile
            dV_profile = dV_profile_airborne_Cluster{ind_Cluster}{k}; %get dV profile
            d10_profile = d10_profile_airborne_Cluster{ind_Cluster}{k}; %get d10 profile
            d50_profile = d50_profile_airborne_Cluster{ind_Cluster}{k}; %get d50 profile
            d90_profile = d90_profile_airborne_Cluster{ind_Cluster}{k}; %get d90 profile
            dbar_profile = dbar_profile_airborne_Cluster{ind_Cluster}{k}; %get dbar profile
            N_samples_profile = length(z_profile); %get number of samples in profile           
           
            %combine airborne samples
            dV_weighted_matrix = zeros(N_samples_profile,N_d_airborne); %initialize matrix of dV's multiplied by weights            
            for m = 1:N_samples_profile %go through each sample to compute flux-weighted profile mean grain-size distribution              
                %get height below for weighting
                z_below = max(z_profile(z_profile<z_profile(m)));
                if isempty(z_below)
%                     z_below = 0;
                    z_below = zq_profile; %set lower limit for weighting zq for profile
                end
                %get height above for weighting
                z_above = min(z_profile(z_profile>z_profile(m)));
                if isempty(z_above)
                    z_above = inf;
                end
                %compute weight for sample - if multiple samples at specific height, distribute weight evenly
                weight_m = (exp(-z_below/zq_profile)-exp(-z_above/zq_profile))/length(find(z_profile==z_profile(m)));

                %compute weighted values for row
                dV_weighted_matrix(m,:) = dV_profile(m,:).*weight_m;
            end

            %compute mean dV profile
            dV_profilebar = sum(dV_weighted_matrix)./sum(sum(dV_weighted_matrix)); %add weighted rows to get flux-weighted grain-size distribution, normalize so that it sums to 1
            dVdlogd_profilebar = dV_profilebar_airborne_Cluster{ind_Cluster}(k,:)./dlogd_airborne; %convert to dVdlogd                

            %add these to cell arrays
            dV_profilebar_airborne_Cluster{ind_Cluster}(k,:) = dV_profilebar;
            dVdlogd_profilebar_airborne_Cluster{ind_Cluster}(k,:) = dVdlogd_profilebar; 

            %get profile mean airborne size distribution reference sizes
            [d10, d50, d90] = ReferenceGrainSizes(dV_profilebar, d_airborne_lower, d_airborne_upper);
            d10_profilebar_airborne_Cluster{ind_Cluster}(k) = d10; %d10 of airborne profile
            d50_profilebar_airborne_Cluster{ind_Cluster}(k) = d50; %d50 of airborne profile
            d90_profilebar_airborne_Cluster{ind_Cluster}(k) = d90; %d90 of airborne profile
            dbar_profilebar_airborne_Cluster{ind_Cluster}(k) = exp(sum(dV_profilebar.*log(d_airborne_mid))); %compute weighted geometric mean of particle size bins
            d10_profilesigma_airborne_Cluster{ind_Cluster}(k) = std(d10_profile); %std of d10 in profile
            d50_profilesigma_airborne_Cluster{ind_Cluster}(k) = std(d50_profile); %std of d50 in profile
            d90_profilesigma_airborne_Cluster{ind_Cluster}(k) = std(d90_profile); %std of d90 in profile
            dbar_profilesigma_airborne_Cluster{ind_Cluster}(k) = std(dbar_profile); %std of dbar in profile
        end

        %compute mean dV and dVdlogd of airborne samples for cluster        
        dV_bar_airborne_Cluster{ind_Cluster} = mean(dV_profilebar_airborne_Cluster{ind_Cluster}(ind_usable_profile_Cluster{ind_Cluster},:),1);
        dVdlogd_bar_airborne_Cluster{ind_Cluster} = dV_bar_airborne_Cluster{ind_Cluster}./dlogd_airborne; 
                
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE UNBINNED SIZE CONDITIONED FLUXS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_surface_Cluster = cell(N_Cluster,1); %fraction of surface sample in bin
f_airborne_profile_Cluster = cell(N_Cluster,1); %fraction of airborne sample in bin
f_airborne_profilebar_Cluster = cell(N_Cluster,1); %fraction of airborne sample in bin
f_airborne_bar_Cluster = cell(N_Cluster,1); %mean fraction of airborne sample in bin
f_ratio_Cluster = cell(N_Cluster,1); %ratio of airborne and surface
f_ratio_bar_Cluster = cell(N_Cluster,1); %mean ratio of airborne and surface
f_ratio_sigma_Cluster = cell(N_Cluster,1); %uncertainty on ratio of airborne and surface
f_ratio_sigmalog_Cluster = cell(N_Cluster,1); %uncertainty on log of ratio of airborne and surface
d_f_lower_Cluster = cell(N_Cluster,1); %lower limit d for cluster
d_f_upper_Cluster = cell(N_Cluster,1); %upper limit d for cluster
d_f_mid_Cluster = cell(N_Cluster,1); %mid d for cluster
qi_Cluster = cell(N_Cluster,1); %size-selective partial flux in cluster
sigma_qi_Cluster = cell(N_Cluster,1); %size-selective partial flux uncertainty in cluster
Qi_Cluster = cell(N_Cluster,1); %size-selective flux in cluster
sigma_Qi_Cluster = cell(N_Cluster,1); %size-selective flux uncertainty in cluster
zqi_Cluster = cell(N_Cluster,1); %size-selective zq in cluster
sigma_zqi_Cluster = cell(N_Cluster,1); %size-selective zq uncertainty in cluster
zqinorm_Cluster = cell(N_Cluster,1); %size-selective zqinorm in cluster
sigma_zqinorm_Cluster = cell(N_Cluster,1); %size-selective zqnorm uncertainty in cluster
Qhat_Cluster = cell(N_Cluster,1); %normalized size-selective flux in cluster
sigma_Qhat_Cluster = cell(N_Cluster,1); %normalized size-selective flux uncertainty in cluster
zqi_bar_Cluster = cell(N_Cluster,1); %mean zqi in Cluster
zqi_sigma_Cluster = cell(N_Cluster,1); %uncertainty on zqi in Cluster
zqinorm_bar_Cluster = cell(N_Cluster,1); %mean zqi/d50 in Cluster
zqinorm_sigma_Cluster = cell(N_Cluster,1); %sigma zqi/d50 in Cluster

for i = 1:N_Cluster
    
    %get grain sizes
    ind_dmin = find(d_surface_mid_Cluster{i} >= d_min); %get indices of values for which d>dmin
    ind_drange = find(d_surface_mid_Cluster{i} >= d_min & d_surface_mid_Cluster{i} <= d_max); %get only values in dmin<d<dmax
    N_drange = length(ind_drange); %get number of d's in this range
    d_f_lower_Cluster{i} = d_surface_lower_Cluster{i}(ind_drange);
    d_f_upper_Cluster{i} = d_surface_upper_Cluster{i}(ind_drange);
    d_f_mid_Cluster{i} = d_surface_mid_Cluster{i}(ind_drange);

    %get fraction of grains in each surface bin
    f_surface_total = sum(dV_bar_surface_Cluster{i}(ind_dmin)); %get total f_surface as fraction among minimum grain size
    f_surface_Cluster{i} = dV_bar_surface_Cluster{i}(ind_drange)./f_surface_total;
    
    %get fraction of grains in each airborne bin
    N_profile_Cluster = length(Q_profile_Cluster{i}); %get number of profiles
    f_airborne_profile_Cluster{i} = cell(N_profile_Cluster,1); %initialize matrix of f_airborne
    f_airborne_profilebar_Cluster{i} = zeros(N_profile_Cluster,N_drange); %initialize matrix of f_airborne_profilebar
    f_ratio_Cluster{i} = zeros(N_profile_Cluster,N_drange); %initialize matrix of f_ratio
    for k = 1:N_profile_Cluster
        N_z = length(z_profile_Cluster{i}{k}); %number of airborne z
        f_airborne_profile_Cluster{i}{k} = zeros(N_z,N_drange);
        for m = 1:N_z
            f_airborne_total = sum(dV_profile_airborne_Cluster{i}{k}(m,ind_dmin)); %get total f_airborne as fraction among minimum grain size
            f_airborne_profile_Cluster{i}{k}(m,:) = dV_profile_airborne_Cluster{i}{k}(m,ind_drange)./f_airborne_total;
        end
        f_airborne_profilebar_total = sum(dV_profilebar_airborne_Cluster{i}(k,ind_dmin)); %get total f_airborne as fraction among minimum grain size
        f_airborne_profilebar_Cluster{i}(k,:) = dV_profilebar_airborne_Cluster{i}(k,ind_drange)./f_airborne_profilebar_total;
        f_ratio_Cluster{i}(k,:) = f_airborne_profilebar_Cluster{i}(k,:)./f_surface_Cluster{i}; %get f_ratio
    end
    
    %get mean fairborne and fratio
    ind_usable = ind_usable_profile_Cluster{i}; %get indices of usable data
    f_airborne_bar_Cluster{i} = mean(f_airborne_profilebar_Cluster{i}(ind_usable,:));
    f_ratio_bar_Cluster{i} = mean(f_ratio_Cluster{i}(ind_usable,:));
    f_ratio_sigma_Cluster{i} = std(f_ratio_Cluster{i}(ind_usable,:))./sqrt(length(ind_usable)); %uncertainty on ratio of airborne and surface
    f_ratio_sigmalog_Cluster{i} = zeros(1,N_drange); %initialize list of log uncertainty
    for k = 1:N_drange
        f_ratio_k = f_ratio_Cluster{i}(ind_usable,k);
        f_ratio_sigmalog_Cluster{i}(k) = std(log(f_ratio_k(f_ratio_k>0)))./sqrt(sum(f_ratio_k>0)); %uncertainty on log of ratio of airborne and surface
    end
        
    %go through each sample in each profile to get size-selective partial fluxes
    qi_Cluster{i} = cell(N_profile_Cluster,1); %size-selective partial flux in cluster
    sigma_qi_Cluster{i} = cell(N_profile_Cluster,1); %size-selective partial flux uncertainty in cluster
    for k = 1:N_profile_Cluster
    	N_z = length(z_profile_Cluster{i}{k});
        qi_Cluster{i}{k} = zeros(N_z,N_drange);
        sigma_qi_Cluster{i}{k} = zeros(N_z,N_drange);
        for m = 1:N_z
            qi_Cluster{i}{k}(m,:) = q_profile_Cluster{i}{k}(m)*f_airborne_profile_Cluster{i}{k}(m,:); %size-selective partial flux in each bin
            sigma_qi_Cluster{i}{k}(m,:) = sigma_q_profile_Cluster{i}{k}(m)*f_airborne_profile_Cluster{i}{k}(m,:); %size-selective partial flux uncertainty in each bin
        end
    end
    
    %compute size-selective total fluxes and zq
    Qi_Cluster{i} = zeros(N_profile_Cluster,N_drange); %size-selective flux in cluster
    sigma_Qi_Cluster{i} = zeros(N_profile_Cluster,N_drange); %size-selective flux uncertainty in cluster
    zqi_Cluster{i} = zeros(N_profile_Cluster,N_drange); %size-selective zq in cluster
    sigma_zqi_Cluster{i} = zeros(N_profile_Cluster,N_drange); %size-selective zq uncertainty in cluster
    zqinorm_Cluster{i} = zeros(N_profile_Cluster,N_drange); %size-selective zq/d in cluster
    sigma_zqinorm_Cluster{i} = zeros(N_profile_Cluster,N_drange); %size-selective zq/d uncertainty in cluster
    Qhat_Cluster{i} = zeros(N_profile_Cluster,N_drange); %normalized size-selective flux in cluster
    sigma_Qhat_Cluster{i} = zeros(N_profile_Cluster,N_drange); %normalized size-selective flux uncertainty in cluster

    %fit size-selective profiles
    for j = 1:N_profile_Cluster
        for k = 1:N_drange
            [~,zqi,Qi,~,sigma_zqi,sigma_Qi] = qz_profilefit_exponential(qi_Cluster{i}{j}(:,k), z_profile_Cluster{i}{j}, sigma_qi_Cluster{i}{j}(:,k), sigma_z_profile_Cluster{i}{j}, zq_profile_Cluster{i}(j)); %perform profile fitting
            Qi_Cluster{i}(j,k) = Qi;
            sigma_Qi_Cluster{i}(j,k) = sigma_Qi;
            zqi_Cluster{i}(j,k) = zqi;
            sigma_zqi_Cluster{i}(j,k) = sigma_zqi;
            zqinorm_Cluster{i}(j,k) = 1000*zqi/d_f_mid_Cluster{i}(k);
            sigma_zqinorm_Cluster{i}(j,k) = 1000*sigma_zqi/d_f_mid_Cluster{i}(k);
        end
        Qhat_Cluster{i}(j,:) = Qi_Cluster{i}(j,:)./f_surface_Cluster{i}; %compute Qhat
        sigma_Qhat_Cluster{i}(j,:) = sigma_Qi_Cluster{i}(j,:)./f_surface_Cluster{i}; %compute sigma_Qhat
    end
            
    %compute mean zq and zqnorm
    zqi_bar_Cluster{i} = zeros(1,N_drange)*NaN; %mean zqi in Cluster
    zqi_sigma_Cluster{i} = zeros(1,N_drange)*NaN; %uncertainty on zqi in Cluster
    zqinorm_bar_Cluster{i} = zeros(1,N_drange)*NaN; %mean zqi/d50 in Cluster
    zqinorm_sigma_Cluster{i} = zeros(1,N_drange)*NaN; %sigma zqi/d50 in Cluster
    ind_usable = ind_usable_profile_Cluster{i}; %get indices of usable data
    for j = 1:N_drange
        ind_mean = intersect(ind_usable,find(~isnan(zqi_Cluster{i}(:,j))));
        if length(ind_mean)>=3 %only compute if there are enough values
            zqi_bar_Cluster{i}(j) = mean(zqi_Cluster{i}(ind_mean,j)); %mean zqi in Cluster
            zqi_sigma_Cluster{i}(j) = std(zqi_Cluster{i}(ind_mean,j))/sqrt(length(ind_mean)); %uncertainty on zqi in Cluster
            zqinorm_bar_Cluster{i}(j) = mean(zqinorm_Cluster{i}(ind_mean,j)); %mean zqi/d50 in Cluster
            zqinorm_sigma_Cluster{i}(j) = std(zqinorm_Cluster{i}(ind_mean,j))/sqrt(length(ind_mean)); %sigma zqi/d50 in Cluster
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP SIZE DISTRIBUTIONS BY TAU/TAUTH BINS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_taunormbar_airborne_mid_Cluster = cell(N_Cluster,1); %grain sizes for this analysis, include only grain size bins with d>d_min and d<2
dV_taunormbar_airborne_Cluster = cell(N_Cluster,1); %all values
dVdlogd_taunormbar_airborne_Cluster = cell(N_Cluster,1); %all values
d_taunormbar_airborne_mid_drange_Cluster{i} = cell(N_Cluster,1); %values restricted to drange
dVdlogd_taunormbar_airborne_drange_Cluster{i} = cell(N_Cluster,1); %values restricted to drange
dbar_taunormbar_airborne_Cluster = cell(N_Cluster,1); %d50 for each taunormbar

for i = 1:N_Cluster    
    
    %get all values
    N_d = length(d_airborne_mid_Cluster{i}); %number of grain size bins
    d_taunormbar_airborne_mid_Cluster{i} = d_airborne_mid_Cluster{i}; %grain sizes for this analysis, include only grain size bins with d>d_min
    dV_taunormbar_airborne_Cluster{i} = zeros(N_taunorm_bins,N_d);
    dVdlogd_taunormbar_airborne_Cluster{i} = zeros(N_taunorm_bins,N_d);
    dbar_taunormbar_airborne_Cluster{i} = zeros(N_taunorm_bins,1);
    
    %restrict to "usable" values
    ind_usable = ind_usable_profile_Cluster{i};
    for j = 1:N_taunorm_bins
        ind_taunorm = find(taunorm_profile_Cluster{i}>=taunorm_min_bins(j) & taunorm_profile_Cluster{i}<=taunorm_max_bins(j));
        ind_average = intersect(ind_taunorm,ind_usable);
        dV_taunormbar_airborne_Cluster{i}(j,:) = mean(dV_profilebar_airborne_Cluster{i}(ind_average,:),1); %normalize by volume fraction > d_min        
        dVdlogd_taunormbar_airborne_Cluster{i}(j,:) = dV_taunormbar_airborne_Cluster{i}(j,:)./dlogd_airborne_Cluster{i}; %normalize by volume fraction > d_min
        dbar_taunormbar_airborne_Cluster{i}(j) = sum(dV_taunormbar_airborne_Cluster{i}(j,:).*d_taunormbar_airborne_mid_Cluster{i}); %mean grain size for each cluster
    end
    
    %get values restricted to drange
    ind_drange = find(d_airborne_lower_Cluster{i} >= d_min & d_airborne_upper_Cluster{i}<=2); %include only grain size bins with d>d_min and d<2
    d_taunormbar_airborne_mid_drange_Cluster{i} = d_taunormbar_airborne_mid_Cluster{i}(ind_drange);
    dVdlogd_taunormbar_airborne_drange_Cluster{i} = dVdlogd_taunormbar_airborne_Cluster{i}(:,ind_drange);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET ZQNORM ASSOCIATED WITH DHAT = FINE, MEDIUM, AND COARSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zqinorm_dhat_fine_Cluster = cell(N_Cluster,1); %zqi/d associated with dhat_fine
sigma_zqinorm_dhat_fine_Cluster = cell(N_Cluster,1); %uncertainty in zqi/d associated with dhat_fine
zqinorm_dhat_medium_Cluster = cell(N_Cluster,1); %zqi/d associated with dhat_medium
sigma_zqinorm_dhat_medium_Cluster = cell(N_Cluster,1); %uncertainty in zqi/d associated with dhat_medium
zqinorm_dhat_coarse_Cluster = cell(N_Cluster,1); %zqi/d associated with dhat_coarse
sigma_zqinorm_dhat_coarse_Cluster = cell(N_Cluster,1); %uncertainty in zqi/d associated with dhat_coarse

dhat_fine = 0.6;
dhat_medium = 0.9;
dhat_coarse = 1.2;

for i = 1:N_Cluster
    %get zqinorm for dhat_fine
    dhat_fine_diff = abs(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i)-dhat_fine);
    ind_dhat_fine = find(dhat_fine_diff==min(dhat_fine_diff));
    zqinorm_dhat_fine_Cluster{i} = zqinorm_Cluster{i}(:,ind_dhat_fine);
    sigma_zqinorm_dhat_fine_Cluster{i} = sigma_zqinorm_Cluster{i}(:,ind_dhat_fine);
    
    %get zqinorm for dhat_medium
    dhat_medium_diff = abs(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i)-dhat_medium);
    ind_dhat_medium = find(dhat_medium_diff==min(dhat_medium_diff));
    zqinorm_dhat_medium_Cluster{i} = zqinorm_Cluster{i}(:,ind_dhat_medium);
    sigma_zqinorm_dhat_medium_Cluster{i} = sigma_zqinorm_Cluster{i}(:,ind_dhat_medium);
    
    %get zqinorm for dhat_coarse
    dhat_coarse_diff = abs(d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i)-dhat_coarse);
    ind_dhat_coarse = find(dhat_coarse_diff==min(dhat_coarse_diff));
    zqinorm_dhat_coarse_Cluster{i} = zqinorm_Cluster{i}(:,ind_dhat_coarse);
    sigma_zqinorm_dhat_coarse_Cluster{i} = sigma_zqinorm_Cluster{i}(:,ind_dhat_coarse);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT TO ZQINORM VS TAUNORM %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize values - fine
a_zqinorm_taunorm_fine = zeros(N_Cluster,1); %fitting intercept
b_zqinorm_taunorm_fine = zeros(N_Cluster,1); %fitting slope
sigma_a_zqinorm_taunorm_fine = zeros(N_Cluster,1); %fitting intercept uncertainty
sigma_b_zqinorm_taunorm_fine = zeros(N_Cluster,1); %fitting slope uncertainty
taunorm_fit_zqinorm_fine = cell(N_Cluster,1); %taunorm values for fitting
zqinorm_pred_fine = cell(N_Cluster,1); %predicted zqinorm values for fitting
sigma_zqinorm_pred_fine = cell(N_Cluster,1); %uncertainty on predicted zqinorm values for fitting

% initialize values - medium
a_zqinorm_taunorm_medium = zeros(N_Cluster,1); %fitting intercept
b_zqinorm_taunorm_medium = zeros(N_Cluster,1); %fitting slope
sigma_a_zqinorm_taunorm_medium = zeros(N_Cluster,1); %fitting intercept uncertainty
sigma_b_zqinorm_taunorm_medium = zeros(N_Cluster,1); %fitting slope uncertainty
taunorm_fit_zqinorm_medium = cell(N_Cluster,1); %taunorm values for fitting
zqinorm_pred_medium = cell(N_Cluster,1); %predicted zqinorm values for fitting
sigma_zqinorm_pred_medium = cell(N_Cluster,1); %uncertainty on predicted zqinorm values for fitting

% initialize values - coarse
a_zqinorm_taunorm_coarse = zeros(N_Cluster,1); %fitting intercept
b_zqinorm_taunorm_coarse = zeros(N_Cluster,1); %fitting slope
sigma_a_zqinorm_taunorm_coarse = zeros(N_Cluster,1); %fitting intercept uncertainty
sigma_b_zqinorm_taunorm_coarse = zeros(N_Cluster,1); %fitting slope uncertainty
taunorm_fit_zqinorm_coarse = cell(N_Cluster,1); %taunorm values for fitting
zqinorm_pred_coarse = cell(N_Cluster,1); %predicted zqinorm values for fitting
sigma_zqinorm_pred_coarse = cell(N_Cluster,1); %uncertainty on predicted zqinorm values for fitting

for i = 1:N_Cluster
    
    %get all values for fitting
    ind_usable = ind_usable_profile_Cluster{i}; %get indices of all usable taunorm
    taunorm_usable = taunorm_profile_Cluster{i}(ind_usable); %get these taunorm
    ind_abovethreshold = find(taunorm_usable>=1); %determine which ones are above threshold
    ind_usable = ind_usable(ind_abovethreshold); %replace ind_usable list with only those above threshold
    taunorm_usable = taunorm_usable(ind_abovethreshold); %replace taunorm_usable list with only those above threshold
    [taunorm_usable,ind_sort] = sort(taunorm_usable); %sort taunorm_usable list
    ind_usable = ind_usable(ind_sort); %replace ind_usable list with list that is sorted by taunorm  
    zqinorm_usable_fine = zqinorm_dhat_fine_Cluster{i}(ind_usable); %get zqinorm_usable_fine values associated with taunorm_usable 
    sigma_zqinorm_usable_fine = sigma_zqinorm_dhat_fine_Cluster{i}(ind_usable); %get sigma_zqinorm_usable_fine values associated with taunorm_usable 
    zqinorm_usable_medium = zqinorm_dhat_medium_Cluster{i}(ind_usable); %get zqinorm_usable_medium values associated with taunorm_usable 
    sigma_zqinorm_usable_medium = sigma_zqinorm_dhat_medium_Cluster{i}(ind_usable); %get sigma_zqinorm_usable_medium values associated with taunorm_usable 
    zqinorm_usable_coarse = zqinorm_dhat_coarse_Cluster{i}(ind_usable); %get zqinorm_usable_coarse values associated with taunorm_usable 
    sigma_zqinorm_usable_coarse = sigma_zqinorm_dhat_coarse_Cluster{i}(ind_usable); %get sigma_zqinorm_usable_coarse values associated with taunorm_usable 
    
    %get indices of non-NaN values for fitting - fine
    ind_fit = find(~isnan(zqinorm_usable_fine)); %get indices of only real zqinorm values for fitting
    taunorm_fit_zqinorm_fine{i} = taunorm_usable(ind_fit);
    zqinorm_fit_fine = zqinorm_usable_fine(ind_fit);
    sigma_zqinorm_fit_fine = sigma_zqinorm_usable_fine(ind_fit);
    
    %get indices of non-NaN values for fitting - medium
    ind_fit = find(~isnan(zqinorm_usable_medium)); %get indices of only real zqinorm values for fitting
    taunorm_fit_zqinorm_medium{i} = taunorm_usable(ind_fit);
    zqinorm_fit_medium = zqinorm_usable_medium(ind_fit);
    sigma_zqinorm_fit_medium = sigma_zqinorm_usable_medium(ind_fit);

    %get indices of non-NaN values for fitting - coarse
    ind_fit = find(~isnan(zqinorm_usable_coarse)); %get indices of only real zqinorm values for fitting
    taunorm_fit_zqinorm_coarse{i} = taunorm_usable(ind_fit);
    zqinorm_fit_coarse = zqinorm_usable_coarse(ind_fit);
    sigma_zqinorm_fit_coarse = sigma_zqinorm_usable_coarse(ind_fit);
    
    %perform fitting - fine
    [a_fine, b_fine, sigma_a_fine, sigma_b_fine, zqinorm_pred_values_fine, sigma_zqinorm_pred_values_fine] = ...
        linearfit(taunorm_fit_zqinorm_fine{i},zqinorm_fit_fine,sigma_zqinorm_fit_fine);

    %perform fitting - medium
    [a_medium, b_medium, sigma_a_medium, sigma_b_medium, zqinorm_pred_values_medium, sigma_zqinorm_pred_values_medium] = ...
        linearfit(taunorm_fit_zqinorm_medium{i},zqinorm_fit_medium,sigma_zqinorm_fit_medium);
    
    %perform fitting - coarse
    [a_coarse, b_coarse, sigma_a_coarse, sigma_b_coarse, zqinorm_pred_values_coarse, sigma_zqinorm_pred_values_coarse] = ...
        linearfit(taunorm_fit_zqinorm_coarse{i},zqinorm_fit_coarse,sigma_zqinorm_fit_coarse);
    
    %assign values to lists - fine
    a_zqinorm_taunorm_fine(i) = a_fine;
    b_zqinorm_taunorm_fine(i) = b_fine;
    sigma_a_zqinorm_taunorm_fine(i) = sigma_a_fine;
    sigma_b_zqinorm_taunorm_fine(i) = sigma_b_fine;
    zqinorm_pred_fine{i} = zqinorm_pred_values_fine;
    sigma_zqinorm_pred_fine{i} = sigma_zqinorm_pred_values_fine;
    
    %assign values to lists - medium
    a_zqinorm_taunorm_medium(i) = a_medium;
    b_zqinorm_taunorm_medium(i) = b_medium;
    sigma_a_zqinorm_taunorm_medium(i) = sigma_a_medium;
    sigma_b_zqinorm_taunorm_medium(i) = sigma_b_medium;
    zqinorm_pred_medium{i} = zqinorm_pred_values_medium;
    sigma_zqinorm_pred_medium{i} = sigma_zqinorm_pred_values_medium;
    
    %assign values to lists - coarse
    a_zqinorm_taunorm_coarse(i) = a_coarse;
    b_zqinorm_taunorm_coarse(i) = b_coarse;
    sigma_a_zqinorm_taunorm_coarse(i) = sigma_a_coarse;
    sigma_b_zqinorm_taunorm_coarse(i) = sigma_b_coarse;
    zqinorm_pred_coarse{i} = zqinorm_pred_values_coarse;
    sigma_zqinorm_pred_coarse{i} = sigma_zqinorm_pred_values_coarse;
end

%%
%%%%%%%%%%%%%%%%
% NEW ANALYSES %
%%%%%%%%%%%%%%%%

znorm_lower_taunorm_bin = [0.75; 1; 1.35; 1.85; 2.45; 3.30; 4.45];
znorm_upper_taunorm_bin = [1; 1.35; 1.85; 2.45; 3.30; 4.45; 6.00];
znorm_mid_taunorm_bin = geomean([znorm_lower_taunorm_bin'; znorm_upper_taunorm_bin'])';
N_znorm_bins = length(znorm_mid_taunorm_bin);

%% get taunorm-conditioned mean grain size profiles at each site
dbar_profile_airborne_taunorm_Cluster = cell(N_Cluster,1); %initialize binned grain size profiles for all clusters
znorm_profile_airborne_taunorm_Cluster = cell(N_Cluster,1); %initialize heights for binned grain size profiles for all clusters

for i = 1:N_Cluster
    
    dbar_profile_airborne_taunorm_Cluster{i} = zeros(N_taunorm_bins,N_znorm_bins)*NaN; %initialize binned grain size profiles for this cluster
    znorm_profile_airborne_taunorm_Cluster{i} = zeros(N_taunorm_bins,N_znorm_bins)*NaN; %initialize heights for binned grain size profiles for all clusters

    %go through each taunorm bin
    ind_usable = ind_usable_profile_Cluster{i}; %get list of "usable" values
    for j = 1:N_taunorm_bins
        ind_taunorm = find(taunorm_profile_Cluster{i}>=taunorm_min_bins(j) & taunorm_profile_Cluster{i}<=taunorm_max_bins(j));
        ind_taunormbin = intersect(ind_taunorm,ind_usable); %get indices of profiles for bin
                
        %go through each height for this bin
        for k = 1:N_znorm_bins
            znorm_taunormbin = []; %intialize list of zs for this bin
            dbar_taunormbin = []; %initialize list of dbars for this bin

            %go through each profile for this height
            for l = 1:length(ind_taunormbin)
                znorm_profile = znorm_profile_Cluster{i}{ind_taunormbin(l)}; %get znorm profiles for u* bin
                ind_znorm = intersect(find(znorm_profile>znorm_lower_taunorm_bin(k)),...
                    find(znorm_profile<znorm_upper_taunorm_bin(k))); %get ind of znorm for height
                znorm_taunormbin = [znorm_taunormbin; znorm_profile(ind_znorm)]; %add znorm to list (if it exists)
                dbar_taunormbin = [dbar_taunormbin; dbar_profile_airborne_Cluster{i}{ind_taunormbin(l)}(ind_znorm)]; %add dbar to list (if it exists)
            end
            if ~isempty(znorm_taunormbin)
                znorm_profile_airborne_taunorm_Cluster{i}(j,k) = mean(znorm_taunormbin(~isnan(dbar_taunormbin))); %get mean znorm for taunorm bin
                dbar_profile_airborne_taunorm_Cluster{i}(j,k) = mean(dbar_taunormbin(~isnan(dbar_taunormbin))); %get mean dbar for taunorm bin
            end
        end
    end
end

%%
%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%

%% save data
AnalysisData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis');
save(AnalysisData_Path,'*Cluster','Cluster*','*bins','binning_type',...
    'dref_type','*min','*max',...
    '*fine','*medium','*coarse','*fit');