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
PrimaryData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis');
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(GrainSizeData_Path); %grain size data
load(SaltationFluxData_Path); %saltation data
load(PrimaryData_Path); %primary size-selective data
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

%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS FOR FITTING %
%%%%%%%%%%%%%%%%%%%%%%%%%%
taunorm_max_Qhat_tau_fit = 2; %maximum normalized shear stress for fitting Qhat versus tau


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP ZQNORM AND FRATIO BY DHAT BINS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zqinorm_bin_Cluster = cell(N_Cluster,1); %separate zq/d into bins
taunorm_zqinorm_bin_Cluster = cell(N_Cluster,1); %get associated taunorm
for i = 1:N_Cluster
    ind_usable = ind_usable_profile_Cluster{i};
    N_profiles = length(ind_usable);
    zqinorm_bin_Cluster{i} = zeros(N_profiles,N_bins)*NaN;
    taunorm_zqinorm_bin_Cluster{i} = taunorm_profile_Cluster{i}(ind_usable);
    
    %get normalized d
    if strcmp(dref_type,'d50')==1
        dhat = d_f_mid_Cluster{i}./d50_bar_surface_Cluster(i);
    elseif strcmp(dref_type,'dmodal')==1
        dhat = d_f_mid_Cluster{i}./dmodal_bar_surface_Cluster(i);
    end

    %compute binned zqinorm
    for j = 1:length(ind_usable)
        for k = 1:N_bins
            ind_drange = find(dhat >= dhat_bin_lower_Cluster(i,k) & dhat<=dhat_bin_upper_Cluster(i,k)); %include only grain size bins with d>d_min and d<2
            zqinorm_bin = zqinorm_Cluster{i}(ind_usable(j),ind_drange); %get zqinorm values in bin
            ind_mean = find(~isnan(zqinorm_bin)); %find the values that are defined
            if length(ind_mean)>=3
                zqinorm_bin_Cluster{i}(j,k) = mean(zqinorm_bin(ind_mean));
            end
        end        
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
    
    %compute size-selective total fluxes
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
C_Qhat_tau_fit = cell(N_Cluster,1);
tauth_Qhat_tau_fit = cell(N_Cluster,1);
sigma_C_Qhat_tau_fit = cell(N_Cluster,1);
sigma_tauth_Qhat_tau_fit = cell(N_Cluster,1);
tau_Qhat_fit = cell(N_Cluster,1);
Qhat_fit = cell(N_Cluster,1);
sigma_Qhat_fit = cell(N_Cluster,1);

for i = 1:N_Cluster
    
    %inialize fit values
    C_Qhat_tau_fit{i} = zeros(N_bins,1);
    tauth_Qhat_tau_fit{i} = zeros(N_bins,1);
    sigma_C_Qhat_tau_fit{i} = zeros(N_bins,1);
    sigma_tauth_Qhat_tau_fit{i} = zeros(N_bins,1);
    tau_Qhat_fit{i} = cell(N_bins,1);
    Qhat_fit{i} = cell(N_bins,1);  
    sigma_Qhat_fit{i} = cell(N_bins,1);  
    
    %go through bins
    for k = 1:N_bins
        ind_fit = intersect(find(~isnan(Qhat_Cluster{i}(:,k))),find(taunorm_profile_Cluster{i}<=taunorm_max_Qhat_tau_fit)); %get indices of non-NaN values and taunorm < taunorm_max for fitting
        taunorm_max_Qhat_tau_fit = 2; %maximum normalized shear stress for fitting Qhat versus tau

        tau_Qhat_fit{i}{k} = tau_profile_Cluster{i}(ind_fit);
        Qhat = Qhat_Cluster{i}(ind_fit,k);
        sigma_Qhat = sigma_Qhat_Cluster{i}(ind_fit,k);
        
        [a, C, sigma_a, sigma_C, Qhat_fit_values, sigma_Qhat_fit_values] = linearfit(tau_Qhat_fit{i}{k},Qhat,sigma_Qhat);

        tauth = -a/C;
        sigma_tauth = -sigma_a/C;
        C_Qhat_tau_fit{i}(k) = C;
        tauth_Qhat_tau_fit{i}(k) = tauth;
        sigma_C_Qhat_tau_fit{i}(k) = sigma_C;
        sigma_tauth_Qhat_tau_fit{i}(k) = sigma_tauth;
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
tau_zqi_fit = cell(N_Cluster,1);
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
    tau_zqi_fit{i} = cell(N_bins,1);
    zqi_fit{i} = cell(N_bins,1);  
    sigma_zqi_fit{i} = cell(N_bins,1);
    zqi_bar_Cluster{i} = zeros(N_bins,1); %mean zqi
    zqi_sigma_Cluster{i} = zeros(N_bins,1); %std dev of zqi
    
    %go through bins
    for k = 1:N_bins
        
        ind_fit = find(~isnan(zqi_Cluster{i}(:,k))); %get indices of non-NaN values for fitting
        tau_zqi_fit{i}{k} = tau_profile_Cluster{i}(ind_fit);
        zqi = zqi_Cluster{i}(ind_fit,k);
        sigma_zqi = sigma_zqi_Cluster{i}(ind_fit,k);
        
        [a, b, sigma_a, sigma_b, zqi_fit_values, sigma_zqi_fit_values] = linearfit(tau_zqi_fit{i}{k},zqi,sigma_zqi);

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
a_zqinorm_taunorm_fit = cell(N_Cluster,1); %fitting intercept
b_zqinorm_taunorm_fit = cell(N_Cluster,1); %fitting slope
sigma_a_zqinorm_taunorm_fit = cell(N_Cluster,1);
sigma_b_zqinorm_taunorm_fit = cell(N_Cluster,1);
taunorm_zqinorm_fit = cell(N_Cluster,1);
zqinorm_fit = cell(N_Cluster,1);
sigma_zqinorm_fit = cell(N_Cluster,1);
zqinorm_bar_Cluster = cell(N_Cluster,1); %mean zqi/d
zqinorm_sigma_Cluster = cell(N_Cluster,1); %std dev of zqi/d

for i = 1:N_Cluster
    
    %inialize fit values
    a_zqinorm_taunorm_fit{i} = zeros(N_bins,1);
    b_zqinorm_taunorm_fit{i} = zeros(N_bins,1);
    sigma_a_zqinorm_taunorm_fit{i} = zeros(N_bins,1);
    sigma_b_zqinorm_taunorm_fit{i} = zeros(N_bins,1);
    taunorm_zqinorm_fit{i} = cell(N_bins,1);
    zqinorm_fit{i} = cell(N_bins,1);  
    sigma_zqinorm_fit{i} = cell(N_bins,1);
    zqinorm_bar_Cluster{i} = zeros(N_bins,1); %mean zqi
    zqinorm_sigma_Cluster{i} = zeros(N_bins,1); %std dev of zqi
    
    %go through bins
    for k = 1:N_bins
        ind_fit = find(~isnan(zqi_Cluster{i}(:,k))); %get indices of non-NaN values for fitting
        taunorm_zqinorm_fit{i}{k} = tau_profile_Cluster{i}(ind_fit)/tauth_Cluster(i);
        zqinorm = 1000*zqi_Cluster{i}(ind_fit,k)/d_bin_mid_Cluster(i,k);
        sigma_zqinorm = 1000*sigma_zqi_Cluster{i}(ind_fit,k)/d_bin_mid_Cluster(i,k);
        
        [a, b, sigma_a, sigma_b, zqinorm_fit_values, sigma_zqinorm_fit_values] = linearfit(taunorm_zqinorm_fit{i}{k},zqinorm,sigma_zqinorm);

        a_zqinorm_taunorm_fit{i}(k) = a;
        b_zqinorm_taunorm_fit{i}(k) = b;
        sigma_a_zqinorm_taunorm_fit{i}(k) = sigma_a;
        sigma_b_zqinorm_taunorm_fit{i}(k) = sigma_b;
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
a_fratio_tau_fit = cell(N_Cluster,1);
b_fratio_tau_fit = cell(N_Cluster,1);
sigma_a_fratio_tau_fit = cell(N_Cluster,1);
sigma_b_fratio_tau_fit = cell(N_Cluster,1);
tau_fratio_fit = cell(N_Cluster,1);
fratio_fit = cell(N_Cluster,1);
fratio_bar_Cluster = cell(N_Cluster,1);
fratio_sigma_Cluster = cell(N_Cluster,1);

for i = 1:N_Cluster
    
    %inialize fit values
    a_fratio_tau_fit{i} = zeros(N_bins,1);
    b_fratio_tau_fit{i} = zeros(N_bins,1);
    sigma_a_fratio_tau_fit{i} = zeros(N_bins,1);
    sigma_b_fratio_tau_fit{i} = zeros(N_bins,1);
    tau_fratio_fit{i} = cell(N_bins,1);
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

        a_fratio_tau_fit{i}(k) = a;
        b_fratio_tau_fit{i}(k) = b;
        sigma_a_fratio_tau_fit{i}(k) = sigma_a;
        sigma_b_fratio_tau_fit{i}(k) = sigma_b;
        
        [tau_fit, ind_sort] = sort(tau_fit);
        tau_fratio_fit{i}{k} = tau_fit;
        fratio_fit{i}{k} = fratio_fit_values(ind_sort);
        
        fratio_bar_Cluster{i}(k) = mean(fratio);
        fratio_sigma_Cluster{i}(k) = std(fratio);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT TO F_AIR/F_BED VS TAU/TAUTH %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize values
a_fratio_taunorm_fit = cell(N_Cluster,1);
b_fratio_taunorm_fit = cell(N_Cluster,1);
sigma_a_fratio_taunorm_fit = cell(N_Cluster,1);
sigma_b_fratio_taunorm_fit = cell(N_Cluster,1);
taunorm_fratio_fit = cell(N_Cluster,1);
fratio_fit = cell(N_Cluster,1);

for i = 1:N_Cluster
    
    %inialize fit values
    a_fratio_taunorm_fit{i} = zeros(N_bins,1);
    b_fratio_taunorm_fit{i} = zeros(N_bins,1);
    sigma_a_fratio_taunorm_fit{i} = zeros(N_bins,1);
    sigma_b_fratio_taunorm_fit{i} = zeros(N_bins,1);
    taunorm_fratio_fit{i} = cell(N_bins,1);
    fratio_fit{i} = cell(N_bins,1);  
    
    %go through bins
    for k = 1:N_bins
        fratio = f_profilebar_airborne_Cluster{i}(:,k)./f_surface_Cluster{i}(k);
        ind_fit = find(~isnan(fratio)); %get indices of non-NaN values for fitting
        taunorm_fit = tau_profile_Cluster{i}(ind_fit)/tauth_Cluster(i);
        fratio = fratio(ind_fit);
        
        [a, b, sigma_a, sigma_b, fratio_fit_values, sigma_fratio_fit_values] = linearfit(taunorm_fit,fratio);

        a_fratio_taunorm_fit{i}(k) = a;
        b_fratio_taunorm_fit{i}(k) = b;
        sigma_a_fratio_taunorm_fit{i}(k) = sigma_a;
        sigma_b_fratio_taunorm_fit{i}(k) = sigma_b;
        
        [taunorm_fit, ind_sort] = sort(taunorm_fit);
        taunorm_fratio_fit{i}{k} = taunorm_fit;
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
        ind_usable = ind_usable_profile_Cluster{i};
        ind_average = intersect(ind_ustnorm,ind_usable);
        dV_ustnormbar_airborne_Cluster{i}(j,:) = mean(dV_profilebar_airborne_Cluster{i}(ind_average,ind_d),1)...
            /sum(mean(dV_profilebar_airborne_Cluster{i}(ind_average,ind_d),1)); %normalize by volume fraction > d_min 
        dVdlogd_ustnormbar_airborne_Cluster{i}(j,:) = dV_ustnormbar_airborne_Cluster{i}(j,:)./dlogd_airborne_Cluster{i}(ind_d);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET ZQNORM ASSOCIATED WITH AIRBORNE D10, D50, AND D90 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zqinorm_d10_Cluster = cell(N_Cluster,1); %zqi/d associated with d10
sigma_zqinorm_d10_Cluster = cell(N_Cluster,1); %uncertainty in zqi/d associated with d10
zqinorm_d50_Cluster = cell(N_Cluster,1); %zqi/d associated with d50
sigma_zqinorm_d50_Cluster = cell(N_Cluster,1); %uncertainty in zqi/d associated with d50
zqinorm_d90_Cluster = cell(N_Cluster,1); %zqi/d associated with d90
sigma_zqinorm_d90_Cluster = cell(N_Cluster,1); %uncertainty in zqi/d associated with d90

for i = 1:N_Cluster
    %get zqinorm for d10
    d10_diff = abs(d_f_mid_Cluster{i}-d10_bar_airborne_Cluster(i));
    ind_d10 = find(d10_diff==min(d10_diff));
    zqinorm_d10_Cluster{i} = zqinorm_Cluster{i}(:,ind_d10);
    sigma_zqinorm_d10_Cluster{i} = sigma_zqinorm_Cluster{i}(:,ind_d10);
    
    %get zqinorm for d50
    d50_diff = abs(d_f_mid_Cluster{i}-d50_bar_airborne_Cluster(i));
    ind_d50 = find(d50_diff==min(d50_diff));
    zqinorm_d50_Cluster{i} = zqinorm_Cluster{i}(:,ind_d50);
    sigma_zqinorm_d50_Cluster{i} = sigma_zqinorm_Cluster{i}(:,ind_d50);
    
    %get zqinorm for d90
    d90_diff = abs(d_f_mid_Cluster{i}-d90_bar_airborne_Cluster(i));
    ind_d90 = find(d90_diff==min(d90_diff));
    zqinorm_d90_Cluster{i} = zqinorm_Cluster{i}(:,ind_d90);
    sigma_zqinorm_d90_Cluster{i} = sigma_zqinorm_Cluster{i}(:,ind_d90);
end

%%
%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%

%% save data
AnalysisData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis_Additional');
save(AnalysisData_Path,'*Cluster','Cluster*','*bins','binning_type',...
    'dref_type','*min','*max',...
    '*fine','*medium','*coarse','*fit');