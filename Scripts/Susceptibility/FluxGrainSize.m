%% initialize
clearvars;

%%
%%%%%%%%%%%%%%
% PARAMETERS %
%%%%%%%%%%%%%%

%% physical parameters
rho_s = 2650*1e3; %sediment density, g/m^3
A_w = 30*0.6*(1e-3)^2; %area of Wenglor, m^2

%% partial flux grain-size fixed bin values
d_bin_lower = [0.1, 0.18, 0.25, 0.3, 0.4, 0.5, 0.7]; %lower edge of bin (mm)
d_bin_upper = [0.18, 0.25, 0.3, 0.4, 0.5, 0.7, 1.0]; %upper edge of bin (mm)
d_bin_mid = geomean([d_bin_lower; d_bin_upper]); %midpoint of bin (mm)
N_d_bin = length(d_bin_lower); %number of bins

%% partial flux bins in terms of CDF values
f_bin_lower = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85]; %CDF values for size bins - lower limit
f_bin_upper = [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]; %CDF values for size bins - upper limit
f_bin_mid = mean([f_bin_lower; f_bin_upper]); %CDF values for size bins - midpoint
N_f_bin = length(f_bin_lower); %number of bins
f_d_surface_CDF = f_bin_upper - f_bin_lower; %fractions of distribution for CDF values

%% bins in terms of values with respect to d50 (i.e., d/d50) -- not currently using
dhat_bin_lower = [0.2, 0.4, 0.6, 0.8, 0.95, 1.05, 1.2, 1.4, 1.6]; %d/d50 values for size bins - lower limit
dhat_bin_upper = [0.4, 0.6, 0.8, 0.95, 1.05, 1.2, 1.4, 1.6, 1.8]; %d/d50 values for size bins - upper limit
dhat_bin_mid = mean([dhat_bin_lower; dhat_bin_upper]); %d/d50 values for size bins - midpoint
N_dhat_bin = length(dhat_bin_lower); %number of bins

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

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%%
%%%%%%%%%%%%%%%%%
% PLOTTING INFO %
%%%%%%%%%%%%%%%%%

%% plotting information
PlotFont = 10; %font for labels
LineWidth_Plot = 1; %width of lines
Marker_Site = {'s','d','o'}; %markers for Sites
Color_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]}; %colors for sites
Marker_bin = {'s','d','o','v','^','<','>','p','h'}; %markers for bins
Color_bin = {[0 0 1],[0.1 0 0.85],[0.25 0 0.7],[0.4 0 0.55],[0.55 0 0.4],[0.7 0 0.25],[0.8 0 0.2],[0.9 0 0.1],[1 0 0]}; %colors for bins

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(GrainSizeData_Path); %grain size data
load(SaltationFluxData_Path); %saltation data
addpath(folder_Functions); %point MATLAB to location of functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - SURFACE SAMPLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize information for surface samples
time_surface = cell(N_Sites,1); %time of surface sample collection
f_d_surface = cell(N_Sites,1); %compute fraction of surface sample in each fixed bin
d_bin_lower_CDF = cell(N_Sites,1); %grain sizes associated with CDF values for size bins - lower limit
d_bin_upper_CDF = cell(N_Sites,1); %grain sizes associated with CDF values for size bins - upper limit
d_bin_mid_CDF = cell(N_Sites,1); %grain sizes associated with CDF values for size bins - midpoint

%% go through sites
for i = 1:N_Sites
    
    %% get data for site
    GrainSize_BSNE = GrainSizeData_all{i}.GrainSize_BSNE; %BSNE grain size data
    GrainSize_Surface = GrainSizeData_all{i}.GrainSize_Surface; %Surface grain size data
    Flux_BSNE = FluxBSNE_all{i}; %saltation flux for BSNEs 
    
    %% get information about data intervals for site
    N_FluxBSNE = length(Flux_BSNE); %number of BSNE flux intervals
    N_Surface = length(GrainSize_Surface); %number of surface grain size samples
    time_surface{i} = [GrainSize_Surface.CollectionTime]'; %time of surface sample collection
    
    %% initialize information for surface samples
    f_d_surface{i} = zeros(N_Surface,N_d_bin); %compute fraction of surface sample in each fixed bin
    d_bin_lower_CDF{i} = zeros(N_Surface,N_f_bin); %grain sizes associated with CDF values for size bins - lower limit
    d_bin_upper_CDF{i} = zeros(N_Surface,N_f_bin); %grain sizes associated with CDF values for size bins - upper limit
    d_bin_mid_CDF{i} = zeros(N_Surface,N_f_bin); %grain sizes associated with CDF values for size bins - midpoint
    
    %% go through each surface sample for calculations
    for j = 1:N_Surface
        %get grain size distribution information for this sample
        d_surface_lower = [GrainSize_Surface(j).gsd.Sizeclass_lower_mm]; %get diameters for lower edges of grain size bins
        d_surface_upper = [GrainSize_Surface(j).gsd.Sizeclass_upper_mm]; %get diameters for upper edges of grain size bins
        dV_surface = [GrainSize_Surface(j).gsd.retained]; %get volume fractions in bins

        %keep only 2nd to end values
        d_surface_lower = d_surface_lower(2:end);
        d_surface_upper = d_surface_upper(2:end);
        dV_surface = dV_surface(2:end)/sum(dV_surface(2:end));
        
        %compute fraction of grain size distribution in each bin for fixed bins
        for k = 1:N_d_bin
            ind_full_bin = find(d_surface_lower>=d_bin_lower(k) & d_surface_upper<=d_bin_upper(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_surface_upper(ind_bottom_bin))-log(d_bin_lower(k)))/...
                (log(d_surface_upper(ind_bottom_bin))-log(d_surface_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper(k))-log(d_surface_lower(ind_top_bin)))/...
                (log(d_surface_upper(ind_top_bin))-log(d_surface_lower(ind_top_bin))); %fraction of partial bin for this bin
            f_d_surface{i}(j,k) = sum(dV_surface(ind_full_bin))+...
                wt_bottom_bin*dV_surface(ind_bottom_bin)+...
                wt_top_bin*dV_surface(ind_top_bin); %compute fraction of grains by volume in this bin
        end
        
        %get reference grain sizes associated with CDF bins
        for k = 1:N_f_bin
            d_bin_lower_CDF{i}(j,k) = ReferenceGrainSizes_arbitrary(dV_surface, d_surface_lower, d_surface_upper, f_bin_lower(k)); %grain sizes associated with CDF values for size bins - lower limit
            d_bin_upper_CDF{i}(j,k) = ReferenceGrainSizes_arbitrary(dV_surface, d_surface_lower, d_surface_upper, f_bin_upper(k)); %grain sizes associated with CDF values for size bins - upper limit
        end
        d_bin_mid_CDF{i}(j,:) = geomean([d_bin_lower_CDF{i}(j,:); d_bin_upper_CDF{i}(j,:)]); %grain sizes associated with CDF values for size bins - midpoint
    end
end
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - AIRBORNE SAMPLES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize information for BSNE intervals - airborne samples
d_airborne_interval = cell(N_Sites,1); %grain sizes for airborne samples for interval
dlogd_airborne_interval = cell(N_Sites); %grain size interval for airborne samples for interval
Q_interval = cell(N_Sites,1); %total saltation flux for interval
ust_interval = cell(N_Sites,1); %shear velocity for interval
tau_interval = cell(N_Sites,1); %shear stress for interval
dV_airborne_interval_bar = cell(N_Sites,1); %volume fraction increment for airborne grain size for interval
dVdlogd_airborne_interval_bar = cell(N_Sites,1); %dV/dlogd for airborne grain size for interval
d10_airborne_interval_bar = cell(N_Sites,1); %10th percentile airborne grain size for interval
d50_airborne_interval_bar = cell(N_Sites,1); %50th percentile airborne grain size for interval
d90_airborne_interval_bar = cell(N_Sites,1); %90th percentile airborne grain size for interval

%% initialize information about grain size variation with height
d10_airborne_z = cell(N_Sites,1); %d10 - by height
d50_airborne_z = cell(N_Sites,1); %d50 - by height
d90_airborne_z = cell(N_Sites,1); %d90 - by height
dV_airborne_z = cell(N_Sites,1); %volume fraction increment - by height
z_airborne = cell(N_Sites,1); %height of grain size
z_airborne_ind = cell(N_Sites,1); %indices of intervals associated with z's
Cqn_airborne = cell(N_Sites,1); %expected calibration coefficient

%% go through sites
for i = 1:N_Sites
    
    %% computations for airborne samples
    %get size bins from first sample
    d_airborne_interval{i} = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_mid_mm];
    d_airborne_lower = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_lower_mm];
    d_airborne_upper = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_upper_mm];
    d_airborne_mid = [GrainSize_BSNE(1).gsd(2:end-1).Sizeclass_mid_mm];
    dlogd_airborne_interval{i} = log(d_airborne_upper) - log(d_airborne_lower);
    N_d = length(d_airborne_interval{i});
    
    %initialize matrices
    Q_interval{i} = zeros(N_FluxBSNE,1); %total fluxes for grain size samples
    ust_interval{i} = zeros(N_FluxBSNE,1); %shear velocities for grain size samples
    tau_interval{i} = zeros(N_FluxBSNE,1); %shear stress for interval
    dV_airborne_interval_bar{i} = zeros(N_FluxBSNE,N_d); %size distributions for grain size samples
    dVdlogd_airborne_interval_bar{i} = zeros(N_FluxBSNE,N_d); %size distributions for grain size samples
    d10_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d10 for grain size samples
    d50_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d50 for grain size samples
    d90_airborne_interval_bar{i} = zeros(N_FluxBSNE,1).*NaN; %weighted mean d90 for grain size samples
      
    %information about grain size variation with height
    d10_airborne_z{i} = [];
    d50_airborne_z{i} = [];
    d90_airborne_z{i} = [];
    dV_airborne_z{i} = [];
    z_airborne{i} = [];
    z_airborne_ind{i} = []; %indices of intervals associated with z's
    Cqn_airborne{i} = []; %expected calibration coefficient
            
    %perform airborne sample calculations
    for j = 1:N_FluxBSNE
        Q_interval{i}(j) = Flux_BSNE(j).Q.Q;
        z_BSNE = Flux_BSNE(j).z.z; %heights of BSNEs
        qz_BSNE = Flux_BSNE(j).qz.qz;
        Name_BSNE = Flux_BSNE(j).name;
        StartTime_BSNE = Flux_BSNE(j).StartTime; %get start time for BSNE interval
        EndTime_BSNE = Flux_BSNE(j).EndTime; %get end time for BSNE interval
        ind_GrainSize = find([GrainSize_BSNE.StartTime]<=StartTime_BSNE & [GrainSize_BSNE.EndTime]>StartTime_BSNE); %get indices within GrainSize_BSNE corresponding to time interval
        
        %get ust for interval
        ind_interval = find(StartTimes_all{i}>=StartTime_BSNE & EndTimes_all{i}<=EndTime_BSNE);
        ust_interval{i}(j)=mean(ustRe_all{i}(ind_interval)); %shear velocity for interval
        tau_interval{i}(j)=mean(tauRe_all{i}(ind_interval)); %shear stress for interval
        
        %initialize matrix of airborne samples
        N_airborne_interval = length(ind_GrainSize);
        dV_airborne_interval = zeros(N_airborne_interval,N_d);
        dVdlogd_airborne_interval = zeros(N_airborne_interval,N_d);
        qz_interval = zeros(N_airborne_interval,1); %partial fluxes for grain size samples
        
        %get each airborne size distribution
        for k = 1:N_airborne_interval
            dV_airborne_interval(k,:) = [GrainSize_BSNE(ind_GrainSize(k)).gsd(2:end-1).retained]/100; %by volume (divide by 100 to convert from % to fraction)
            dVdlogd_airborne_interval(k,:) = dV_airborne_interval(k,:)./dlogd_airborne_interval{i}; %by normalized volume
            try
                ind_BSNE = find(strcmp(Name_BSNE,GrainSize_BSNE(ind_GrainSize(k)).NameBSNE)); %get index for BSNE
                qz_interval(k) = qz_BSNE(ind_BSNE); %get flux associated with BSNE
            catch
                ind_BSNE = [];
                qz_interval(k) = 0;
            end
            
            %get reference sizes, add to list
            if ~isempty(ind_BSNE)
                [d10, d50, d90] = ReferenceGrainSizes(dV_airborne_interval(k,:), d_airborne_lower, d_airborne_upper);
                d10_airborne_z{i} = [d10_airborne_z{i}; d10];
                d50_airborne_z{i} = [d50_airborne_z{i}; d50];
                d90_airborne_z{i} = [d90_airborne_z{i}; d90];
                dV_airborne_z{i} = [dV_airborne_z{i}; dV_airborne_interval(k,:)];
                z_interval = z_BSNE(ind_BSNE); %get height associated with BSNE
                z_airborne{i} = [z_airborne{i}; z_interval];
                z_airborne_ind{i} = [z_airborne_ind{i}; j];
                d_bar = sum(dV_airborne_interval(k,:).*d_airborne_mid)/sum(dV_airborne_interval(k,:)); %get mean grain size, mm
                Cqn = (pi/6)*(1e-3*d_bar).^3*(rho_s)/A_w; %get expected calibration coefficient
                Cqn_airborne{i} = [Cqn_airborne{i}; Cqn]; %expected calibration coefficient
            end
        end
        
        %get mean of airborne sizes
        dV_airborne_interval_bar{i}(j,:) = sum(dV_airborne_interval.*qz_interval)/sum(qz_interval);
        dVdlogd_airborne_interval_bar{i}(j,:) = sum(dVdlogd_airborne_interval.*qz_interval)/sum(qz_interval);
                
        %compute benchmark grain sizes
        if(~isnan(dV_airborne_interval_bar{i}(j,1)))
            [d10, d50, d90] = ReferenceGrainSizes(dV_airborne_interval_bar{i}(j,:), d_airborne_lower, d_airborne_upper);
            d10_airborne_interval_bar{i}(j) = d10;
            d50_airborne_interval_bar{i}(j) = d50;
            d90_airborne_interval_bar{i}(j) = d90;
        end
    end
end

%% INFORMATION ABOUT FRACTIONS OF PARTICLES IN BINS

%% initialize fractions in bins - fixed bins
f_d_airborne = cell(N_Sites,1); %fraction of airborne sample in each bin
f_d_surface_airborne = cell(N_Sites,1); %fraction of surface sample in each bin - associated with airborne sample time interval
f_d_airborne_surface_ratio = cell(N_Sites,1); %fraction of airborne sample in size bin versus surface sample

%% initialize fractions in bins - CDF bins
f_d_airborne_CDF = cell(N_Sites,1); %compute fraction of airborne sample in each bin
f_d_airborne_surface_ratio_CDF = cell(N_Sites,1); %fraction of airborne sample in size bin versus surface sample

%% go through sites
for i = 1:N_Sites
    
    %initialize matrices
    f_d_airborne{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %compute fraction of airborne sample in each bin - fixed bins
    f_d_airborne_CDF{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %compute fraction of airborne sample in each bin - CDF bins
    
    %initialize list of volume fractions for surface samples - fixed bins
    f_d_surface_airborne{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %fraction of surface sample in each bin associated with airborne time interval
    f_d_airborne_surface_ratio{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %compute fraction of airborne sample in each bin versus surface fraction

    %initialize list of volume fractions for surface samples - CDF bins
    f_d_airborne_surface_ratio_CDF{i} = zeros(N_FluxBSNE,N_f_bin).*NaN; %compute fraction of airborne sample in each bin versus surface fraction
    
    %% go through BSNE intervals
    for j = 1:N_FluxBSNE

        %compute fraction of airborne sample in each bin - fixed bins
        for k = 1:N_d_bin
            ind_full_bin = find(d_airborne_lower>=d_bin_lower(k) & d_airborne_upper<=d_bin_upper(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_airborne_upper(ind_bottom_bin))-log(d_bin_lower(k)))/...
                (log(d_airborne_upper(ind_bottom_bin))-log(d_airborne_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper(k))-log(d_airborne_lower(ind_top_bin)))/...
                (log(d_airborne_upper(ind_top_bin))-log(d_airborne_lower(ind_top_bin))); %fraction of partial bin for this bin
            f_d_airborne{i}(j,k) = sum(dV_airborne_interval_bar{i}(j,ind_full_bin))+...
                wt_bottom_bin*dV_airborne_interval_bar{i}(j,ind_bottom_bin)+...
                wt_top_bin*dV_airborne_interval_bar{i}(j,ind_top_bin); %compute fraction of grains by volume in this bin
        end

        %compute fraction of airborne sample in each bin - CDF bins
        for k = 1:N_f_bin
            ind_full_bin = find(d_airborne_lower>=d_bin_lower_CDF{i}(k) & d_airborne_upper<=d_bin_upper_CDF{i}(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_airborne_upper(ind_bottom_bin))-log(d_bin_lower_CDF{i}(k)))/...
                (log(d_airborne_upper(ind_bottom_bin))-log(d_airborne_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper_CDF{i}(k))-log(d_airborne_lower(ind_top_bin)))/...
                (log(d_airborne_upper(ind_top_bin))-log(d_airborne_lower(ind_top_bin))); %fraction of partial bin for this bin
            f_d_airborne_CDF{i}(j,k) = sum(dV_airborne_interval_bar{i}(j,ind_full_bin))+...
                wt_bottom_bin*dV_airborne_interval_bar{i}(j,ind_bottom_bin)+...
                wt_top_bin*dV_airborne_interval_bar{i}(j,ind_top_bin); %compute fraction of grains by volume in this bin
        end
        
        %get associated surface sample information for airborne interval - fixed bins
        deltat_surface_airborne = abs(time_surface{i}-mean([StartTime_BSNE EndTime_BSNE])); %get difference between time of surface sample and BSNE sample      
        ind_Surface = find(deltat_surface_airborne == min(deltat_surface_airborne)); %get index of surface samples closest to airborne sample time
        if length(ind_Surface)==1
            f_d_surface_airborne{i}(j,:) = f_d_surface{i}(ind_Surface,:); %get closest sample      
        else
            f_d_surface_airborne{i}(j,:) = mean(f_d_surface{i}(ind_Surface,:)); %take mean of all closest samples        
        end
        f_d_airborne_surface_ratio{i}(j,:) = f_d_airborne{i}(j,:)./f_d_surface_airborne{i}(j,:); %ratio of fraction in airborne versus surface sample for size range
        
        %compute airborne surface ratio for CDF bins
        f_d_airborne_surface_ratio_CDF{i}(j,:) = f_d_airborne_CDF{i}(j,:)./f_d_surface_CDF; %ratio of fraction in airborne versus surface sample for size range
    end
end

%% SIZE-CONDITIONED BSNE FLUXES

%% initialize information about size-conditioned BSNE fluxes - fixed bins
q_profile_bin = cell(N_Sites,1); %size conditioned flux profile
Q_bin = cell(N_Sites,1); %size conditioned total flux
sigma_Q_bin = cell(N_Sites,1); %size conditioned total flux uncertainty
Qhat_bin = cell(N_Sites,1); %size conditioned total flux - normalized by surface fraction
sigma_Qhat_bin = cell(N_Sites,1); %size conditioned total flux uncertainty - normalized by surface fraction

%% initialize information about size-conditioned BSNE fluxes - CDF bins
q_profile_bin_CDF = cell(N_Sites,1); %size conditioned flux profile
Q_bin_CDF = cell(N_Sites,1); %size conditioned total flux
sigma_Q_bin_CDF = cell(N_Sites,1); %size conditioned total flux uncertainty
Qhat_bin_CDF = cell(N_Sites,1); %size conditioned total flux - normalized by surface fraction
sigma_Qhat_bin_CDF = cell(N_Sites,1); %size conditioned total flux uncertainty - normalized by surface fraction

%% go through sites
for i = 1:N_Sites
    
    %% information about size-conditioned BSNE fluxes - for fixed bins
    q_profile_bin{i} = cell(N_FluxBSNE,N_d_bin); %size conditioned flux profile
    Q_bin{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %size conditioned total flux
    sigma_Q_bin{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %size conditioned total flux uncertainty
    Qhat_bin{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %size conditioned total flux - normalized by surface fraction
    sigma_Qhat_bin{i} = zeros(N_FluxBSNE,N_d_bin).*NaN; %size conditioned total flux uncertainty - normalized by surface fraction

    %% information about size-conditioned BSNE fluxes - for CDF bins
    q_profile_bin_CDF{i} = cell(N_FluxBSNE,N_f_bin); %size conditioned flux profile
    Q_bin_CDF{i} = zeros(N_FluxBSNE,N_f_bin).*NaN; %size conditioned total flux
    sigma_Q_bin_CDF{i} = zeros(N_FluxBSNE,N_f_bin).*NaN; %size conditioned total flux uncertainty
    Qhat_bin_CDF{i} = zeros(N_FluxBSNE,N_f_bin).*NaN; %size conditioned total flux - normalized by surface fraction
    sigma_Qhat_bin_CDF{i} = zeros(N_FluxBSNE,N_f_bin).*NaN; %size conditioned total flux uncertainty - normalized by surface fraction
    
    %go through each BSNE time interval
    for j = 1:N_FluxBSNE
        
        %grain size information for this interval
        ind_interval = find(z_airborne_ind{i} == j);
        z_interval = z_airborne{i}(ind_interval);
        dV_interval = dV_airborne_z{i}(ind_interval,:);
        
        %BSNE flux profile for this interval
        z_BSNE = Flux_BSNE(j).z.z;
        sigma_z_BSNE = Flux_BSNE(j).z.sigma_z;
        qz_BSNE = Flux_BSNE(j).qz.qz;
        sigma_qz_BSNE = Flux_BSNE(j).qz.sigma;
        
        %compute size conditioned flux - fixed bins
        for k = 1:N_d_bin
            ind_full_bin = find(d_airborne_lower>=d_bin_lower(k) & d_airborne_upper<=d_bin_upper(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_airborne_upper(ind_bottom_bin))-log(d_bin_lower(k)))/...
                (log(d_airborne_upper(ind_bottom_bin))-log(d_airborne_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper(k))-log(d_airborne_lower(ind_top_bin)))/...
                (log(d_airborne_upper(ind_top_bin))-log(d_airborne_lower(ind_top_bin))); %fraction of partial bin for this bin
            
            %initialize q_profile_bin for specific BSNE interval
            q_profile_bin{i}{j,k} = zeros(size(qz_BSNE)); %profile of partial fluxes
            sigma_q_profile_bin = zeros(size(qz_BSNE)); %uncertainties on these values
            
            %go through each BSNE height, compute partial fluxes and their uncertainties
            for l = 1:length(z_BSNE)
                dz = abs(z_BSNE(l)-z_interval); %difference between BSNE height and heights for grain size information
                ind_closest = find(dz==min(dz)); %closest interval with grain size information
                if ~isempty(ind_closest)
                    f_d = sum(dV_interval(ind_closest,ind_full_bin))+...
                        wt_bottom_bin*dV_interval(ind_closest,ind_bottom_bin)+...
                        wt_top_bin*dV_interval(ind_closest,ind_top_bin); %compute fraction of grains by volume in this bin
                    q_profile_bin{i}{j,k}(l) = qz_BSNE(l)*f_d; %compute partial flux
                    sigma_q_profile_bin(l) = sigma_qz_BSNE(l)*f_d; %compute uncertainty in partial flux
                end
            end
            
            %fit profile
            [~,~,Q,~,~,sigma_Q] = qz_profilefit_exponential(q_profile_bin{i}{j,k}, z_BSNE, sigma_q_profile_bin, sigma_z_BSNE); %perform profile fitting
            Q_bin{i}(j,k) = Q; %size conditioned total flux
            sigma_Q_bin{i}(j,k) = sigma_Q; %size conditioned total flux uncertainty
            
            %normalize total flux by fraction of grains by volume in surface
            Qhat_bin{i}(j,k) = Q./f_d_surface_airborne{i}(j,k); %compute normalized total flux
            sigma_Qhat_bin{i}(j,k) = sigma_Q./f_d_surface_airborne{i}(j,k); %compute normalized total flux uncertainty
        end
        
        %compute size conditioned flux - CDF bins
        for k = 1:N_f_bin
            ind_full_bin = find(d_airborne_lower>=d_bin_lower_CDF{i}(k) & d_airborne_upper<=d_bin_upper_CDF{i}(k)); %grain size bins completely in bin
            ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
            wt_bottom_bin = (log(d_airborne_upper(ind_bottom_bin))-log(d_bin_lower_CDF{i}(k)))/...
                (log(d_airborne_upper(ind_bottom_bin))-log(d_airborne_lower(ind_bottom_bin))); %fraction of partial bin for this bin
            ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
            wt_top_bin = (log(d_bin_upper_CDF{i}(k))-log(d_airborne_lower(ind_top_bin)))/...
                (log(d_airborne_upper(ind_top_bin))-log(d_airborne_lower(ind_top_bin))); %fraction of partial bin for this bin
            
            %initialize q_profile_bin for specific BSNE interval
            q_profile_bin_CDF{i}{j,k} = zeros(size(qz_BSNE));
            sigma_q_profile_bin_CDF = zeros(size(qz_BSNE));
            
            %go through each BSNE height, compute partial fluxes and their uncertainties
            for l = 1:length(z_BSNE)
                dz = abs(z_BSNE(l)-z_interval); %difference between BSNE height and heights for grain size information
                ind_closest = find(dz==min(dz)); %closest interval with grain size information
                if ~isempty(ind_closest)
                    f_d = sum(dV_interval(ind_closest,ind_full_bin))+...
                        wt_bottom_bin*dV_interval(ind_closest,ind_bottom_bin)+...
                        wt_top_bin*dV_interval(ind_closest,ind_top_bin); %compute fraction of grains by volume in this bin
                    q_profile_bin_CDF{i}{j,k}(l) = qz_BSNE(l)*f_d; %compute partial flux
                    sigma_q_profile_bin_CDF(l) = sigma_qz_BSNE(l)*f_d; %compute uncertainty in partial flux
                end
            end
            
            %fit profile
            [~,~,Q,~,~,sigma_Q] = qz_profilefit_exponential(q_profile_bin_CDF{i}{j,k}, z_BSNE, sigma_q_profile_bin_CDF, sigma_z_BSNE); %peform profile fit
            Q_bin_CDF{i}(j,k) = Q; %size conditioned total flux
            sigma_Q_bin_CDF{i}(j,k) = sigma_Q; %size conditioned total flux uncertainty
            
            %normalize total flux by fraction of grains by volume in surface
            Qhat_bin_CDF{i}(j,k) = Q./f_d_surface_CDF(k); %compute normalized total flux
            sigma_Qhat_bin_CDF{i}(j,k) = sigma_Q./f_d_surface_CDF(k); %compute normalized total flux uncertainty
        end
    end    
end
% 
% %% perform fits to Qhat versus tau
% 
% % initialize values - fixed bins
% C_Qhat_tau = cell(N_Sites,1);
% tauth_Qhat_tau = cell(N_Sites,1);
% sigma_C_Qhat_tau = cell(N_Sites,1);
% sigma_tauth_Qhat_tau = cell(N_Sites,1);
% tau_fit = cell(N_Sites,1);
% Qhat_fit = cell(N_Sites,1);
% 
% % initialize values - CDF bins
% C_Qhat_tau_CDF = cell(N_Sites,1);
% tauth_Qhat_tau_CDF = cell(N_Sites,1);
% sigma_C_Qhat_tau_CDF = cell(N_Sites,1);
% sigma_tauth_Qhat_tau_CDF = cell(N_Sites,1);
% tau_fit_CDF = cell(N_Sites,1);
% Qhat_fit_CDF = cell(N_Sites,1);  
% 
% for i = 1:N_Sites
%     %inialize fit values - fixed bins
%     C_Qhat_tau{i} = zeros(N_d_bin,1);
%     tauth_Qhat_tau{i} = zeros(N_d_bin,1);
%     sigma_C_Qhat_tau{i} = zeros(N_d_bin,1);
%     sigma_tauth_Qhat_tau{i} = zeros(N_d_bin,1);
%     tau_fit{i} = cell(N_d_bin,1);
%     Qhat_fit{i} = cell(N_d_bin,1);  
%     
%     %fixed bins
%     for k = 1:N_d_bin
%         %remove NaN values prior to fitting
%         ind_isreal = find(~isnan(Qhat_bin{i}(:,k)));
%         tau_fit{i}{k} = tau_interval{i}(ind_isreal);
%         Qhat = Qhat_bin{i}(ind_isreal,k);
%         sigma_Qhat = sigma_Qhat_bin{i}(ind_isreal,k);
%         
%         %[a, C, sigma_a, sigma_C, Qhat_fit_values, sigma_Qhat_fit_values] = linearfit(tau_fit{i}{k},Qhat,sigma_Qhat);
%         [a, C, sigma_a, sigma_C, Qhat_fit_values, sigma_Qhat_fit_values] = linearfit(tau_fit{i}{k},Qhat);
% 
%         tauth = -a/C;
%         sigma_tauth = -sigma_a/C;
%         C_Qhat_tau{i}(k) = C;
%         tauth_Qhat_tau{i}(k) = tauth;
%         sigma_C_Qhat_tau{i}(k) = sigma_C;
%         sigma_tauth_Qhat_tau{i}(k) = sigma_tauth;
%         Qhat_fit{i}{k} = Qhat_fit_values;
%     end
% 
%     %inialize fit values - CDF bins
%     C_Qhat_tau_CDF{i} = zeros(N_f_bin,1);
%     tauth_Qhat_tau_CDF{i} = zeros(N_f_bin,1);
%     sigma_C_Qhat_tau_CDF{i} = zeros(N_f_bin,1);
%     sigma_tauth_Qhat_tau_CDF{i} = zeros(N_f_bin,1);
%     tau_fit_CDF{i} = cell(N_f_bin,1);
%     Qhat_fit_CDF{i} = cell(N_f_bin,1);  
%     
%     %CDF bins
%     for k = 1:N_f_bin
%         %remove NaN values prior to fitting
%         ind_isreal = find(~isnan(Qhat_bin_CDF{i}(:,k)));
%         tau_fit_CDF{i}{k} = tau_interval{i}(ind_isreal);
%         Qhat = Qhat_bin_CDF{i}(ind_isreal,k);
%         sigma_Qhat = sigma_Qhat_bin_CDF{i}(ind_isreal,k);
%         
%         %[a, C, sigma_a, sigma_C, Qhat_fit_values, sigma_Qhat_fit_values] = linearfit(tau_fit{i}{k},Qhat,sigma_Qhat);
%         [a, C, sigma_a, sigma_C, Qhat_fit_values, sigma_Qhat_fit_values] = linearfit(tau_fit_CDF{i}{k},Qhat);
% 
%         tauth = -a/C;
%         sigma_tauth = -sigma_a/C;
%         C_Qhat_tau_CDF{i}(k) = C;
%         tauth_Qhat_tau_CDF{i}(k) = tauth;
%         sigma_C_Qhat_tau_CDF{i}(k) = sigma_C;
%         sigma_tauth_Qhat_tau_CDF{i}(k) = sigma_tauth;
%         Qhat_fit_CDF{i}{k} = Qhat_fit_values;
%     end
% end
% 
% %% plot variation in reference grain sizes with saltation flux
% figure(1); clf;
% for i = 1:N_Sites
%     subplot(1,N_Sites,i); hold on;
%     h90 = plot(Q_interval{i},d90_airborne_interval_bar{i},'^');
%     h50 = plot(Q_interval{i},d50_airborne_interval_bar{i},'o');
%     h10 = plot(Q_interval{i},d10_airborne_interval_bar{i},'v');
%     c90 = get(h90,'Color');
%     c50 = get(h50,'Color');
%     c10 = get(h10,'Color');
%     plot([0 max(Q_interval{i})],d90_surface_site(i)*[1 1],'Color',c90);
%     plot([0 max(Q_interval{i})],d50_surface_site(i)*[1 1],'Color',c50);
%     plot([0 max(Q_interval{i})],d10_surface_site(i)*[1 1],'Color',c10);
%     ylim([0 0.9]);
%     xlabel('saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','Latex')
%     ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     title(Sites{i});
% end
% legend('d_{90}','d_{50}','d_{10}');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
% print([folder_Plots,'FluxGrainSize.png'],'-dpng');
% 
% %%plot variation in reference grain sizes with shear velocity
% figure(2); clf;
% for i = 1:N_Sites
%     subplot(1,N_Sites,i); hold on;
%     h90 = plot(ust_interval{i},d90_airborne_interval_bar{i},'^');
%     h50 = plot(ust_interval{i},d50_airborne_interval_bar{i},'o');
%     h10 = plot(ust_interval{i},d10_airborne_interval_bar{i},'v');
%     c90 = get(h90,'Color');
%     c50 = get(h50,'Color');
%     c10 = get(h10,'Color');
%     plot([0 max(ust_interval{i})],d90_surface_site(i)*[1 1],'Color',c90);
%     plot([0 max(ust_interval{i})],d50_surface_site(i)*[1 1],'Color',c50);
%     plot([0 max(ust_interval{i})],d10_surface_site(i)*[1 1],'Color',c10);
%     ylim([0 0.9]);
%     xlabel('shear velocity, $$u_{*}$$ (m s$$^{-1}$$)','Interpreter','Latex')
%     ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     title(Sites{i});
% end
% legend('d_{90}','d_{50}','d_{10}');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
% print([folder_Plots,'ShearVelocityGrainSize.png'],'-dpng');
% 
% %%plot variation in partial fluxes with total flux
% figure(3); clf;
% for k = 1:N_d_bin
%     subplot(1,N_d_bin,k); hold on;
%     for i = 1:N_Sites
%         plot(Q_interval{i},f_d_airborne_surface_ratio{i}(:,k),Marker_Site{i},'Color',Color_Site{i})
%         title(['d = ',num2str(d_bin_lower(k)),' - ',num2str(d_bin_upper(k)),' mm']);
%     end
%     xlabel('saltation flux, $$Q$$','interpreter','latex');
%     ylabel('fraction of size in air versus surface, $$f_{d,air}/f_{d,surf}$$','interpreter','latex');
%     set(gca,'YScale','Log','YTickLabel',{'0.01','0.1','1','10'},'Box','On');
%     ylim([1e-2 1e1]);
% end
% legend(Sites,'Location','NorthEast');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 4]);
% print([folder_Plots,'Flux_GrainSizeRelativeFlux.png'],'-dpng');
% 
% %%plot variation in partial fluxes with shear velocity
% figure(4); clf;
% for k = 1:N_d_bin
%     subplot(1,N_d_bin,k); hold on;
%     for i = 1:N_Sites
%         plot(ust_interval{i},f_d_airborne_surface_ratio{i}(:,k),Marker_Site{i},'Color',Color_Site{i})
%         title(['d = ',num2str(d_bin_lower(k)),' - ',num2str(d_bin_upper(k)),' mm']);
%     end
%     xlabel('shear velocity, $$u_{*}$$ (m s$$^{-1}$$)','Interpreter','Latex')
%     ylabel('fraction of size in air versus surface, $$f_{d,air}/f_{d,surf}$$','interpreter','latex');
%     set(gca,'YScale','Log','YTickLabel',{'0.01','0.1','1','10'},'Box','On');
%     ylim([1e-2 1e1]);
% end
% legend(Sites,'Location','NorthEast');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 4]);
% print([folder_Plots,'ShearVelocity_GrainSizeRelativeFlux.png'],'-dpng');
% 
% %%plot variation in airborne grain sizes with height
% figure(5); clf;
% for i = 1:N_Sites
%     %subplot(1,N_Sites,i); hold on;
%     subplot('Position',[-0.26+0.325*i 0.13 0.27 0.82]); hold on;
%     h90 = plot(z_airborne{i},d90_airborne_z{i},'^');
%     h50 = plot(z_airborne{i},d50_airborne_z{i},'o');
%     h10 = plot(z_airborne{i},d10_airborne_z{i},'v');
%     c90 = get(h90,'Color');
%     c50 = get(h50,'Color');
%     c10 = get(h10,'Color');
%     s90 = plot([0 max(z_airborne{i})],d90_surface_site(i)*[1 1],'Color',c90);
%     s50 = plot([0 max(z_airborne{i})],d50_surface_site(i)*[1 1],'Color',c50);
%     s10 = plot([0 max(z_airborne{i})],d10_surface_site(i)*[1 1],'Color',c10);
%     ylim([0 0.9]);
%     xlabel('height, $$z$$ (m)','Interpreter','Latex')
%     if i==1
%         ylabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
%     elseif i==2
%         legend([h90,h50,h10],'airborne d_{90}','airborne d_{50}','airborne d_{10}','Location','NorthEast');
%     elseif i==3
%         legend([s90,s50,s10],'surface d_{90}','surface d_{50}','surface d_{10}','Location','NorthEast');
%     end
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On','FontSize',PlotFont);
%     title(Sites{i});
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9 4],'PaperPosition',[0 0 9 4],'PaperPositionMode','Manual');
% print([folder_Plots,'HeightGrainSize.png'],'-dpng');
% 
% %%plot variation in expected calibration factor with height
% figure(6); clf; 
% subplot('Position',[0.08 0.115 0.4 0.86]); hold on;
% for i = 1:N_Sites
%     for j = 1:length(zW_all{i})
%         if ~isnan(zW_all{i}{j})
%             plot(zW_all{i}{j},Cqnbar_all{i}{j},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
%         end
%     end
% end
% xlabel('Wenglor height, $$z$$ (m)','Interpreter','Latex')
% ylabel('Observed calibration factor, $$C_{qn}$$ (g m$$^{-2}$$)','Interpreter','Latex');
% text(0.012,80,'(a)','FontSize',PlotFont);
% xlim([1e-2, 1]); 
% ylim([1e-1 1e2]);
% set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont);
% 
% subplot('Position',[0.58 0.115 0.4 0.86]); hold on;
% for i = 1:N_Sites
%     plot(z_airborne{i},Cqn_airborne{i},Marker_Site{i},'Color',Color_Site{i},'LineWidth',LineWidth_Plot);
% end
% xlabel('Height of grain size measurement, $$z$$ (m)','Interpreter','Latex')
% ylabel('Expected calibration factor, $$C_{qn,pred}$$ (g m$$^{-2}$$)','Interpreter','Latex');
% text(0.012,80,'(b)','FontSize',PlotFont);
% legend(Sites,'Location','SouthEast');
% xlim([1e-2, 1]); 
% ylim([1e-1 1e2]);
% set(gca,'xscale','log','yscale','log','box','on','FontSize',PlotFont)
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperSize',[8 4.5],'PaperPosition',[0 0 8 4.5],'PaperPositionMode','Manual');
% print([folder_Plots,'HeightCalibration.png'],'-dpng');
% 
% %% PLOT size-conditioned Q VS u*
% figure(7); clf;
% 
% 
% %create legend
% legend_items = cell(N_d_bin,1);
% for k = 1:N_d_bin
%     legend_items{k} = ['d_',int2str(k),'=',num2str(d_bin_lower(k)),'-',num2str(d_bin_upper(k)),' mm'];
% end
% 
% for i = 1:N_Sites
%     %subplot(1,N_Sites,i); hold on;
%     subplot('Position',[-0.24+0.32*i 0.13 0.27 0.82]); hold on;
% 
%     for k = 1:N_d_bin
%         %plot data
%         plot(tau_interval{i},Q_bin{i}(:,k),Marker_bin{k},'Color',Color_bin{k});
%     end
%     
% 
%     
% %     %plot error bars
% %     for j = 1:length(tauRe_all{i})
% %         plot(ones(2,1)*tauRe_all{i}(j),Q_all{i}(j)+[-1 1]*sigma_Q_all{i}(j),'Color',Color_Site{i},'LineWidth',LineWidth_plot); %y error
% %         plot(tauRe_all{i}(j)+[1 -1]*sigma_tauRe_all{i}(j),ones(2,1)*Q_all{i}(j),'Color',Color_Site{i},'LineWidth',LineWidth_plot); %x error
% %     end
%     
%     %format plot
%     %xlim([0 0.45]);
%     xlim([0 ceil(max(tau_interval{i}/0.05))*0.05]);
%     %ylim([0 65]);
%     xlims = xlim;
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
%     if i==1
%         ylabel('Size-specific mass flux, $$Q_{i}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
%         legend(legend_items,'Location','NorthWest');
%     end
%     %text(0.05*xlims(2), 62, panel_labels{i},'FontSize',9,'FontWeight','Bold');
%     title(SiteNames{i});
%     set(gca,'FontSize',PlotFont);
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'Q_d.png'],'-dpng'); %for draft
% %print([folder_Plots,'Q_d.tif'],'-dtiff','-r600'); %for publication
% 
% 
% %% PLOT size-conditioned normalized Q VS u* - fixed bins
% figure(8); clf;
% 
% %create legend
% legend_items = cell(N_d_bin,1);
% for k = 1:N_d_bin
%     legend_items{k} = ['d_',int2str(k),'=',num2str(d_bin_lower(k)),'-',num2str(d_bin_upper(k)),' mm'];
% end
% 
% for i = 1:N_Sites
%     %subplot(1,N_Sites,i); hold on;
%     subplot('Position',[-0.24+0.32*i 0.13 0.27 0.82]); hold on;
% 
%     for k = 1:N_d_bin
%         %plot data
%         plot(tau_interval{i},Qhat_bin{i}(:,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3); %plot values
%         %plot(tau_interval{i},Qhat_bin{i}(:,k),'o');
%     end
%     
%     %plot fit
%     for k = 1:N_d_bin
%         plot([tauth_Qhat_tau{i}(k); tau_fit{i}{k}],[0; Qhat_fit{i}{k}],'Color',Color_bin{k}); %plot fit
%     end   
%     
%     %format plot
%     %xlim([0 0.45]);
%     xlim([0 ceil(max(tau_interval{i}/0.05))*0.05]);
%     %ylim([0 65]);
%     xlims = xlim;
%     ylims = ylim;
%     ylim([0 ylims(2)]);
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
%     if i==1
%         ylabel('Bed fraction norm. size-spec. flux, $$\hat{Q_{i}}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
% %     elseif i==2
%         h_legend = legend(legend_items,'Location','NorthWest');
%         set(h_legend,'FontSize',6);
%     end
%     %text(0.05*xlims(2), 62, panel_labels{i},'FontSize',9,'FontWeight','Bold');
%     title(SiteNames{i});
%     set(gca,'FontSize',PlotFont);
%     %set(gca,'yscale','log');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'Qhat_d.png'],'-dpng'); %for draft
% %print([folder_Plots,'Qhat_d_30minute.tif'],'-dtiff','-r600'); %for publication
% 
% 
% %% PLOT size-conditioned normalized Q VS u* - CDF bins
% figure(9); clf;
% 
% %create legend
% legend_items = cell(N_d_bin,1);
% for k = 1:N_d_bin
%     legend_items{k} = ['d_',int2str(k),'=',num2str(100*f_bin_lower(k),2),'-',num2str(100*f_bin_upper(k),2),' pctile'];
% end
% 
% for i = 1:N_Sites
%     %subplot(1,N_Sites,i); hold on;
%     subplot('Position',[-0.24+0.32*i 0.13 0.27 0.82]); hold on;
% 
%     for k = 1:N_f_bin
%         %plot data
%         plot(tau_interval{i},Qhat_bin_CDF{i}(:,k),Marker_bin{k},'Color',Color_bin{k},'MarkerSize',3);
%         %plot(tau_interval{i},Qhat_bin_CDF{i}(:,k),'o');
%     end
%     
%     %plot fit
%     for k = 1:N_f_bin
%         plot([tauth_Qhat_tau_CDF{i}(k); tau_fit_CDF{i}{k}],[0; Qhat_fit_CDF{i}{k}],'Color',Color_bin{k}); %plot fit
%     end   
%     
%     %format plot
%     %xlim([0 0.45]);
%     xlim([0 ceil(max(tau_interval{i}/0.05))*0.05]);
%     %ylim([0 65]);
%     xlims = xlim;
%     ylims = ylim;
%     ylim([0 ylims(2)]);
%     set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
%     xlabel('Wind shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
%     if i==1
%         ylabel('Bed fraction norm. size-spec. flux, $$\hat{Q_{i}}$$ (gm$$^{-1}$$s$$^{-1}$$)','Interpreter','Latex');
% %     elseif i==2
%         h_legend = legend(legend_items,'Location','NorthWest');
%         set(h_legend,'FontSize',6);
%     end
%     %text(0.05*xlims(2), 62, panel_labels{i},'FontSize',9,'FontWeight','Bold');
%     title(SiteNames{i});
%     set(gca,'FontSize',PlotFont);
%     %set(gca,'yscale','log');
% end
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'Qhat_d_CDF.png'],'-dpng'); %for draft
% %print([folder_Plots,'Qhat_d_30minute.tif'],'-dtiff','-r600'); %for publication
% 
% 
% %% PLOT fit values - fixed bins
% figure(10); clf;
% 
% %plotting info
% Marker_Site = {'s','d','o'};
% Color_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
% 
% %plot C values and error bars
% subplot(1,2,1); hold on;
% for i = 1:N_Sites
%     plot(d_bin_mid,C_Qhat_tau{i},Marker_Site{i},'Color',Color_Site{i}); %values
% end
% for i = 1:N_Sites
%     for k = 1:N_d_bin
%         plot(d_bin_mid(k)*[1 1],C_Qhat_tau{i}(k)+[-1 1]*sigma_C_Qhat_tau{i}(k),'Color',Color_Site{i}); %error bars
%     end
% end
% xlabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
% ylabel('flux coefficient, $$C$$ (s)','Interpreter','Latex');
% legend(SiteNames);
% set(gca,'Yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
% 
% %plot tauth values
% subplot(1,2,2); hold on;
% for i = 1:N_Sites
%     plot(d_bin_mid,tauth_Qhat_tau{i},Marker_Site{i},'Color',Color_Site{i});
% end
% for i = 1:N_Sites
%     for k = 1:N_d_bin
%         plot(d_bin_mid(k)*[1 1],tauth_Qhat_tau{i}(k)+[-1 1]*sigma_tauth_Qhat_tau{i}(k),'Color',Color_Site{i}); %error bars
%     end
% end
% xlabel('grain diameter, $$d$$ (mm)','Interpreter','Latex');
% ylabel('threshold, $$\tau_{th}$$ (Pa)','Interpreter','Latex');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% legend(SiteNames,'Location','SouthWest');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'C_tauth_d.png'],'-dpng'); %for draft
% 
% %% PLOT fit values - CDF bins
% figure(11); clf;
% 
% Color_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
% 
% %plot C values and error bars
% subplot(1,2,1); hold on;
% for i = 1:N_Sites
%     plot(100*f_bin_mid,C_Qhat_tau_CDF{i},Marker_Site{i},'Color',Color_Site{i}); %values
% end
% for i = 1:N_Sites
%     for k = 1:N_d_bin
%         plot(100*f_bin_mid(k)*[1 1],C_Qhat_tau_CDF{i}(k)+[-1 1]*sigma_C_Qhat_tau_CDF{i}(k),'Color',Color_Site{i}); %error bars
%     end
% end
% xlabel('grain size percentile');
% ylabel('flux coefficient, $$C$$ (s)','Interpreter','Latex');
% legend(SiteNames,'Location','SouthEast');
% set(gca,'Yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
% 
% %plot tauth values
% subplot(1,2,2); hold on;
% for i = 1:N_Sites
%     plot(100*f_bin_mid,tauth_Qhat_tau_CDF{i},Marker_Site{i},'Color',Color_Site{i});
% end
% for i = 1:N_Sites
%     for k = 1:N_d_bin
%         plot(100*f_bin_mid(k)*[1 1],tauth_Qhat_tau_CDF{i}(k)+[-1 1]*sigma_tauth_Qhat_tau_CDF{i}(k),'Color',Color_Site{i}); %error bars
%     end
% end
% xlabel('grain size percentile');
% ylabel('threshold, $$\tau_{th}$$ (Pa)','Interpreter','Latex');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% legend(SiteNames,'Location','SouthWest');
% 
% %print plot
% set(gcf,'PaperUnits','inches','PaperSize',[9.5 4],'PaperPosition',[0 0 9.5 4],'PaperPositionMode','Manual');
% print([folder_Plots,'C_tauth_d_CDF.png'],'-dpng'); %for draft