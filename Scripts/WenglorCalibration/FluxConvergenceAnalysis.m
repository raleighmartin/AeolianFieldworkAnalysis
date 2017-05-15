%% INITIALIZATION
%initialize
clearvars;
close all;

%% PARAMETERS AND INPUTS
zq_max = 0.5; %maximum value for zq (m) - otherwise, it's an outlier
zq_Site = [0.097, 0.107, 0.055]; %saltation height at sites from flux law paper

%% LOAD DATA AND FUNCTIONS
%folders for loading data, saving data, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Functions = '../Functions/'; %folder with functions

%paths for loading and saving data
LoadData_Path = strcat(folder_LoadData,'DataFullSubwindowCalcs_30min_Restricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'DataFullSubwindowAnalysis_30min_Restricted'); %path for 30 minute data

%load data
load(LoadData_Path); %load window data

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%get time information
T_subwindow_s = seconds(T_subwindow); %duration of subwindows

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - PROFILE FITTING %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute fraction of Wenglor heights used in analysis
f_zW_analysis = cell(N_Sites,1); %fraction of Wenglors used
for i = 1:N_Sites
    f_zW_analysis{i} = cell(N_T_subwindow,1);
    for m = 1:N_T_subwindow      
        f_zW_analysis{i}{m} = N_z_fit_subwindow{i}{m}./N_z_subwindow{i}{m};
    end
end

%% get indices for subsets of data
ind_real = cell(N_Sites,1); %indices of profile fits that are not NaN
ind_positive_zq = cell(N_Sites,1); %indices of profile fits that give positive zq
ind_nooutlier_zq = cell(N_Sites,1); %indices of profile fits that do not give outlier zq
ind_positive_Qsum = cell(N_Sites,1); %indices of profiles that give positive Qsum
ind_full = cell(N_Sites,1); %indices of profile fits for which all Wenglor heights are represented
for i = 1:N_Sites
    ind_real{i} = cell(N_T_subwindow,1);
    ind_positive_zq{i} = cell(N_T_subwindow,1);
    ind_nooutlier_zq{i} = cell(N_T_subwindow,1);
    ind_positive_Qsum{i} = cell(N_T_subwindow,1);
    ind_full{i} = cell(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_real{i}{m} = find(~isnan(Qfit_subwindow{i}{m}));
        ind_positive_zq{i}{m} = find(zq_subwindow{i}{m}>0);
        ind_nooutlier_zq{i}{m} = find(zq_subwindow{i}{m}>0 & zq_subwindow{i}{m}<=zq_max);
        ind_positive_Qsum{i}{m} = find(Qsum_subwindow{i}{m}>0);
        ind_full{i}{m} = find(f_zW_analysis{i}{m}==1);
    end
end

%% compute median fraction of Wenglor heights used in profiles (and fraction with 100% used)
median_f_zW_analysis = cell(N_Sites,1); %median fraction of Wenglors used
lowerquartile_f_zW_analysis = cell(N_Sites,1); %25th percentile of Wenglors used
upperquartile_f_zW_analysis = cell(N_Sites,1); %75th percentile of Wenglors used
fraction_f_zW_full_analysis = cell(N_Sites,1); %fraction of time intervals with 100% usage
for i = 1:N_Sites
    median_f_zW_analysis{i} = zeros(N_T_subwindow,1);
    lowerquartile_f_zW_analysis{i} = zeros(N_T_subwindow,1); %25th percentile of Wenglors used
    upperquartile_f_zW_analysis{i} = zeros(N_T_subwindow,1); %75th percentile of Wenglors used
    fraction_f_zW_full_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        f_zW = f_zW_analysis{i}{m}(ind_positive_Qsum{i}{m}); %use only intervals with positive Qsum
        median_f_zW_analysis{i}(m) = median(f_zW);
        lowerquartile_f_zW_analysis{i}(m) = prctile(f_zW,25);
        upperquartile_f_zW_analysis{i}(m) = prctile(f_zW,75);
        fraction_f_zW_full_analysis{i}(m) = length(find(f_zW==1))/length(f_zW);
    end
end

%% compute fraction of positive flux intervals with outlier values for q(z) fits
f_outlier_analysis = cell(N_Sites,1); %fraction of profile fits that are well defined but have outliers
f_full_outlier_analysis = cell(N_Sites,1); %fraction of full profile fits that are well defined but have outliers
for i = 1:N_Sites
    f_outlier_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        length_positive_Qsum = length(ind_positive_Qsum{i}{m}); %number of positive Qsum values
        length_notoutlier_positive_Qsum = length(intersect(ind_positive_Qsum{i}{m},ind_nooutlier_zq{i}{m})); %number of positive Qsum values that are not outlier zq's
        length_full_positive_Qsum = length(intersect(ind_positive_Qsum{i}{m},ind_full{i}{m})); %number of positive Qsum values that have full q(z) profiles
        length_full_notoutlier_positive_Qsum = length(intersect(intersect(ind_positive_Qsum{i}{m},ind_full{i}{m}),ind_nooutlier_zq{i}{m})); %number of positive Qsum values that have full q(z) profiles and are not outlier zq's
        f_outlier = (length_positive_Qsum - length_notoutlier_positive_Qsum) / length_positive_Qsum; %compute fraction of positive Qsum that are not outliers for fits
        f_outlier_analysis{i}(m) = f_outlier; %add to list
        f_full_outlier = (length_full_positive_Qsum - length_full_notoutlier_positive_Qsum) / length_full_notoutlier_positive_Qsum; %compute fraction of full q(z) profiles that are not outliers for fits
        f_full_outlier_analysis{i}(m) = f_full_outlier; %add to list
    end
end

%% compute median chi2nu for q(z) profile fits
Chi2nu_analysis = cell(N_Sites,1);
median_Chi2nu_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    Chi2nu_analysis{i} = cell(N_T_subwindow,1);
    median_Chi2nu_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_calc = ind_positive_zq{i}{m}; %use only intervals with positive zq for fit (i.e., positive flux)
        Chi2nu = Chi2_Qfit_subwindow{i}{m}(ind_calc)./df_Qfit_subwindow{i}{m}(ind_calc); %compute chi2nu
        Chi2nu_analysis{i}{m} = Chi2nu; %add these to list
        median_Chi2nu_analysis{i}(m) = median(Chi2nu); %compute median and add to list
    end
end

%% subanalysis of 0 values in profiles

%set parameters
T_subwindow_s_f0 = 5;
zW_zq_f0 = 2;
u_min = 5.5:1:10.5;
u_max = 6.5:1:11.5;
u_f0_analysis = mean([u_min;u_max]);
N_u = length(u_min);
ind_T = find(T_subwindow_s == T_subwindow_s_f0);

%initialize values
f0_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    
    %go through subwindows, find q's of nearest z's
    N_subwindow = length(qbar_subwindow{i}{ind_T});
    q_zW_zq = zeros(N_subwindow,1);
    for j = 1:N_subwindow
        zW_zq_diff = abs(zW_subwindow{i}{ind_T}{j}/zq_Site(i)-zW_zq_f0);
        ind_z = find(zW_zq_diff==min(zW_zq_diff),1);
        q_zW_zq(j) = qbar_subwindow{i}{ind_T}{j}(ind_z);
    end
    
    %examine only intervals for which Qsum is positive
    ind_Qsum = ind_positive_Qsum{i}{ind_T};
    u = ubar_subwindow{i}{ind_T}(ind_Qsum);
    q = q_zW_zq(ind_Qsum);
    
    %compute f0 for subclasses of u
    f0_analysis{i} = zeros(N_u,1);
    for j = 1:N_u
        ind_u = intersect(find(u>=u_min(j)),find(u<u_max(j)));
        f0_analysis{i}(j) = length(find(q(ind_u)==0))/length(ind_u);
    end
end

%%
% %%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - WIND %
% %%%%%%%%%%%%%%%%%%%%%
ubar_all_analysis = ubar_subwindow; %get all ubar values
ubar_real_analysis = cell(N_Sites,1); %values of u for all real profiles
ubar_full_analysis = cell(N_Sites,1); %values for u for all full profiles
ubar_full_notoutlier_analysis = cell(N_Sites,1); %values of u for all non-outlier profiles
for i = 1:N_Sites
    ubar_real_analysis{i} = cell(N_T_subwindow,1); %values for all real profiles
    ubar_full_analysis{i} = cell(N_T_subwindow,1); %values for u for all full profiles
    ubar_full_notoutlier_analysis{i} = cell(N_T_subwindow,1); %values of u for all non-outlier profiles
    for m = 1:N_T_subwindow
        ubar_real_analysis{i}{m} = ubar_subwindow{i}{m}(ind_real{i}{m}); %values for u for all real profiles
        ubar_full_analysis{i}{m} = ubar_subwindow{i}{m}(ind_full{i}{m}); %values for u for all full profiles
        ind_full_notoutlier = intersect(ind_nooutlier_zq{i}{m},ind_full{i}{m}); %get full, non-outlier intervals
        ubar_full_notoutlier_analysis{i}{m} = ubar_subwindow{i}{m}(ind_full_notoutlier); %values of u for only full, non-outlier intervals for fit
    end
end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - SALTATION HEIGHT %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute median zq and sigma_zq for q(z) profile fits

%real intervals
zq_real_analysis = cell(N_Sites,1);
sigma_zq_real_analysis = cell(N_Sites,1);
mean_zq_real_analysis = cell(N_Sites,1);
median_zq_real_analysis = cell(N_Sites,1);
median_sigma_zq_real_analysis = cell(N_Sites,1);
std_zq_real_analysis = cell(N_Sites,1);

%full intervals
median_zq_full_analysis = cell(N_Sites,1);
median_sigma_zq_full_analysis = cell(N_Sites,1);
std_zq_full_analysis = cell(N_Sites,1);

%intervals with no outliers
zq_nooutlier_analysis = cell(N_T_subwindow,1);
std_zq_nooutlier_analysis = cell(N_Sites,1);

%full intervals with no outliers
zq_full_nooutlier_analysis = cell(N_T_subwindow,1);
median_zq_full_nooutlier_analysis = cell(N_Sites,1);
std_zq_full_nooutlier_analysis = cell(N_Sites,1);

for i = 1:N_Sites
    
    %real intervals
    zq_real_analysis{i} = cell(N_T_subwindow,1);
    sigma_zq_real_analysis{i} = cell(N_T_subwindow,1);
    mean_zq_real_analysis{i} = zeros(N_T_subwindow,1);
    median_zq_real_analysis{i} = zeros(N_T_subwindow,1);
    median_sigma_zq_real_analysis{i} = zeros(N_T_subwindow,1);
    std_zq_real_analysis{i} = zeros(N_T_subwindow,1);

    %full intervals
    median_zq_full_analysis{i} = zeros(N_T_subwindow,1);
    median_sigma_zq_full_analysis{i} = zeros(N_T_subwindow,1);
    std_zq_full_analysis{i} = zeros(N_T_subwindow,1);
    
    %intervals with no outliers
    zq_nooutlier_analysis{i} = cell(N_T_subwindow,1); %values for all profiles excluding outliers
    std_zq_nooutlier_analysis{i} = zeros(N_T_subwindow,1); %value for profiles excluding outliers

    %full profiles with no outliers
    zq_full_nooutlier_analysis{i} = cell(N_T_subwindow,1); %values for all full profiles excluding outliers
    median_zq_full_nooutlier_analysis{i} = zeros(N_T_subwindow,1); %median value for full profiles excluding outliers 
    std_zq_full_nooutlier_analysis{i} = zeros(N_T_subwindow,1); %value for full profiles excluding outliers 

    for m = 1:N_T_subwindow
        
        %real intervals
        zq_real = zq_subwindow{i}{m}(ind_real{i}{m}); %get zq - real 
        sigma_zq_real = sigma_zq_subwindow{i}{m}(ind_real{i}{m}); %get sigma_zq - real
        zq_real_analysis{i}{m} = zq_real; %add these to list
        sigma_zq_real_analysis{i}{m} = sigma_zq_real; %add these to list
        mean_zq_real_analysis{i}(m) = mean(zq_real); %compute median and add to list
        median_zq_real_analysis{i}(m) = median(zq_real); %compute median and add to list
        median_sigma_zq_real_analysis{i}(m) = median(sigma_zq_real); %compute median and add to list      
        std_zq_real_analysis{i}(m) = std(zq_real);

        %full intervals
        zq_full = zq_subwindow{i}{m}(ind_full{i}{m}); %get zq - full
        sigma_zq_full = sigma_zq_subwindow{i}{m}(ind_full{i}{m}); %get sigma_zq - full
        median_zq_full_analysis{i}(m) = median(zq_full); %compute median and add to list
        median_sigma_zq_full_analysis{i}(m) = median(sigma_zq_full); %compute median and add to list 
        std_zq_full_analysis{i}(m) = std(zq_full);
        
        %intervals with no outliers
        zq_nooutlier = zq_subwindow{i}{m}(ind_nooutlier_zq{i}{m});
        zq_nooutlier_analysis{i}{m} = zq_nooutlier;
        std_zq_nooutlier_analysis{i}(m) = std(zq_nooutlier);

        %full profiles with no outliers
        zq_full_nooutlier = zq_subwindow{i}{m}(intersect(ind_nooutlier_zq{i}{m},ind_full{i}{m}));
        zq_full_nooutlier_analysis{i}{m} = zq_full_nooutlier;
        median_zq_full_nooutlier_analysis{i}(m) = median(zq_full_nooutlier); %median value for full profiles excluding outliers 
        std_zq_full_nooutlier_analysis{i}(m) = std(zq_full_nooutlier);
    end
end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - TOTAL FLUX %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%


%% compute relative (fractional) differences between Qsum and Qfit estimates
Qfit_all_analysis = Qfit_subwindow; %get all total fluxes by fit
Qsum_all_analysis = Qsum_subwindow; %get all total fluxes by summation
Qfit_real_analysis = cell(N_Sites,1); %total fluxes by fitting
Qsum_real_analysis = cell(N_Sites,1); %total fluxes by summation
Qrel_real_analysis = cell(N_Sites,1); %relative total fluxes
median_Qrel_real_analysis = cell(N_Sites,1); %median of relative total fluxes
median_abs_Qrel_real_analysis = cell(N_Sites,1); %median of absolute value of relative total fluxes 
for i = 1:N_Sites
    Qfit_real_analysis{i} = cell(N_T_subwindow,1);
    Qsum_real_analysis{i} = cell(N_T_subwindow,1);
    Qrel_real_analysis{i} = cell(N_T_subwindow,1);
    median_Qrel_real_analysis{i} = zeros(N_T_subwindow,1);
    median_abs_Qrel_real_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        Qfit_real = Qfit_subwindow{i}{m}(ind_real{i}{m}); %get Qfit
        Qsum_real = Qsum_subwindow{i}{m}(ind_real{i}{m}); %get Qsum
        Qrel_real = (Qfit_real-Qsum_real)./Qsum_real; %compute Qrel
        Qfit_real_analysis{i}{m} = Qfit_real; %add Qfit to list
        Qsum_real_analysis{i}{m} = Qsum_real; %add Qfit to list
        Qrel_real_analysis{i}{m} = Qrel_real; %add these to list
        median_Qrel_real_analysis{i}(m) = median(Qrel_real); %compute median and add to list
        median_abs_Qrel_real_analysis{i}(m) = median(abs(Qrel_real)); %compute median and add to list
    end
end

%% compute relative (fractional) differences between sigma_Qsum and sigma_Qfit estimates
sigma_Qfit_real_analysis = cell(N_Sites,1); %sigmas for total fluxes by fit
sigma_Qsum_real_analysis = cell(N_Sites,1); %sigmas for total fluxes by summation
sigma_Qratio_real_analysis = cell(N_Sites,1); %sigmas for ratio of flux calculation methods
median_sigma_Qratio_real_analysis = cell(N_Sites,1); %median of ratios of sigmas for total fluxes
for i = 1:N_Sites
    sigma_Qfit_real_analysis{i} = cell(N_T_subwindow,1);
    sigma_Qsum_real_analysis{i} = cell(N_T_subwindow,1);
    sigma_Qratio_real_analysis{i} = cell(N_T_subwindow,1);
    median_sigma_Qratio_real_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        sigma_Qfit_real = sigma_Qfit_subwindow{i}{m}(ind_real{i}{m}); %get sigma_Qfit
        sigma_Qsum_real = sigma_Qsum_subwindow{i}{m}(ind_real{i}{m}); %get sigma_Qsum
        sigma_Qratio_real = sigma_Qfit_real./sigma_Qsum_real; %compute sigma_Qratio
        sigma_Qfit_real_analysis{i}{m} = sigma_Qfit_real; %add Qfit to list
        sigma_Qsum_real_analysis{i}{m} = sigma_Qsum_real; %add Qsum to list
        sigma_Qratio_real_analysis{i}{m} = sigma_Qratio_real; %add Qratio to list
        median_sigma_Qratio_real_analysis{i}(m) = median(sigma_Qratio_real(~isnan(sigma_Qratio_real))); %compute median and add to list
    end
end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATIONS - FLUX LAW %
% %%%%%%%%%%%%%%%%%%%%%%%%%

%% perform fit to Q versus u2
Chi2nu_Q_u2_analysis = cell(N_Sites,1);
uth_Q_u2_analysis = cell(N_Sites,1);
C_Q_u2_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    Chi2nu_Q_u2_analysis{i} = zeros(N_T_subwindow,1);
    uth_Q_u2_analysis{i} = zeros(N_T_subwindow,1);
    C_Q_u2_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        u2_fit = ubar_real_analysis{i}{m}.^2;
        Q_fit = Qfit_real_analysis{i}{m};
        [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(u2_fit, Q_fit);
        uth = sqrt(-a/b); C = b;
        Q_residuals = Q_pred - Q_fit; %residuals between observed and predicted q
        Chi2 = sum((Q_residuals./sigma_Q_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
        df = length(Q_fit)-2; %compute degrees of freedom for Qfit
        Chi2nu_Q_u2_analysis{i}(m) = Chi2/df;
        uth_Q_u2_analysis{i}(m) = uth;
        C_Q_u2_analysis{i}(m) = C;
    end
end

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE ANALYSIS VARIABLES %
% %%%%%%%%%%%%%%%%%%%%%%%%%

% add list of values
StartTime_all_analysis = StartTime_subwindow; %get all start times


%% save analysis data
save(SaveData_Path,'N_T_subwindow','T_subwindow_s','SiteNames','Sites','N_Sites','*analysis'); 