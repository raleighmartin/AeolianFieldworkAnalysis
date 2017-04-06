%% INITIALIZATION
%initialize
clearvars;
close all;

%% PARAMETERS AND INPUTS
kappa = 0.4; %von Karman parameter
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
g = 9.8; %gravity m/s^2

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
%%%%%%%%%%%%%%%%
% CALCULATIONS %
%%%%%%%%%%%%%%%%

%% compute relative differences of between Qsum and Qfit estimates
Qrel_analysis = cell(N_Sites,1);
median_Qrel_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    Qrel_analysis{i} = cell(N_T_subwindow,1);
    median_Qrel_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        Qrel_analysis{i}{m} = abs(Qfit_subwindow{i}{m}-Qsum_subwindow{i}{m})./Qsum_subwindow{i}{m};
        median_Qrel_analysis{i}(m) = median(Qrel_analysis{i}{m}(~isnan(Qrel_analysis{i}{m})));
    end
end

%% compute median chi2nu for fit q(z) profile fits
Chi2nu_analysis = cell(N_Sites,1);
median_Chi2nu_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    Chi2nu_analysis{i} = cell(N_T_subwindow,1);
    median_Chi2nu_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        Chi2nu_analysis{i}{m} = Chi2_Qfit_subwindow{i}{m}./df_Qfit_subwindow{i}{m};
        median_Chi2nu_analysis{i}(m) = median(Chi2nu_analysis{i}{m}(~isnan(Chi2nu_analysis{i}{m})));
    end
end

%% compute fraction NaN obtained for q(z) fits
fNaN_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    fNaN_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        fNaN_analysis{i}(m) = length(find(isnan(Chi2nu_analysis{i}{m})))/length(Chi2nu_analysis{i}{m});
    end
end

%% perform fit to Q versus u2
Chi2nu_Q_u2_analysis = cell(N_Sites,1);
uth_Q_u2_analysis = cell(N_Sites,1);
C_Q_u2_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    Chi2nu_Q_u2_analysis{i} = zeros(N_T_subwindow,1);
    uth_Q_u2_analysis{i} = zeros(N_T_subwindow,1);
    C_Q_u2_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_fit = find(Qsum_subwindow{i}{m}>1);
        u2_fit = ubar_subwindow{i}{m}(ind_fit).^2;
        Q_fit = Qsum_subwindow{i}{m}(ind_fit);
        [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(u2_fit, Q_fit);
        uth = sqrt(-a/b); C = b;
        Q_residuals = Q_pred - Q_fit; %residuals between observed and predicted q
        Chi2 = sum((Q_residuals./sigma_Q_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
        df = length(ind_fit)-2; %compute degrees of freedom for Qfit
        Chi2nu_Q_u2_analysis{i}(m) = Chi2/df;
        uth_Q_u2_analysis{i}(m) = uth;
        C_Q_u2_analysis{i}(m) = C;
    end
end

%% compute standard deviation of fitted zq values
std_zq_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    std_zq_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        std_zq_analysis{i}(m) = std(zq_subwindow{i}{m}(~isnan(zq_subwindow{i}{m})));
    end
end

%% compute median zq
median_zq_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    median_zq_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        zq_all = zq_subwindow{i}{m}(~isnan(zq_subwindow{i}{m}));
        median_zq_analysis{i}(m) = median(zq_all);
    end
end

%% compute relative uncertainty of q(z) values
sigmaq_rel_analysis = cell(N_Sites,1);
median_sigmaq_rel_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    sigmaq_rel_analysis{i} = cell(N_T_subwindow,1);
    median_sigmaq_rel_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        N_subwindow = length(sigma_qbar_subwindow{i}{m});
        sigmaq_rel_analysis{i}{m} = zeros(N_subwindow,1);
        for j = 1:N_subwindow
            sigmaq_rel_analysis{i}{m}(j) = mean(sigma_qbar_subwindow{i}{m}{j}./qbar_subwindow{i}{m}{j});
        end
        median_sigmaq_rel_analysis{i}(m) = median(sigmaq_rel_analysis{i}{m}(~isnan(sigmaq_rel_analysis{i}{m})));
    end
end

%% compute median fraction of Wenglors used in profiles (and fraction with 100% used)
fW_analysis = cell(N_Sites,1); %fraction of Wenglors used
median_fW_analysis = cell(N_Sites,1); %median fraction of Wenglors used
fraction_fW_full_analysis = cell(N_Sites,1); %fraction of time intervals with 100% usage
for i = 1:N_Sites
    fW_analysis{i} = cell(N_T_subwindow,1);
    median_fW_analysis{i} = zeros(N_T_subwindow,1);
    fraction_fW_full_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        fW = N_z_fit_subwindow{i}{m}./N_Wenglor_subwindow{i}{m};
        fW_analysis{i}{m} = fW;
        ind_calc = find(Qsum_subwindow{i}{m}>0);
        median_fW_analysis{i}(m) = median(fW(ind_calc));
        fraction_fW_full_analysis{i}(m) = length(find(fW(ind_calc)==1))/length(fW(ind_calc));
    end
end

%% compute standard deviation of fitted zq values, but using only those with all Wenglor heights represented
std_zq_fW1_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    std_zq_fW1_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_calc = find(fW_analysis{i}{m}==1);
        std_zq_fW1_analysis{i}(m) = std(zq_subwindow{i}{m}(ind_calc));
    end
end

%% compute standard deviation of fitted zq values, but using only those with all Wenglor heights represented
std_zq_fW1_positive_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    std_zq_fW1_positive_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_calc = intersect(find(fW_analysis{i}{m}==1),find(zq_subwindow{i}{m}>0));
        std_zq_fW1_positive_analysis{i}(m) = std(zq_subwindow{i}{m}(ind_calc));
    end
end

%% compute standard deviation of fitted zq values, but using only those with all Wenglor heights represented
std_zq_fW1_nooutliers_analysis = cell(N_Sites,1);
for i = 1:N_Sites
    std_zq_fW1_nooutliers_analysis{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_calc = intersect(find(fW_analysis{i}{m}==1),intersect(find(zq_subwindow{i}{m}>0),find(zq_subwindow{i}{m}<0.5)));
        std_zq_fW1_nooutliers_analysis{i}(m) = std(zq_subwindow{i}{m}(ind_calc));
    end
end

%% rename additional variables before saving
StartTime_analysis = StartTime_subwindow;
ubar_analysis = ubar_subwindow;
qbar_analysis = qbar_subwindow;
sigma_qbar_analysis = sigma_qbar_subwindow;
zW_analysis = zW_subwindow;
sigma_zW_analysis = sigma_zW_subwindow;
Qsum_analysis = Qsum_subwindow;
Qfit_analysis = Qfit_subwindow;
zq_analysis = zq_subwindow;
sigma_zq_analysis = sigma_zq_subwindow;

%% save analysis data
save(SaveData_Path,'N_T_subwindow','T_subwindow_s','SiteNames','Sites','N_Sites','*analysis'); 