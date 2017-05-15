%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
%Q_min = 0.05; %detection limit for Q, set to zero if below this value
%zq_Q_min = 0.10; %assumed saltation height for detection limit for exponential profile for detection limit for individual Wenglor
zW_min = 3; %limit on number of unique Wenglor heights in profile

%% folders for loading and saving data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

%% paths for loading and saving data - restricted
LoadData_Path = strcat(folder_LoadData,'DataFullSubwindows_30min_Restricted'); %path for 30 minute data - for thresholds analysis
SaveData_Path = strcat(folder_SaveData,'DataFullSubwindowCalcs_30min_Restricted'); %path for 30 minute data - for thresholds analysis

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize flux arrays
N_Wenglor_subwindow = cell(N_Sites,1); %number of Wenglors in subwindow
N_z_subwindow = cell(N_Sites,1); %number of unique Wenglor heights for subwindow
N_z_fit_subwindow = cell(N_Sites,1); %number of Wenglor heights for fitting in subwindow
Qfit_subwindow = cell(N_Sites,1); %total flux by fitting for subwindow
sigma_Qfit_subwindow = cell(N_Sites,1); %total flux by fitting uncertainty for subwindow
zq_subwindow = cell(N_Sites,1); %flux height by fitting for subwindow
sigma_zq_subwindow = cell(N_Sites,1); %uncertainty in flux height by fitting for subwindow
Chi2_Qfit_subwindow = cell(N_Sites,1); %goodness of Q fit
df_Qfit_subwindow = cell(N_Sites,1); %degrees of freedom for fitting
Qsum_subwindow = cell(N_Sites,1); %total flux by summing for subwindow
sigma_Qsum_subwindow = cell(N_Sites,1); %total flux by summing uncertainty for subwindow

%% GO THROUGH SITES
for i = 1:N_Sites
   
    %% initialize variable arrays
    N_Wenglor_subwindow{i} = cell(N_T_subwindow,1); %number of Wenglors in subwindow
    N_z_subwindow{i} = cell(N_T_subwindow,1); %number of unique Wenglor heights for subwindow
    N_z_fit_subwindow{i} = cell(N_T_subwindow,1); %number of Wenglor heights for fitting in subwindow
    Qfit_subwindow{i} = cell(N_T_subwindow,1); %total flux by fitting for subwindow
    sigma_Qfit_subwindow{i} = cell(N_T_subwindow,1); %total flux by fitting uncertainty for subwindow
    Chi2_Qfit_subwindow{i} = cell(N_T_subwindow,1); %goodness of Q fit
    df_Qfit_subwindow{i} = cell(N_T_subwindow,1); %degrees of freedom for fitting
    zq_subwindow{i} = cell(N_T_subwindow,1); %flux height by fitting for subwindow
    sigma_zq_subwindow{i} = cell(N_T_subwindow,1); %uncertainty in flux height by fitting for subwindow
    Qsum_subwindow{i} = cell(N_T_subwindow,1); %total flux by summing for subwindow
    sigma_Qsum_subwindow{i} = cell(N_T_subwindow,1); %total flux by summing uncertainty for subwindow
    
    %% GO THROUGH MEASUREMENT INTERVALS
    for m = 1:N_T_subwindow

        %% get information about subwindows
        N_subwindows = length(StartTime_subwindow{i}{m}); %get number of subwindows
        
        %% display processing status
        processing_status = [Sites{i},', ',int2str(m),' of ',int2str(N_T_subwindow),', ',datestr(now)]

        %% initialize flux values
        N_Wenglor_subwindow{i}{m} = zeros(N_subwindows,1); %number of Wenglors in subwindow
        N_z_subwindow{i}{m} = zeros(N_subwindows,1); %number of unique Wenglor heights for subwindow
        N_z_fit_subwindow{i}{m} = zeros(N_subwindows,1); %number of Wenglor heights for fitting in subwindow
        Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by fitting for subwindow
        sigma_Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by fitting uncertainty for subwindow
        Chi2_Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %goodness of Q fit
        df_Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %degrees of freedom for fitting
        zq_subwindow{i}{m} = zeros(N_subwindows,1); %flux height by fitting for subwindow
        sigma_zq_subwindow{i}{m} = zeros(N_subwindows,1); %uncertainty in flux height by fitting for subwindow
        Qsum_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by summing for subwindow
        sigma_Qsum_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by summing uncertainty for subwindow
          
        %% GO THROUGH SUBWINDOWS
        for k = 1:N_subwindows
                                
            %% get flux data in subwindow
            zq_BSNE = zq_BSNE_subwindow{i}{m}(k); %zq for BSNE
            sigma_zq_BSNE = sigma_zq_BSNE_subwindow{i}{m}(k); %sigma_zq for BSNE
            zW = zW_subwindow{i}{m}{k}; %flux heights for subwindow
            sigma_zW = sigma_zW_subwindow{i}{m}{k}; %uncertainty on flux heights for subwindow
            qW = qbar_subwindow{i}{m}{k}; %partial fluxes
            sigma_qW = sigma_qbar_subwindow{i}{m}{k}; %partial flux uncertainty for subwindow
            nbar = nbar_subwindow{i}{m}{k}; %counts rates
            sigma_nbar = sigma_nbar_subwindow{i}{m}{k}; %particle counts uncertainty for subwindow
            Cqnbar = Cqnbar_subwindow{i}{m}{k}; %calibration coefficient
            sigma_Cqnbar = sigma_Cqnbar_subwindow{i}{m}{k}; %calibration coefficient

            %% perform exponential fit to get Q
            %[q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,Chi2_Qfit,df_Qfit,qW_unique,zW_unique,sigma_qW_unique,sigma_zW_unique,N_z_fit] = qz_fit_Wenglor(qW, zW, sigma_qW, sigma_zW, Q_min, zq_Q_min, zW_min, zq_BSNE);
            [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,Chi2_Qfit,df_Qfit,qW_unique,zW_unique,sigma_qW_unique,sigma_zW_unique,N_z_fit] = qz_fit_Wenglor(qW, zW, sigma_qW, sigma_zW, zW_min, zq_BSNE);
            
            %% Perform summation to get Q
            [Qsum,sigma_Qsum] = qz_summation_Wenglor(qW_unique, zW_unique, sigma_qW_unique, zq_BSNE, sigma_zq_BSNE);
            
            %% add values to arrays
            N_Wenglor_subwindow{i}{m}(k) = length(zW); %number of Wenglors for calculation
            N_z_subwindow{i}{m}(k) = length(zW_unique); %number of unique Wenglor heights for calculation
            N_z_fit_subwindow{i}{m}(k) = N_z_fit; %number of Wenglor heights for fit
            Qfit_subwindow{i}{m}(k) = Q; %total flux by fitting for subwindow
            sigma_Qfit_subwindow{i}{m}(k) = sigma_Q; %total flux by fitting uncertainty for subwindow
            Chi2_Qfit_subwindow{i}{m}(k) = Chi2_Qfit; %total flux by fitting uncertainty for subwindow
            df_Qfit_subwindow{i}{m}(k) = df_Qfit; %degrees of freedom for fitting
            zq_subwindow{i}{m}(k) = zq; %flux height by fitting for subwindow
            sigma_zq_subwindow{i}{m}(k) = sigma_zq; %uncertainty in flux height by fitting for subwindow
            Qsum_subwindow{i}{m}(k) = Qsum; %total flux by summing for subwindow
            sigma_Qsum_subwindow{i}{m}(k) = sigma_Qsum; %total flux by summing uncertainty for subwindow
        end           
    end
end

%% save analysis data
save(SaveData_Path,'SiteNames','Sites','N_Sites','*subwindow'); 