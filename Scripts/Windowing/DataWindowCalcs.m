%% SCRIPT TO CALCULATE SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
kappa = 0.4; %von Karman parameter
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
g = 9.8; %gravity m/s^2
deltat_fQ = duration(0,0,1); %sampling interval for flux activity analysis
fD_min = 0.005; %minimum detection rate for zero flux
Q_min = 0.05; %detection limit for Q, set to zero if below this value
zq_Q_min = 0.10; %assumed saltation height for detection limit for exponential profile for detection limit for individual Wenglor
zW_limit = 3; %limit on number of unique Wenglor heights in profile

%% dates with known zero flux at each site (to set Q=fQ=0 if Wenglors are inactive)
Dates_zeroflux = {[datetime(2014,11,18)],... %Jeri
    [datetime(0,0,0)],... %none for Rancho Guadalupe
    [datetime(2015,5,17),datetime(2015,5,20),datetime(2015,5,21)]}; %Oceano

%% folders for loading and saving data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions

% %% paths for loading and saving data - restricted
% LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Restricted'); %path for 30 minute data - for flux law analysis
% SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data - for flux law analysis

% %% paths for loading and saving data - restricted - with alternate base anemometer
% LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Restricted_alt'); %path for 30 minute data - for flux law analysis
% SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_30min_Restricted_alt'); %path for 30 minute data - for flux law analysis

%% paths for loading and saving data - unrestricted
LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis
SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_30min_Unrestricted'); %path for 30 minute data - for thresholds analysis

%% paths for loading and saving data - Yue
% rho_a = [1.22]; %air density kg/m^3 (assumes ~15 C at Oceano)
% LoadData_Path = strcat(folder_LoadData,'DataWindows_Oceano_Yue_1'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_Oceano_Yue_1'); %path for 30 minute data
% rho_a = [1.22]; %air density kg/m^3 (assumes ~15 C at Oceano)
% LoadData_Path = strcat(folder_LoadData,'DataWindows_Oceano_Yue_2'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_Oceano_Yue_2'); %path for 30 minute data
% rho_a = [1.22]; %air density kg/m^3 (assumes ~15 C at Oceano)
% LoadData_Path = strcat(folder_LoadData,'DataWindows_Oceano_Yue_3'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_Oceano_Yue_3'); %path for 30 minute data
% rho_a = [1.22]; %air density kg/m^3 (assumes ~15 C at Oceano)
% LoadData_Path = strcat(folder_LoadData,'DataWindows_Oceano_Yue_4'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'DataWindowCalcs_Oceano_Yue_4'); %path for 30 minute data

%% load data and functions
load(LoadData_Path); %load data
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize variable lists
%initialize lists of flux values
Q_all = cell(N_Sites,1); %total flux
sigma_Q_all = cell(N_Sites,1); %total flux uncertainty
zq_all = cell(N_Sites,1); %characteristic flux height
sigma_zq_all = cell(N_Sites,1); %characteristic flux height uncertainty
qbar_all = cell(N_Sites,1); %mean partial flux
sigma_qbar_all = cell(N_Sites,1); %mean partial flux uncertainty
nbar_all = cell(N_Sites,1); %mean counts rate
sigma_nbar_all = cell(N_Sites,1); %mean uncertainty in counts rate
Cqnbar_all = cell(N_Sites,1); %mean calibration factor list
sigma_Cqnbar_all = cell(N_Sites,1); %uncertainty in calibration factor list
Chi2_Qfit_all = cell(N_Sites,1); %calculate Chi2 value for flux profiles
df_Qfit_all = cell(N_Sites,1); %calculate degrees of freedom for fitting flux profiles

%initialize lists of wind values for lowest anemometer
zL_all = cell(N_Sites,1); %stability parameter
ustRe_all = cell(N_Sites,1); %u* for Reynolds stress
sigma_ustRe_all = cell(N_Sites,1); %uncertainty in u*
tauRe_all = cell(N_Sites,1); %tau for Reynolds stress
sigma_tauRe_all = cell(N_Sites,1); %uncertainty in tau
zsRe_all = cell(N_Sites,1); %observed roughness height from Reynolds stress
ubar_all = cell(N_Sites,1); %mean wind velocity from lowest anemometer
uwbar_all = cell(N_Sites,1); %mean wind product from lowest anemometer
  
%initiate lists of activity values
fD_all = cell(N_Sites,1); %Wenglor detection frequency list
fQ_all = cell(N_Sites,1); %Wenglor transport frequency matrix
lambda_all = cell(N_Sites,1); %Total particle arrival rate per sampling interval for flux frequency correction

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% get number of windows
    WindowStartTimes = StartTime_window{i};
    WindowEndTimes = EndTime_window{i};
    N_Windows = length(WindowStartTimes);
    
    %% initialize lists of values
    %flux values
    Q_all{i} = zeros(N_Windows,1)*NaN; %total flux
    sigma_Q_all{i} = zeros(N_Windows,1)*NaN; %total flux uncertainty
    zq_all{i} = zeros(N_Windows,1)*NaN; %characteristic flux height
    sigma_zq_all{i} = zeros(N_Windows,1)*NaN; %characteristic flux height uncertainty
    qbar_all{i} = cell(N_Windows,1); %mean partial flux
    sigma_qbar_all{i} = cell(N_Windows,1); %mean partial flux uncertainty
    nbar_all{i} = cell(N_Windows,1); %mean counts rate
    sigma_nbar_all{i} = cell(N_Windows,1); %uncertainty in mean counts rate
    Cqnbar_all{i} = cell(N_Windows,1); %calibration factor
    sigma_Cqnbar_all{i} = cell(N_Windows,1); %uncertainty in calibration factor
    Chi2_Qfit_all{i} = zeros(N_Windows,1)*NaN; %calculate Chi2 for flux profiles
    df_Qfit_all{i} = zeros(N_Windows,1)*NaN; %calculate degrees of freedom for fitting flux profiles

    %initiate lists of activity values
    fD_all{i} = zeros(N_Windows,1)*NaN; %Wenglor detection frequency - 1 s
    fQ_all{i} = zeros(N_Windows,1)*NaN; %Wenglor total transport frequency - 1 s
    lambda_all{i} = zeros(N_Windows,1)*NaN; %Wenglor particle arrival rate - 1 s
    
    %wind values - lowest anemometer
    zL_all{i} = zeros(N_Windows,1)*NaN; %stability parameter
    ustRe_all{i} = zeros(N_Windows,1)*NaN; %u* for Reynolds stress
    sigma_ustRe_all{i} = zeros(N_Windows,1)*NaN; %uncertainty in u*
    tauRe_all{i} = zeros(N_Windows,1)*NaN; %tau for Reynolds stress
    sigma_tauRe_all{i} = zeros(N_Windows,1)*NaN; %uncertainty in tau
    zsRe_all{i} = zeros(N_Windows,1)*NaN; %observed roughness height from lowest anemometer
    ubar_all{i} = zeros(N_Windows,1)*NaN; %observed wind velocity from lowest anemometer
    uwbar_all{i} = zeros(N_Windows,1)*NaN; %mean wind product from lowest anemometer
        
    %% go through time blocks
    for j = 1:N_Windows

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]

        %get specific start and end time
        StartTime = WindowStartTimes(j);
        EndTime = WindowEndTimes(j);

        %get duration of interval in seconds
        T_interval = seconds(EndTime-StartTime);
        
        %% FLUX CALCULATIONS FOR INTERVAL

        %only perform calculation if there is flux data in window
        if ~isnan(zW_window{i}{j})
            
            %% total flux calcs
           
            %get Wenglor heights
            zW = zW_window{i}{j};
            sigma_zW = sigma_zW_window{i}{j};

            %get raw q, Cqn, sigma_Cqn and n values
            q = q_window{i}{j};
            Cqn = Cqn_window{i}{j};
            sigma_Cqn = sigma_Cqn_window{i}{j};
            n = n_window{i}{j};

            %get mean q profile information for fitting
            qbar = mean(q); %compute mean qz profile
            nbar = mean(n)/dt_flux_window(i); %compute mean particle counts per second for each qz
            sigma_nbar = sqrt(nbar/T_interval); %compute uncertainty in mean particle counts
            Cqnbar = mean(Cqn); %compute mean calibration coefficient for each qz
            sigma_Cqnbar = mean(sigma_Cqn); %compute mean calibration coefficient uncertainty for each qz
            sigma_qbar = sqrt((sigma_Cqnbar.*nbar).^2+(sigma_nbar.*Cqnbar).^2); %compute uncertainty in qz profile, from counts and calibration coefficient

            %get associated Wenglor IDs only for heights included in profile
            W_ID = W_ID_window{i};

            %deal with repeated values
            zW_unique = unique(zW);
            N_zW_unique = length(zW_unique);
            qbar_unique = zeros(size(zW_unique));
            sigma_qbar_unique = zeros(size(zW_unique));
            sigma_zW_unique = zeros(size(zW_unique));
            for k = 1:N_zW_unique
                ind_zW = find(zW==zW_unique(k));
                if length(ind_zW)>1 %compute mean and uncertainty for repeated heights
                    q_min = (Q_min/zW_unique(k))*exp(-zW_unique(k)/zq_Q_min); %get min q detection limit for this height
                    qbar_z = mean(qbar(ind_zW)); %get mean q for this height
                    if qbar_z < q_min %if mean q at height is below detection limit
                        qbar_z = 0; sigma_qbar_z = 0; %then just set values to zero
                    else %otherwise, use script "MeanUncertainty.m"
                        [qbar_z, sigma_qbar_z] = MeanUncertainty(qbar(ind_zW), sigma_qbar(ind_zW)); %compute flux mean and uncertainy for repeated heights
                    end
                    sigma_zW_unique(k) = mean(sigma_zW(ind_zW)); %no reduction in uncertainty for height, because all uncertainties are correlated, so just take mean of values
                    qbar_unique(k) = qbar_z;
                    sigma_qbar_unique(k) = sigma_qbar_z;
                else %otherwise, if only Wenglor at given height, just use existing values
                    qbar_unique(k) = qbar(ind_zW);
                    sigma_qbar_unique(k) = sigma_qbar(ind_zW);
                    sigma_zW_unique(k) = sigma_zW(ind_zW);
                end
            end

            %Remove 0 values for fitting
            ind_fit = find(qbar_unique>0);
            N_fit = length(ind_fit);
            qbar_fit = qbar_unique(ind_fit);
            zW_fit = zW_unique(ind_fit);
            sigma_qbar_fit = sigma_qbar_unique(ind_fit);
            sigma_zW_fit = zeros(1,N_fit); %neglect uncertainty in Wenglor height, which is already accounted for by calibration

            %Perform profile fit to get q0, zq, and Q if sufficient points for fitting
            if N_fit>=zW_limit
                [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,q_pred,sigma_q_pred] = qz_profilefit(qbar_fit,zW_fit,sigma_qbar_fit,sigma_zW_fit);
                q_residuals = q_pred - qbar_fit; %residuals between observed and predicted q
                Chi2_Qfit = sum((q_residuals./sigma_q_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
                df_Qfit = N_fit-2; %compute degrees of freedom for Qfit
            else %otherwise, set to NaN
                Q=NaN;
                sigma_Q = NaN;
                zq=NaN;
                sigma_q0=NaN;
                sigma_zq=NaN;
                Chi2_Qfit=NaN;
                df_Qfit=NaN;
            end

            %Get window-averaged counts values for determining flux frequency
            ntotal_int = ntotal_int_window{i}{j}; %interpolated total counts rate (counts/s)
            t_flux = t_flux_int_window{i}{j}; %times for interpolated flux
            ind_flux_err = ind_flux_err_window{i}{j}; %list of error time indices for flux
            flux_err_binary = zeros(length(ntotal_int),1); %initialize timeseries with 0's and 1's for error points
            flux_err_binary(ind_flux_err) = 1; %set error indices to 1
            [ntotal_int_deltat, t_deltat_n] = window_average(ntotal_int, t_flux, deltat_fQ); %window-averaged counts rate timeseries (counts/s)
            flux_err_binary_deltat = window_average(flux_err_binary, t_flux, deltat_fQ); %binary flux errors
            ind_noerr_deltat = find(flux_err_binary_deltat==0); %get indices with no error in n
            ntotal_deltat = ntotal_int_deltat(ind_noerr_deltat); %get total window-averaged counts rate (counts/s) only for non-error points
 
            %calculate flux activity
            [fD,fQ,lambda] = CalculateFluxActivity(ntotal_deltat,seconds(deltat_fQ),fD_min);
            
            %convert to 0 if Q<Q_min or if Q=NaN AND expected Q<Q_min g/m/s based on lowest Wenglor and zq = 10 cm
            if Q<Q_min||(isnan(Q)&&...
                    ((zq_Q_min*qbar_unique(1))/exp(-zW_unique(1)/zq_Q_min)<Q_min))
                Q=0;
                sigma_Q=0;
                zq=0;
                sigma_q0=0;
                sigma_zq=0;
                fQ=0;
                sigma_fQ=0;
            end

            %add to lists of values
            Q_all{i}(j) = Q; %total flux
            sigma_Q_all{i}(j) = sigma_Q; %uncertainty in total flux
            zq_all{i}(j) = zq; %characteristic flux height
            sigma_zq_all{i}(j) = sigma_zq; %uncertainty in flux height
            qbar_all{i}{j} = qbar; %partial flux
            sigma_qbar_all{i}{j} = sigma_qbar; %partial flux uncertainty
            nbar_all{i}{j} = nbar; %counts rate
            sigma_nbar_all{i}{j} = sigma_nbar; %counts rate uncertainty
            Cqnbar_all{i}{j} = Cqnbar; %calibration factor
            Chi2_Qfit_all{i}(j) = Chi2_Qfit; %normalized Chi2 for flux profile fit
            df_Qfit_all{i}(j) = df_Qfit; %degrees of freedom for flux profile fit
            fD_all{i}(j) = fD; %flux detection rate
            fQ_all{i}(j) = fQ; %flux frequency
            lambda_all{i}(j) = lambda; %particle arrivals rate
        
        %if it's a designated day when we know there is no flux, set values accordingly
        elseif ~isempty(find(Dates_zeroflux{i}==datetime(year(StartTime),month(StartTime),day(StartTime))))
            fQ_all{i}(j) = 0;
            Q_all{i}(j) = 0;
            sigma_Q_all{i}(j) = 0;
        end
            
        %% WIND CALCULATIONS FOR INTERVAL - BASE ANEMOMETER

        %get anemometer height
        zU = zU_base_window{i}(j);
        
        %get wind speed
        u = u_window{i}{j}; %error-free wind
        w = w_window{i}{j}; %error-free wind
        Temp = Temp_window{i}{j}; %error-free temperatures
        
        %make computations
        ubar = mean(u);
        wbar = mean(w);
        Tempbar = mean(Temp);
        uw = (u-ubar).*(w-wbar); %u'w' product
        tauRe_kernal = mean(uw);
        if tauRe_kernal<=0
            ustRe = sqrt(-tauRe_kernal);
            tauRe = -rho_a(i)*tauRe_kernal;
        else
            tauRe = NaN;
            ustRe = NaN;
        end
        heatflux = mean((Temp-Tempbar).*(w-wbar));
        zL = (-(g./(Tempbar+273.15)).*heatflux)./(ustRe.^3./(kappa*zU));
        
        %compute stress uncertainty using Marcelo's code - use only error-free data
        [sigma_uw_bar, ~] = random_error(uw, length(uw), 1/dt_wind_window(i), zU, ubar, 0, 0, 0, 1); %use script from Salesky et al (2012)
        sigma_tauRe_all{i}(j) = rho_a(i)*sigma_uw_bar; %uncertainty in tau
        sigma_ustRe_all{i}(j) = sigma_uw_bar/(2*ustRe); %uncertainty in u*    
               
        %add velocity values to list
        ubar_all{i}(j) = ubar;
        uwbar_all{i}(j) = mean(u.*w);
        
        %add stress values to list - using raw values
        ustRe_all{i}(j) = ustRe;
        tauRe_all{i}(j) = tauRe;

        %add zs to list - using raw values
        zsRe_all{i}(j) = zU*exp(-kappa*ubar/ustRe); %observed roughness height from lowest anemometer

        %add stability value to list
        zL_all{i}(j) = zL;
    end
end

%calculate adjusted theta (relative to dominant wind during saltation at site)
theta_adjusted_all = cell(N_Sites,1);
theta_all = theta_window; %rename variable for wind angles
for i = 1:N_Sites
    theta_adjusted_all{i}=theta_all{i}-mean(theta_all{i}(Q_all{i}>0));
end

%%RENAME ADDITIONAL VALUES NEEDED FOR FLUX LAW ANALYSIS
W_ID_all = W_ID_window; %Wenglor IDs
Date_all = Date_window; %dates
StartTimes_all = StartTime_window; %start times
EndTimes_all = EndTime_window; %start times
zU_all = zU_base_window; %anemometer heights
zW_all = zW_window; %Wenglor height
sigma_zW_all = sigma_zW_window; %uncertainty in Wenglor height
zW_min_all = zW_min_window; %mininum Wenglor height
zW_max_all = zW_max_window; %maximum Wenglor height

% SAVE DATA
save(SaveData_Path,'Site*','N_Sites','*_all'); %save reduced file in GitHub folder