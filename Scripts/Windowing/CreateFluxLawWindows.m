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
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};

%% load window data and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
LoadData_Path = strcat(folder_LoadData,'DataWindows_30min'); %path for 30 minute data
load(LoadData_Path);
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% information about where to save data and plots
folder_SaveData = '../../AnalysisData/FluxLaw/'; %folder for outputs of this analysis
SaveData_Path = strcat(folder_SaveData,'FluxLawWindows_30min'); %path for 30 minute data
folder_WenglorPlots = '../../PlotOutput/WenglorFluxProfiles/'; %folder for plots

%% initialize variable lists

%initialize lists of flux values
Q_all = cell(N_Sites,1); %total flux
sigma_Q_all = cell(N_Sites,1); %total flux uncertainty
zq_all = cell(N_Sites,1); %characteristic flux height
sigma_zq_all = cell(N_Sites,1); %characteristic flux height uncertainty
fD_all = cell(N_Sites,1); %Wenglor detection frequency list - 1 second
fQ_all = cell(N_Sites,1); %Wenglor transport frequency matrix - 1 second
lambda_all = cell(N_Sites,1); %Wenglor arrival rate for flux frequency correction
qbar_all = cell(N_Sites,1); %mean partial flux
sigma_qbar_all = cell(N_Sites,1); %mean partial flux uncertainty
nbar_all = cell(N_Sites,1); %mean counts rate
sigma_nbar_all = cell(N_Sites,1); %mean uncertainty in counts rate
Cqnbar_all = cell(N_Sites,1); %mean calibration factor list
sigma_Cqnbar_all = cell(N_Sites,1); %uncertainty in calibration factor list
Qbinary_avg_all = cell(N_Sites,1); %1 second wind average flux occurrence timeseries
Chi2_Qfit_all = cell(N_Sites,1); %calculate Chi2 value for flux profiles
df_Qfit_all = cell(N_Sites,1); %calculate degrees of freedom for fitting flux profiles

%initialize lists of wind values for lowest anemometer
zL_all = cell(N_Sites,1); %stability parameter
ustRe_all = cell(N_Sites,1); %u* for Reynolds stress
sigma_ustRe_all = cell(N_Sites,1); %uncertainty in u*
tauRe_all = cell(N_Sites,1); %tau for Reynolds stress
sigma_tauRe_all = cell(N_Sites,1); %uncertainty in tau
zs_all = cell(N_Sites,1); %observed roughness height from lowest anemometer

%initialize lists of other values
ubar_all = cell(N_Sites,1); %mean wind velocity from lowest anemometer
uw_all = cell(N_Sites,1); %mean wind product from lowest anemometer

%initialize plot
close all;
h = figure;
set(h,'visible','off');
set(gca,'FontSize',16,'XMinorTick','On','YMinorTick','On','Box','On','YScale','log');
set(gca,'LooseInset',get(gca,'TightInset'));
xlabel('Wenglor height, $$z_i$$ (m)','Interpreter','Latex');
ylabel('30-minute partial saltation flux, $$\tilde{q}_i$$ (g m$$^{-2}$$ s$$^{-1}$$)','Interpreter','Latex');
set(gcf, 'PaperPosition',[0 0 8 6]);

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
    fD_all{i} = zeros(N_Windows,1)*NaN; %Wenglor detection frequency - 1 s
    fQ_all{i} = zeros(N_Windows,1)*NaN; %Wenglor total transport frequency - 1 s
    lambda_all{i} = zeros(N_Windows,1)*NaN; %Wenglor particle arrival rate - 1 s
    qbar_all{i} = cell(N_Windows,1); %mean partial flux
    sigma_qbar_all{i} = cell(N_Windows,1); %mean partial flux uncertainty
    nbar_all{i} = cell(N_Windows,1); %mean counts rate
    sigma_nbar_all{i} = cell(N_Windows,1); %uncertainty in mean counts rate
    Cqnbar_all{i} = cell(N_Windows,1); %calibration factor
    sigma_Cqnbar_all{i} = cell(N_Windows,1); %uncertainty in calibration factor
    Qbinary_avg_all{i} = cell(N_Windows,1); %1 second wind average flux occurrence timeseries
    Chi2_Qfit_all{i} = zeros(N_Windows,1); %calculate Chi2 for flux profiles
    df_Qfit_all{i} = zeros(N_Windows,1); %calculate degrees of freedom for fitting flux profiles
    
    %wind values
    zL_all{i} = zeros(N_Windows,1)*NaN; %stability parameter
    ustRe_all{i} = zeros(N_Windows,1)*NaN; %u* for Reynolds stress
    sigma_ustRe_all{i} = zeros(N_Windows,1)*NaN; %uncertainty in u*
    tauRe_all{i} = zeros(N_Windows,1)*NaN; %tau for Reynolds stress
    sigma_tauRe_all{i} = zeros(N_Windows,1)*NaN; %uncertainty in tau
    zs_all{i} = zeros(N_Windows,1)*NaN; %observed roughness height from lowest anemometer
    
    %other values
    ubar_all{i} = zeros(N_Windows,1)*NaN; %observed wind velocity from lowest anemometer
    uw_all{i} = zeros(N_Windows,1)*NaN; %mean wind product from lowest anemometer
    
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

        %get times
        t_flux = t_flux_window{i}{j};
        
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

            %plot flux profile fit
            cla; hold on;
            errorbar(zW_fit,qbar_fit,sigma_qbar_fit,'b+','MarkerSize',10);
            plot(zW_fit,q_pred,'k');
            legend('data','fit','Location','NorthEast');
            %legend('data',['fit, \chi^2_{\nu} = ',num2str((Chi2_Qfit/df_Qfit),'%.2f')],'Location','NorthEast');
            title([SiteNames{i},', ',datestr(StartTime, 'yyyy-mm-dd HH:MM'),' - ',datestr(EndTime, 'HH:MM')]);
            print([folder_WenglorPlots,'WenglorFluxProfile_',Sites{i},'_',int2str(j),'.png'],'-dpng');
        else %otherwise, set to NaN
            Q=NaN;
            sigma_Q = NaN;
            zq=NaN;
            sigma_q0=NaN;
            sigma_zq=NaN;
            Chi2_Qfit=NaN;
            df_Qfit=NaN;
        end
        
        %convert to 0 if Q<Q_min or if Q=NaN AND expected Q<Q_min g/m/s based on lowest Wenglor and zq = 10 cm
        if Q<Q_min||(isnan(Q)&&...
                ((zq_Q_min*qbar_unique(1))/exp(-zW_unique(1)/zq_Q_min)<Q_min))
            Q=0;
            sigma_Q=0;
            zq=0;
            sigma_q0=0;
            sigma_zq=0;
        end

        %Determine flux frequencies
        N = sum(n,2); %interpolated total counts rate (counts/s)
        [n_deltat, t_deltat_n] = window_average(n, t_flux, deltat_fQ); %window-averaged counts rate timeseries (counts/s)
        flux_err_binary = zeros(length(N),1); %initialize timeseries with 0's and 1's for error points
        flux_err_binary(ind_flux_err_window{i}{j}) = 1; %set error indices to 1
        flux_err_binary_deltat = window_average(flux_err_binary, t_flux, deltat_fQ); %binary flux errors
        ind_noerr_deltat = find(flux_err_binary_deltat==0); %non-error indices are window-averaged points == 0
        N_noerr_deltat = n_deltat(ind_noerr_deltat); %get counts for non-error window averaged points
        Qbinary = N_noerr_deltat~=0; %flux occurence timeseries
        fD = sum(Qbinary>0)/length(Qbinary); %detection rate
        N_bar = mean(N_noerr_deltat); %mean total counts rate
        %estimate particle arrival rate per averaging window
        if fD<=fD_min %set to zero if below detection limit
            lambda = 0;
        else %otherwise estimate arrival rate based on fD
            lambda = N_bar*seconds(deltat_fQ)/fD;
        end
        %estimate flux activity from flux detection rate and particle arrival rate
        if lambda==0 %if no (or negligible) flux, set frequencies to zero
            fQ = 0;
        elseif fD==1 %if fD = 1, set fQ to 1
            fQ = 1;
        else %otherwise, estimate fQ based on other parameters
            fQ = fD/(1-exp(-lambda)); %calculate fQ
            if fQ>1
                fQ = 1; %if correction gives fQ>1, just set as 1
            end
        end
        
        %add to lists of all fD, fQ, and lambda
        fD_all{i}(j) = fD;
        fQ_all{i}(j) = fQ;
        lambda_all{i}(j) = lambda;

        %add to list
        Q_all{i}(j) = Q; %total flux
        sigma_Q_all{i}(j) = sigma_Q; %uncertainty in total flux
        zq_all{i}(j) = zq; %characteristic flux height
        sigma_zq_all{i}(j) = sigma_zq; %uncertainty in flux height
        qbar_all{i}{j} = qbar; %partial flux
        sigma_qbar_all{i}{j} = sigma_qbar; %partial flux uncertainty
        nbar_all{i}{j} = nbar; %counts rate
        sigma_nbar_all{i}{j} = sigma_nbar; %counts rate uncertainty
        Cqnbar_all{i}{j} = Cqnbar; %calibration factor
        Qbinary_avg_all{i}{j} = Qbinary; %flux occurence timeseries
        Chi2_Qfit_all{i}(j) = Chi2_Qfit; %normalized Chi2 for flux profile fit
        df_Qfit_all{i}(j) = df_Qfit; %degrees of freedom for flux profile fit
        
        %% WIND CALCULATIONS FOR INTERVAL - BASE ANEMOMETER

        %get anemometer height
        zU = zU_window{i}(j);
        
        %get wind speed
        u = u_window{i}{j}; %uninterpolated wind
        u_int = u_int_window{i}{j}; %interpolated wind
        w = w_window{i}{j}; %uninterpolated wind
        w_int = w_int_window{i}{j}; %interpolated wind
        Temp = Temp_window{i}{j};
        
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
        
%         %compute stress uncertainty using Marcelo's code - use interpolated data for continuity
%         uw_int = (u_int-mean(u_int)).*(w_int-mean(w_int));
%         [sigma_uw_bar, ~] = random_error(uw_int, length(uw_int), 1/dt_wind_window(i), zU, ubar, 0, 0, 0, 1); %use script from Salesky et al (2012)
%         sigma_tauRe_all{i}(j) = rho_a(i)*sigma_uw_bar; %uncertainty in tau
%         sigma_ustRe_all{i}(j) = sigma_uw_bar/(2*ustRe); %uncertainty in u*
               
        %add velocity values to list
        ubar_all{i}(j) = ubar;
        uw_all{i}(j) = mean(u.*w);
        
        %add stress values to list - using raw values
        ustRe_all{i}(j) = ustRe;
        tauRe_all{i}(j) = tauRe;

        %add zs to list - using raw values
        zs_all{i}(j) = zU*exp(-kappa*ubar/ustRe); %observed roughness height from lowest anemometer

        %add stability value to list
        zL_all{i}(j) = zL;
    end
end

%%RENAME ADDITIONAL VALUES NEEDED FOR FLUX LAW ANALYSIS
theta_all = theta_window; %wind angles
W_ID_all = W_ID_window; %Wenglor IDs
Date_all = Date_window; %dates
zU_all = zU_window; %anemometer heights
zW_min_all = zW_min_window; %mininum Wenglor height
zW_max_all = zW_max_window; %maximum Wenglor height

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','*_all'); %save reduced file in GitHub folder