%% SCRIPT TO GENERATE WINDOWS OF FLUX ACTIVITY, SALTATION FLUX, AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
rho_a = 1.18; %air density kg/m^3
g = 9.8; %gravity m/s^2
kappa = 0.4; %von Karman parameter
z0 = 1e-4; %aerodynamic roughness length (m)
Q_min = 0.05; %detection limit for Q, set to zero if below this value
zq_Q_min = 0.10; %assumed saltation height for detection limit for exponential profile for detection limit for individual Wenglor
zW_min = 0.018; %minimum Wenglor height (m) = 1.5*height of instrument (to allow one full instrument height between bottom of instrument and bed)
u_sigma_max = 5; %maximum standard deviation in total wind for error detection
zW_limit = 3; %limit on number of unique Wenglor heights in profile
dt_W = 0.04; %time interval for Wenglor (s)

%% set time interval for computing wind/flux/activity windows
dt_fQ = duration(0,0,1); %time interval for window averaging to compute fQ

%% information about where to load/save data, plots, and functions
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_WindowingData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_WenglorPlots = '../../PlotOutput/WenglorFluxProfiles/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for windowing data
TimeWindow_Path = strcat(folder_WindowingData,'TimeWindows_30min_Thresholds'); %path for loading time windows
SaveData_Path = strcat(folder_WindowingData,'ThresholdWindows_Analysis'); %path for saving output data

%% load time windows
load(TimeWindow_Path);
SiteNames = {'Jericoacoara';'Rancho Guadalupe';'Oceano'};
%SiteNames = {'Jericoacoara'};
%N_Sites = 1;

%% load processed data for each site, add to cell arrays of all sites 
WindData = cell(N_Sites,1);
FluxData = cell(N_Sites,1);
WeatherData = cell(N_Sites,1);
FluxBSNE = cell(N_Sites,1);
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path); %load processed data
    WindData{i} = ProcessedData.(AnemometerType{i}).(BaseAnemometer{i}); %data only for base anemometer
    FluxData{i} = ProcessedData.FluxWenglor; %Wenglor flux data
    WeatherData{i} = ProcessedData.Weather.WS; %all weather station data
    FluxBSNE{i} = ProcessedData.FluxBSNE; %BSNE flux data
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
end

%% initialize variable lists

%initialize lists of flux values
Q_all = cell(N_Sites,1); %total flux
sigma_Q_all = cell(N_Sites,1); %total flux uncertainty
zq_all = cell(N_Sites,1); %characteristic flux height
sigma_zq_all = cell(N_Sites,1); %characteristic flux height uncertainty
fD_all = cell(N_Sites,1); %Wenglor detection frequency list - 1 second
fQ_all = cell(N_Sites,1); %Wenglor transport frequency matrix - 1 second
q_all = cell(N_Sites,1); %partial flux
sigma_q_all = cell(N_Sites,1); %partial flux uncertainty
n_all = cell(N_Sites,1); %counts rate
sigma_n_all = cell(N_Sites,1); %uncertainty in counts rate
W_ID_all = cell(N_Sites,1); %get IDs of Wenglors in profiles
zW_all = cell(N_Sites,1); %Wenglor flux height list
sigma_zW_all = cell(N_Sites,1); %Wenglor flux height uncertainty list
N_zW_unique_all = cell(N_Sites,1); %number of unique Wenglor heights
zW_min_all = cell(N_Sites,1); %Wenglor min height list
zW_max_all = cell(N_Sites,1); %Wenglor min height list
Cqn_all = cell(N_Sites,1); %calibration factor list
Qbinary_avg_all = cell(N_Sites,1); %1 second wind average flux occurrence timeseries
Chi2_Qfit_all = cell(N_Sites,1); %calculate Chi2 value for flux profiles
df_Qfit_all = cell(N_Sites,1); %calculate degrees of freedom for fitting flux profiles

%initialize lists of wind values for lowest anemometer
zU_all = cell(N_Sites,1); %height of anemometer
theta_all = cell(N_Sites,1); %mean theta
zL_all = cell(N_Sites,1); %stability parameter
ustRe_all = cell(N_Sites,1); %u* for Reynolds stress
sigma_ustRe_all = cell(N_Sites,1); %uncertainty in u*
tauRe_all = cell(N_Sites,1); %tau for Reynolds stress
sigma_tauRe_all = cell(N_Sites,1); %uncertainty in tau
uth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of u_th
tauth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of tau_th
zs_all = cell(N_Sites,1); %observed roughness height from lowest anemometer

%initialize lists of other values
RH_all = cell(N_Sites,1); %relative humidity for weather station
ubar_all = cell(N_Sites,1); %mean wind velocity from lowest anemometer
uw_all = cell(N_Sites,1); %mean wind product from lowest anemometer


%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    
    %% get start times and end times from file
    BlockStartTimes = StartTime_all{i};
    BlockEndTimes = EndTime_all{i};
    N_Blocks = length(BlockStartTimes);
    
    
    %% initialize lists of values
    
    %flux values
    Q_all{i} = zeros(N_Blocks,1)*NaN; %total flux
    sigma_Q_all{i} = zeros(N_Blocks,1)*NaN; %total flux uncertainty
    zq_all{i} = zeros(N_Blocks,1)*NaN; %characteristic flux height
    sigma_zq_all{i} = zeros(N_Blocks,1)*NaN; %characteristic flux height uncertainty
    fD_all{i} = zeros(N_Blocks,1)*NaN; %Wenglor detection frequency - 1 s
    fQ_all{i} = zeros(N_Blocks,1)*NaN; %Wenglor total transport frequency - 1 s
    q_all{i} = cell(N_Blocks,1); %partial flux
    sigma_q_all{i} = cell(N_Blocks,1); %partial flux uncertainty
    n_all{i} = cell(N_Blocks,1); %counts rate
    sigma_n_all{i} = cell(N_Blocks,1); %uncertainty in counts rate
    W_ID_all{i} = cell(N_Blocks,1); %get names of Wenglors in profiles
    zW_all{i} = cell(N_Blocks,1); %Wenglor height
    sigma_zW_all{i} = cell(N_Blocks,1); %Wenglor height uncertainty
    N_zW_unique_all{i} = zeros(N_Blocks,1); %number of unique Wenglor heights
    zW_min_all{i} = zeros(N_Blocks,1); %minimum Wenglor height
    zW_max_all{i} = zeros(N_Blocks,1); %maximum Wenglor height
    Cqn_all{i} = cell(N_Blocks,1); %calibration factor
    Qbinary_avg_all{i} = cell(N_Blocks,1); %1 second wind average flux occurrence timeseries
    Chi2_Qfit_all{i} = zeros(N_Blocks,1); %calculate Chi2 for flux profiles
    df_Qfit_all{i} = zeros(N_Blocks,1); %calculate degrees of freedom for fitting flux profiles
    
    %wind values
    zU_all{i} = zeros(N_Blocks,1)*NaN; %height of anemometer
    theta_all{i} = zeros(N_Blocks,1)*NaN; %mean theta
    zL_all{i} = zeros(N_Blocks,1)*NaN; %stability parameter
    ustRe_all{i} = zeros(N_Blocks,1)*NaN; %u* for Reynolds stress
    sigma_ustRe_all{i} = zeros(N_Blocks,1)*NaN; %uncertainty in u*
    tauRe_all{i} = zeros(N_Blocks,1)*NaN; %tau for Reynolds stress
    sigma_tauRe_all{i} = zeros(N_Blocks,1)*NaN; %uncertainty in tau
    zs_all{i} = zeros(N_Blocks,1)*NaN; %observed roughness height from lowest anemometer
    
    %other values
    RH_all{i} = zeros(N_Blocks,1)*NaN; %relative humidity
    ubar_all{i} = zeros(N_Blocks,1)*NaN; %observed wind velocity from lowest anemometer
    uw_all{i} = zeros(N_Blocks,1)*NaN; %mean wind product from lowest anemometer
    
    %threshold values
    uth_TFEM_all{i} = zeros(N_Blocks,1)*NaN; %1 second TFEM calc of threshold
    tauth_TFEM_all{i} = zeros(N_Blocks,1)*NaN; %1 second TFEM calc of threshold stress
    
    %% go through time blocks
    for j = 1:N_Blocks

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Blocks),', ',datestr(now)]

        %get specific start and end time
        StartTime = BlockStartTimes(j);
        EndTime = BlockEndTimes(j);

        %get duration of interval in seconds
        T_interval = seconds(EndTime-StartTime);
        
        %% FLUX CALCULATIONS FOR INTERVAL

        %extract time interval - get interval number and indices within interval for analysis
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(FluxData{i},StartTime,EndTime,'t','t','t');

        %use only longest interval
        ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd))); %index of longest interval
        IntervalN = IntervalN(ind_longest); %interval number for longest interval
        IntervalInd = IntervalInd{ind_longest}; %indices for longest interval

        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(FluxData{i}(IntervalN).t.t,FluxData{i}(IntervalN).t.err); %indices of error times
        IntervalInd = setdiff(IntervalInd,ErrInd); %indices for longest interval excluding error times
        t_Interval = FluxData{i}(IntervalN).t.t(IntervalInd); %times for longest interval excluding error times

        %get Wenglor heights
        zW_profile = mean(FluxData{i}(IntervalN).qz.z(IntervalInd,:)); %compute mean Wenglor heights
        ind_zW_positive = find(zW_profile>zW_min); %find indices of points definitely above zero (based on zW_min)
        zW_profile = mean(FluxData{i}(IntervalN).qz.z(IntervalInd,ind_zW_positive)); %compute mean Wenglor heights (but now only positive values)
        sigma_zW_profile = mean(FluxData{i}(IntervalN).qz.sigma_z(IntervalInd,ind_zW_positive)); %compute uncertainty in Wenglor heights (but now only positive values)
        
        %get raw qz, Cqn, sigma_Cqn and n values
        qz = FluxData{i}(IntervalN).qz.qz(IntervalInd,ind_zW_positive); %partial fluxes
        Cqn = FluxData{i}(IntervalN).qz.Cqn(IntervalInd,ind_zW_positive); %calibration factors
        sigma_Cqn = FluxData{i}(IntervalN).qz.sigma_Cqn(IntervalInd,ind_zW_positive); %uncertainty in calibration factors
        n = FluxData{i}(IntervalN).qz.n(IntervalInd,ind_zW_positive); %particle counts per increment (not per second)
        N = sum(sum(n)); %total detected particles for all Wenglors

        %get profile information for fitting
        q_profile = mean(qz); %compute mean qz profile
        n_profile = mean(n)/dt_W; %compute mean particle counts per second for each qz
        sigma_n_profile = sqrt(n_profile/T_interval); %compute uncertainty in mean particle counts
        Cqn_profile = mean(Cqn); %compute mean calibration coefficient for each qz
        sigma_Cqn_profile = mean(sigma_Cqn); %compute mean calibration coefficient uncertainty for each qz
        sigma_q_profile = sqrt((sigma_Cqn_profile.*n_profile).^2+(sigma_n_profile.*Cqn_profile).^2); %compute uncertainty in qz profile, from both counts and calibration coefficient
        
        %get associated Wenglor IDs only for heights included in profile
        W_ID = FluxData{i}(IntervalN).qz.WenglorID(ind_zW_positive);
        
        %deal with repeated values
        zW_unique = unique(zW_profile); %unique heights
        N_zW_unique = length(zW_unique); %number of unique heights
        q_unique = zeros(size(zW_unique)); %initialize partial fluxes corresponding to unique heights 
        sigma_q_unique = zeros(size(zW_unique)); %initialize uncertainties in partial fluxes corresponding to unique heights
        sigma_zW_unique = zeros(size(zW_unique)); %initialize uncertainties in unique heights
        for k = 1:N_zW_unique; %go through each unique height
            ind_zW = find(zW_profile==zW_unique(k)); %find all Wenglors associated with unique height
            if length(ind_zW)>1 %if more than one Wenglor for unique height, compute mean and uncertainty for repeated heights
                q_min = (Q_min/zW_unique(k))*exp(-zW_unique(k)/zq_Q_min); %estimate q associated with min detection limit for this height
                q_bar = mean(q_profile(ind_zW)); %get mean q for this height
                if q_bar < q_min; %if mean q at height is below detection limit
                    q = 0; sigma_q = 0; %then just set values to zero
                else %otherwise, use script "MeanUncertainty.m"
                    [q, sigma_q] = MeanUncertainty(q_profile(ind_zW), sigma_q_profile(ind_zW)); %compute flux mean and uncertainy for repeated heights
                end
                sigma_zW = mean(sigma_zW_profile(ind_zW)); %no reduction in uncertainty for height, because all uncertainties are correlated, so just take mean of values
                q_unique(k) = q; %set mean q as partial flux for unique height
                sigma_q_unique(k) = sigma_q; %set sigma_q as uncertainty in partial flux for unique height
                sigma_zW_unique(k) = sigma_zW; %set sigma_zW as uncertainty in unique height
            else %otherwise, if only Wenglor at given height, just use existing values
                q_unique(k) = q_profile(ind_zW); %unique q is single partial flux value
                sigma_q_unique(k) = sigma_q_profile(ind_zW); %unique q uncertainty is single sigma_q value
                sigma_zW_unique(k) = sigma_zW_profile(ind_zW); %unique zW uncertainty is single zW value
            end
        end
        
        %Fit to q profile, remove 0 values for fitting
        ind_fit = find(q_unique>0); %get indices of partial flux greater than 0
        N_fit = length(ind_fit); %get number of available values for fit
        q_fit = q_unique(ind_fit); %get nonzero partial fluxes for fit
        zW_fit = zW_unique(ind_fit); %get corresponding heights for fit
        sigma_q_fit = sigma_q_unique(ind_fit); %get corresponding nonzero partial flux uncertainties for fit
        sigma_zW_fit = zeros(1,N_fit); %when fitting, neglect uncertainty in Wenglor height, which is already accounted for by calibration
        
        %Perform profile fit to get q0, zq, and Q if sufficient points for fitting
        if N_fit>=zW_limit %if number of nonzero points exceeds zW_limit, perform fit 
            [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,q_pred,sigma_q_pred] = qz_profilefit(q_fit,zW_fit,sigma_q_fit,sigma_zW_fit); %apply qz_profilefit routine
            q_residuals = q_pred - q_fit; %residuals between observed and predicted q
            Chi2_Qfit = sum((q_residuals./sigma_q_pred).^2); %compute Chi2 based on residuals (Bevington and Robinson, Eq. 8.4)
            df_Qfit = N_fit-2; %compute degrees of freedom for Qfit
        else %if insufficient values for fitting, set to NaN
            Q=NaN;
            sigma_Q = NaN;
            zq=NaN;
            sigma_q0=NaN;
            sigma_zq=NaN;
            Chi2_Qfit=NaN;
            df_Qfit=NaN;
        end
        
        %convert to 0 if Q<Q_min or if Q=NaN AND expected Q<Q_min g/m/s based on lowest Wenglor and zq = 10 cm
        if Q<Q_min||(isnan(Q)&&... %condition of small or NaN Q 
                ((zq_Q_min*q_unique(1))/exp(-zW_unique(1)/zq_Q_min)<Q_min)); %condition of lowest partial flux very small
            Q=0; %set Q to zero
            sigma_Q=0; %set sigma_Q to zero
            q0=0; %set q0 to zero
            zq=0; %set zq to zero
            sigma_q0=0; %set q0 uncertainty to zero
            sigma_zq=0; %set zq uncertainty to zero
        end

        %add to list
        Q_all{i}(j) = Q; %total flux
        sigma_Q_all{i}(j) = sigma_Q; %uncertainty in total flux
        zq_all{i}(j) = zq; %characteristic flux height
        sigma_zq_all{i}(j) = sigma_zq; %uncertainty in flux height
        q_all{i}{j} = q_profile; %partial flux
        sigma_q_all{i}{j} = sigma_q_profile; %partial flux uncertainty
        n_all{i}{j} = n_profile; %counts rate
        sigma_n_all{i}{j} = sigma_n_profile; %counts rate uncertainty
        W_ID_all{i}{j} = W_ID; %get names of Wenglors in profiles
        zW_all{i}{j} = zW_profile; %Wenglor flux height
        sigma_zW_all{i}{j} = sigma_zW_profile; %uncertainty in Wenglor flux height
        N_zW_unique_all{i}(j) = N_zW_unique; %number of unique Wenglor heights
        zW_min_all{i}(j) = min(zW_profile); %minimum Wenglor height
        zW_max_all{i}(j) = max(zW_profile); %maximum Wenglor height
        Cqn_all{i}{j} = Cqn_profile; %calibration factor
                
        %FLUX FREQUENCY CALCULATIONS FOR INTERVAL
        
        %whole window (30-minute) flux frequencies
        nsum = sum(n')'; %create vector nsum, which is sum of all n's for each time increment
        nsum_avg = window_average(nsum, t_Interval, dt_fQ); %compute 1s window average
        Qbinary_avg = nsum_avg~=0; %flux occurence timeseries
        if Q==0 %if no flux, set frequencies to zero
            fD_all{i}(j) = 0;
            fQ_all{i}(j) = 0;
        else %otherwise, do calculations
            T = length(nsum_avg); %number of timesteps

            %get frequency of detection
            fD = (sum(nsum_avg>0))/T;

            %estimate fQ iteratively 
            fQ_prev = 0; %previous estimate of fQ (arbitrarily set to 0 for while loop functioning) 
            fQ = fD; %start by assuming fQ=fD
            while abs(fQ-fQ_prev)>0.001 %compare current and previous estimations of fQ for while loop
                fQ_prev = fQ; %update fQ_prev with last fQ
                fQ = fD/(1-exp(-sum(n_profile)/(fQ*dt_W))); %calculate new fQ
            end

            %remove points with fQ>1
            if fQ>1
                fQ = NaN;
            end

            %add to lists of all fD and fQ
            fD_all{i}(j) = fD;
            fQ_all{i}(j) = fQ;
        end

        %partial window (25-second) flux frequencies
        
        %add to list
        Qbinary_avg_all{i}{j} = Qbinary_avg; %flux occurence timeseries
        Chi2_Qfit_all{i}(j) = Chi2_Qfit; %normalized Chi2 for flux profile fit
        df_Qfit_all{i}(j) = df_Qfit; %degrees of freedom for flux profile fit
        
        %% WIND CALCULATIONS FOR INTERVAL - BASE ANEMOMETER
        %extract time interval
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData{i},StartTime,EndTime,'u','int','int');
        
        %use only first interval
        IntervalN = IntervalN(1);
        IntervalInd = IntervalInd{1};

        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(WindData{i}(IntervalN).t.int,WindData{i}(IntervalN).t.err);
        IntervalInd_NoErr = setdiff(IntervalInd,ErrInd);

        %get anemometer height
        zU = WindData{i}(IntervalN).z.z;

        %get velocity values using no error points
        u = WindData{i}(IntervalN).u.int(IntervalInd_NoErr);
        v = WindData{i}(IntervalN).v.int(IntervalInd_NoErr);
        w = WindData{i}(IntervalN).w.int(IntervalInd_NoErr);
        Temp = WindData{i}(IntervalN).T.int(IntervalInd_NoErr);

        %remove additional error points based on large deviations in wind velocity
        u_total = sqrt(u.^2+v.^2+w.^2); %total wind
        u_total_max = mean(u_total)+u_sigma_max*std(u_total); %maximum total wind based on multiple of std dev
        ind_good_wind = find(u_total<=u_total_max); %get indices of points with total wind below upper limit
        u = u(ind_good_wind); %keep only good u
        v = v(ind_good_wind); %keep only good v
        w = w(ind_good_wind); %keep only good w
        Temp = Temp(ind_good_wind); %keep only good temperatures
                
        %compute wind angle
        theta_all{i}(j) = atan(mean(v./u))*180/pi;

        %rotate instrument, call these 'rot' values
        [u_rot, ~, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

        %make computations
        u_bar = mean(u_rot);
        w_bar = mean(w_rot);
        T_bar = mean(Temp);
        uw = (u_rot-u_bar).*(w_rot-w_bar); %u'w' product
        tauRe_kernal = mean(uw);
        if tauRe_kernal<=0
            ustRe = sqrt(-tauRe_kernal);
            tauRe = -rho_a*tauRe_kernal;
        else
            tauRe = NaN;
            ustRe = NaN;
        end
        heatflux = mean((Temp-T_bar).*(w_rot-w_bar));
        zL = (-(g./(T_bar+273.15)).*heatflux)./(ustRe.^3./(kappa*zU));

        %compute stress uncertainty using Marcelo's code
        [sigma_uw_bar, ~] = random_error(uw, length(uw), 1/dt_u_s(i), zU, u_bar, 0, 0, 0, 1); %use script from Salesky et al (2012)
        sigma_tauRe_all{i}(j) = rho_a*sigma_uw_bar; %uncertainty in tau
        sigma_ustRe_all{i}(j) = sigma_uw_bar/(2*ustRe); %uncertainty in u*
        % Q: Would it be better to use interpolated data for random_error.m, since it examines temporal features of the dataset?
        % see Vickers and Mahrt (1997) for best procedure for data despiking
               
        %add velocity values to list
        zU_all{i}(j) = zU; %height of anemometer
        ubar_all{i}(j) = u_bar;
        uw_all{i}(j) = mean(u_rot.*w_rot);
        
        %add stress values to list - using raw values
        ustRe_all{i}(j) = ustRe;
        tauRe_all{i}(j) = tauRe;

        %add zs to list - using raw values
        zs_all{i}(j) = zU*exp(-kappa*u_bar/ustRe); %observed roughness height from lowest anemometer

        %add stability value to list
        zL_all{i}(j) = zL;

        %% create averaged / low-pass wind timeseries, get velocity values including error points
        u_witherror = WindData{i}(IntervalN).u.int(IntervalInd);
        v_witherror = WindData{i}(IntervalN).v.int(IntervalInd);
        w_witherror = WindData{i}(IntervalN).w.int(IntervalInd);
        t_witherror = WindData{i}(IntervalN).t.int(IntervalInd);
        [u_rot_witherror, ~, ~] = reorient_anemometers_vanboxel2004(u_witherror, v_witherror, w_witherror); %rotate instrument
        u_avg = window_average(u_rot_witherror, t_witherror, dt_fQ); %compute 1 second window average wind timeseries - using raw values

        
        %% HUMIDITY VALUE FOR INTERVAL
        %extract humidity values
        [H_Interval, ~, ~, ~] = ExtractVariableTimeInterval(WeatherData{i},StartTime,EndTime,'H','int','int');
        RH_all{i}(j) = mean(H_Interval);
        
        
        %% TFEM threshold calcs
        if fQ_all{i}(j)~=1 %calculation only valid if 0<fQ<1
            uth_TFEM = prctile(u_avg,100*(1-fQ_all{i}(j))); %generate uth from TFEM method based on fQ and 1 s window averaged u
            tauth_TFEM = rho_a*kappa^2*uth_TFEM^2/log(zU/z0)^2; %generate tauth from TFEM method based on uth_TFEM 
        end

        %add thresholds to lists
        uth_TFEM_all{i}(j) = uth_TFEM; %uth from TFEM
        tauth_TFEM_all{i}(j) = tauth_TFEM; %tauth from TFEM
    end
end

% SAVE DATA
save(SaveData_Path,'Sites','SiteNames','*all'); %save reduced file in GitHub folder