%% SCRIPT TO GENERATE WINDOWS OF SALTATION FLUX AND STRESS VALUES FOR ANALYSIS
% FUNCTION DEPENDENCIES: ExtractVariableTimeInterval.m, 
% qz_profilefit, reorient_anemometers_vanboxel2004, window_average

%% initialize
clearvars;

%% parameter values
rho_a = 1.18; %air density kg/m^3
g = 9.8; %gravity m/s^2
kappa = 0.4; %von Karman parameter
z0 = [1e-4 1e-4 1e-4]; %aerodynamic roughness length (m)

%% information about sites for analysis
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
AnemometerType = {'Ultrasonic';'Ultrasonic';'Sonic'};
BaseAnemometer = {'U1';'U1';'S1'};
dt_u_s = [0.04; 0.04; 0.02]; %time interval of sonic (s)
N_Sites = length(Sites);


%% set time interval for computing wind/flux windows
dt_fQ = duration(0,0,1); %time interval for window averaging to compute fQ


%% information about where to load/save data, plots, and functions
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
TimeWindow_Path = strcat(folder_AnalysisData,'TimeWindows'); %path for loading time windows
SaveData_Path = strcat(folder_AnalysisData,'StressFluxWindows'); %path for saving output data
BSNEData_Path = strcat(folder_AnalysisData,'FluxBSNE'); %path for saving output data
addpath(folder_Functions); %point MATLAB to location of functions


%% load time windows
load(TimeWindow_Path);


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
zW_all = cell(N_Sites,1); %Wenglor flux height list
sigma_zW_all = cell(N_Sites,1); %Wenglor flux height uncertainty list
qcal_all = cell(N_Sites,1); %calibration factor list
Qbinary_avg_all = cell(N_Sites,1); %1 second wind average flux occurrence timeseries

%initialize lists of wind values for lowest anemometer
zU_all = cell(N_Sites,1); %height of anemometer
theta_all = cell(N_Sites,1); %mean theta
zL_all = cell(N_Sites,1); %stability parameter
ustRe_all = cell(N_Sites,1); %u* for Reynolds stress
tauRe_all = cell(N_Sites,1); %tau for Reynolds stress
uth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of u_th
tauth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of tau_th
zs_all = cell(N_Sites,1); %observed roughness height from lowest anemometer

%initialize lists of other values
RH_all = cell(N_Sites,1);
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
    zW_all{i} = cell(N_Blocks,1); %Wenglor height
    sigma_zW_all{i} = cell(N_Blocks,1); %Wenglor height uncertainty
    qcal_all{i} = cell(N_Blocks,1); %calibration factor
    Qbinary_avg_all{i} = cell(N_Blocks,1); %1 second wind average flux occurrence timeseries
    
    %wind values
    zU_all{i} = zeros(N_Blocks,1)*NaN; %height of anemometer
    theta_all{i} = zeros(N_Blocks,1)*NaN; %mean theta
    zL_all{i} = zeros(N_Blocks,1)*NaN; %stability parameter
    ustRe_all{i} = zeros(N_Blocks,1)*NaN; %calibrated u* for Reynolds stress
    tauRe_all{i} = zeros(N_Blocks,1)*NaN; %calibrated tau for Reynolds stress
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


    %% FLUX CALCULATIONS FOR INTERVAL

    %extract time interval - get interval number and indices within interval for analysis
    [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(FluxData{i},StartTime,EndTime,'t','t','t');

        %use only longest interval
        ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
        IntervalN = IntervalN(ind_longest);
        IntervalInd = IntervalInd{ind_longest};
        StartTime_Extraction = FluxData{i}(IntervalN).t.t(IntervalInd(1));
        EndTime_Extraction = FluxData{i}(IntervalN).t.t(IntervalInd(end));

        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(FluxData{i}(IntervalN).t.t,FluxData{i}(IntervalN).t.err);
        IntervalInd = setdiff(IntervalInd,ErrInd);
        t_Interval = FluxData{i}(IntervalN).t.t(IntervalInd);

        %get raw qz and qcal values, determine n from this
        qz = FluxData{i}(IntervalN).qz.qz(IntervalInd,:);
        qcal = FluxData{i}(IntervalN).qz.qzPerCount(IntervalInd,:);
        sigma_qcal = FluxData{i}(IntervalN).qz.sigma_qzPerCount(IntervalInd,:);
        n = FluxData{i}(IntervalN).qz.n(IntervalInd,:); %particle arrival rates
        N = sum(sum(n)); %total detected particles

        %get profile information for fitting
        q_profile = mean(qz); %compute mean qz profile
        sigma_q_n_profile = abs(q_profile.*(1./sqrt(sum(n)))); %contribution of uncertainty in counting particles
        sigma_q_qcal_profile = abs(q_profile.*(mean(sigma_qcal)./mean(qcal))); %contribution of uncertainty in calibration coefficient
        sigma_q_profile = sqrt(sigma_q_n_profile.^2+sigma_q_qcal_profile.^2); %compute uncertainty in qz profile
        %sigma_q_profile = mean(sigma_qcal.*n); %compute uncertainty in qz profile
        zW_profile = mean(FluxData{i}(IntervalN).qz.z(IntervalInd,:)); %compute mean Wenglor heights
        sigma_zW_profile = mean(FluxData{i}(IntervalN).qz.sigma_z(IntervalInd,:)); %compute uncertainty in Wenglor heights

        %deal with repeated values
        zW_unique = unique(zW_profile);
        q_unique = zeros(size(zW_unique));
        sigma_q_unique = zeros(size(zW_unique));
        sigma_zW_unique = zeros(size(zW_unique));
        for k = 1:length(zW_unique);
            ind_zW = find(zW_profile==zW_unique(k));
            q_unique(k) = mean(q_profile(ind_zW));
            sigma_q_unique(k) = mean(sigma_q_profile(ind_zW))/sqrt(length(ind_zW)); %reduce error based on number of values
            sigma_zW_unique(k) = mean(sigma_zW_profile(ind_zW))/sqrt(length(ind_zW)); %reduce error based on number of values
        end

        %Perform profile fit to get q0, zq, and Q
        ind_fit = intersect(find(q_unique>0),find(zW_unique>0)); %only use values with q>0 and zW>0
        [q0,zq,sigma_q0,sigma_zq] = qz_profilefit(q_unique(ind_fit),zW_unique(ind_fit),sigma_q_unique(ind_fit),sigma_zW_unique(ind_fit));
        Q = q0*zq; %get total flux [g/m/s]
        sigma_Q = sqrt((sigma_q0*zq)^2+(sigma_zq*q0)^2); %estimate uncertainty in total flux

        %convert to 0 if Q<0.05 or if Q=NaN AND expected Q<0.05 g/m/s based on lowest Wenglor and zq = 5 cm
        if Q<0.05||(isnan(Q)&&((0.1*q_unique(1))/exp(-zW_unique(1)/0.1)<0.05));
            Q=0;
            zq=0;
            sigma_q0=0;
            sigma_zq=0;
        end

        %Determine flux frequencies for different window averaging times
        qsum = sum(qz')'; %create vector qsum, which is sum of all q's for each time increment
        if Q==0 %if no flux, set frequencies to zero
            fD_all{i}(j) = 0;
            fQ_all{i}(j) = 0;
        else %otherwise, do calculations
            qsum_avg = window_average(qsum, t_Interval, dt_fQ); %compute window average
            T = length(qsum_avg); %number of timesteps

            %get frequency of detection
            fD = (sum(qsum_avg>0))/T;

            %estimate fQ iteratively 
            fQ_prev = 0; %previous estimate of fQ (arbitrarily set to 0 for while loop functioning) 
            fQ = fD; %start by assuming fQ=fD
            while abs(fQ-fQ_prev)>0.001 %compare current and previous estimations of fQ for while loop
                fQ_prev = fQ; %update fQ_prev with last fQ
                fQ = fD/(1-exp(-N/(fQ_prev*T))); %calculate new fQ
            end

            %remove points with fQ>1
            if fQ>1
                fQ = NaN;
            end

            %add to lists of all fD and fQ
            fD_all{i}(j) = fD;
            fQ_all{i}(j) = fQ;
        end

        %Timesteps with flux for 1 second wind average
        qsum_avg = window_average(qsum, t_Interval, dt_fQ);
        Qbinary_avg = qsum_avg~=0; %flux occurence timeseries

        %get calibration values
        qcal_profile = mean(qcal);

        %add to list
        Q_all{i}(j) = Q; %total flux
        sigma_Q_all{i}(j) = sigma_Q; %uncertainty in total flux
        zq_all{i}(j) = zq; %characteristic flux height
        sigma_zq_all{i}(j) = sigma_zq; %uncertainty in flux height
        q_all{i}{j} = q_profile; %partial flux
        sigma_q_all{i}{j} = sigma_q_profile; %partial flux uncertainty
        zW_all{i}{j} = zW_profile; %Wenglor flux height
        sigma_zW_all{i}{j} = sigma_zW_profile; %uncertainty in Wenglor flux height
        qcal_all{i}{j} = qcal_profile; %calibration factor
        Qbinary_avg_all{i}{j} = Qbinary_avg; %flux occurence timeseries

        
        %% WIND CALCULATIONS FOR INTERVAL - BASE ANEMOMETER
        %extract time interval
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData{i},StartTime,EndTime,'u','int','int');
        
        %use only first interval
        IntervalN = IntervalN(1);
        IntervalInd = IntervalInd{1};
        StartTime_Extraction = WindData{i}(IntervalN).t.int(IntervalInd(1));
        EndTime_Extraction = WindData{i}(IntervalN).t.int(IntervalInd(end));

        %further reduce IntervalInd based on eliminating error times
        [~, ErrInd, ~] = intersect(WindData{i}(IntervalN).t.int,WindData{i}(IntervalN).t.err);
        IntervalInd_NoErr = setdiff(IntervalInd,ErrInd);

        %get anemometer height
        zU = WindData{i}(IntervalN).z.z;

        %get velocity values using no error points
        u = WindData{i}(IntervalN).u.int(IntervalInd_NoErr);
        v = WindData{i}(IntervalN).v.int(IntervalInd_NoErr);
        w = WindData{i}(IntervalN).w.int(IntervalInd_NoErr);
        t = WindData{i}(IntervalN).t.int(IntervalInd_NoErr);
        T_raw = WindData{i}(IntervalN).T.int(IntervalInd_NoErr);

        %compute wind angle
        theta_all{i}(j) = atan(mean(v./u))*180/pi;

        %rotate instrument, call these 'raw' values
        [u_raw, ~, w_raw] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

        %make computations - raw values
        u_bar_raw = mean(u_raw);
        w_bar_raw = mean(w_raw);
        T_bar_raw = mean(T_raw);
        ustRe_kernal = mean((u_raw-u_bar_raw).*(w_raw-w_bar_raw));
        if ustRe_kernal<=0
            ustRe_raw = sqrt(-ustRe_kernal);
        else
            ustRe_raw = NaN;
        end
        heatflux_raw = mean((T_raw-T_bar_raw).*(w_raw-w_bar_raw));
        zL = (-(g./(T_bar_raw+273.15)).*heatflux_raw)./(ustRe_raw.^3./(kappa*zU));

        %add velocity values to list
        zU_all{i}(j) = zU; %height of anemometer
        ubar_all{i}(j) = u_bar_raw;
        uw_all{i}(j) = mean(u_raw.*w_raw);
        
        %add stress values to list - using raw values
        ustRe_all{i}(j) = ustRe_raw;
        tauRe_all{i}(j) = rho_a*ustRe_raw.^2;

        %add zs to list - using raw values
        zs_all{i}(j) = zU*exp(-kappa*u_bar_raw/ustRe_raw); %observed roughness height from lowest anemometer

        %add stability value to list
        zL_all{i}(j) = zL;

        %% create averaged / low-pass wind timeseries, get velocity values including error points
        u_witherror = WindData{i}(IntervalN).u.int(IntervalInd);
        v_witherror = WindData{i}(IntervalN).v.int(IntervalInd);
        w_witherror = WindData{i}(IntervalN).w.int(IntervalInd);
        t_witherror = WindData{i}(IntervalN).t.int(IntervalInd);
        [u_raw_witherror, ~, ~] = reorient_anemometers_vanboxel2004(u_witherror, v_witherror, w_witherror); %rotate instrument
        u_avg = window_average(u_raw_witherror, t, dt_fQ); %compute 1 second window average wind timeseries - using raw values

        
        %% HUMIDITY VALUE FOR INTERVAL
        %extract humidity values
        [H_Interval, ~, ~, ~] = ExtractVariableTimeInterval(WeatherData{i},StartTime,EndTime,'H','int','int');
        RH_all{i}(j) = mean(H_Interval);
        
        
        %% TFEM threshold calcs
        if fQ_all{i}(j)~=1 %calculation only valid if 0<fQ<1
            uth_TFEM = prctile(u_avg,100*(1-fQ_all{i}(j))); %generate uth from TFEM method based on fQ and 1 s window averaged u
            tauth_TFEM = rho_a*kappa^2*uth_TFEM^2/log(zU/z0(i))^2; %generate tauth from TFEM method based on uth_TFEM 
        end

        %add thresholds to lists
        uth_TFEM_all{i}(j) = uth_TFEM; %uth from TFEM
        tauth_TFEM_all{i}(j) = tauth_TFEM; %tauth from TFEM
    end
end

% SAVE DATA
save(SaveData_Path,'Sites','*all'); %save reduced file in GitHub folder
save(BSNEData_Path,'FluxBSNE'); %save BSNE data to avoid having to open full fill in future