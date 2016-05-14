%% SCRIPT TO GENERATE WINDOWS OF SALTATION FLUX AND STRESS VALUES FOR ANALYSIS
% SCRIPT DEPENDENCIES: CreateStressFluxWindows.m, CreateTimeBlocks.m, 
% ExtractVariableTimeInterval.m, IntersectingTimeIntervals.m,
% qz_profilefit, reorient_anemometers_vanboxel2004, window_average


%% initialize
clearvars;

%% parameter values
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
g = 9.8; %gravity m/s^2
kappa = 0.39; %von Karman parameter
tau_it = [0.1317 0.0951 0.0841]; %impact threshold stress (Pa)
tau_ft = [0.1787 0.1631 0.1299]; %fluid threshold stress (Pa)
% tauit = [0.1 0.1 0.0864]; %impact threshold stress (Pa)
% tauft = [0.15 0.15 0.1308]; %fluid threshold stress (Pa)
ust_it = sqrt(tau_it./rho_a); %impact threshold shear velocity (m/s)
ust_ft = sqrt(tau_ft./rho_a); %fluid threshold shear velocity (m/s)
k_zs = 0.004; %proportionality constant for zs (m^2 s^2 kg^-1)
z0 = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m) at threshold

%% information about sites for analysis
%Sites = {'Jericoacoara'};
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% set time interval for computing wind/flux windows
WindowTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
N_RunningPerProfile = floor(WindowTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times


%% set dt of flux timeseries
dt_q_s = 0.04; %time interval of flux (s)


%% set durations of window average for flux frequencies
dt_fQ_min_s = 1;
dt_fQ_max_s = 1;
N_dt_fQ = 1;
% dt_fQ_min_s = 0.04; %mininum window average time
% dt_fQ_max_s = 600; %maximum window average time
% N_dt_fQ = 30; %number of durations of window average
dt_fQ_s = unique([round(logspace(0,log10(dt_fQ_max_s/dt_fQ_min_s),N_dt_fQ))*dt_fQ_min_s, 1]); %create window average dt's, make sure to include 1 second in list
N_dt_fQ = length(dt_fQ_s); %recalculate number of window average dt's after removing repeats
dt_fQ_windowaverage = duration(0,0,dt_fQ_s);
ind_dt_fQ_1s = find(dt_fQ_s==1); %index in list of dt's for 1 second


%% set number of points to use for rating curve of zs versus u
N_zs_ratingcurve = 1000;


%% information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../AnalysisData/'; %folder for storing outputs of this analysis
SaveFullData_Path = strcat(folder_ProcessedData,'StressFluxWindows_all');
SaveData_Path = strcat(folder_AnalysisData,'StressFluxWindows_all');


%% load processed and metadata for each site, add to structured arrays of all data and metadata
Data = cell(N_Sites,1);
Metadata = cell(N_Sites,1);
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Sites{i});
    Data{i} = load(ProcessedData_Path); %load processed data
    Metadata{i} = load(Metadata_Path); %load metadata
end


%% initialize variable lists

%initialize lists of grain sizes, dates, start times, and end times
d50_all = cell(N_Sites,1); %lists of surface grain sizes corresponding to calculated values
d10_all = cell(N_Sites,1);
d90_all = cell(N_Sites,1);
date_all = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_all = cell(N_Sites,1); %lists of block start times
EndTime_all = cell(N_Sites,1); %lists of block end times

%initialize lists of flux values
Q_all = cell(N_Sites,1); %total flux
sigma_Q_all = cell(N_Sites,1); %total flux uncertainty
zq_all = cell(N_Sites,1); %characteristic flux height
sigma_zq_all = cell(N_Sites,1); %characteristic flux height uncertainty
fD_all = cell(N_Sites,1); %Wenglor detection frequency matrix
fD1_all = cell(N_Sites,1); %Wenglor detection frequency list - 1 second
fQ_all = cell(N_Sites,1); %Wenglor transport frequency matrix
fQ1_all = cell(N_Sites,1); %Wenglor transport frequency matrix - 1 second
fQalt_all = cell(N_Sites,1); %Wenglor transport frequency matrix - alternate (non-iterative) calc
q_all = cell(N_Sites,1); %partial flux
sigma_q_all = cell(N_Sites,1); %partial flux uncertainty
zW_all = cell(N_Sites,1); %Wenglor flux height list
sigma_zW_all = cell(N_Sites,1); %Wenglor flux height uncertainty list
qcal_all = cell(N_Sites,1); %calibration factor list

%initialize lists of wind values for lowest anemometer
zU_all = cell(N_Sites,1); %height of anemometer
u_bar_all = cell(N_Sites,1); %mean u values
u_std_all = cell(N_Sites,1); %standard deviation u values
u2_all = cell(N_Sites,1); %mean u^2
theta_all = cell(N_Sites,1); %mean theta
zL_all = cell(N_Sites,1); %stability parameter
ustRe_all = cell(N_Sites,1); %u* for Reynolds stress
tauRe_all = cell(N_Sites,1); %tau for Reynolds stress
f_uaboveft_all = cell(N_Sites,1); %frequency of u above fluid threshold
f_ubelowit_all = cell(N_Sites,1); %frequency of u below impact threshold
f_ubetween_all = cell(N_Sites,1); %frequency of u between thresholds
f_ubetween_frombelow_all = cell(N_Sites,1); %frequency of u between thresholds approaching from below it
f_ubetween_fromabove_all = cell(N_Sites,1); %frequency of u between thresholds approaching from above ft
uth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of u_th
tauth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of tau_th
zs_all = cell(N_Sites,1); %observed roughness height from lowest anemometer

%stress partition values
eta_inactive_1s_Q0_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, times of zero flux - based on 1s winds, using fQ
eta_inactive_fQ_1s_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, constant tau_it - based on 1s winds, using fQ
eta_inactive_fQ_1s_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on 1s winds, using fQ
eta_inactive_fQ_LP_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, constant tau_it - based on LP winds, using fQ
eta_inactive_fQ_LP_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on LP winds, using fQ
eta_active_fQ_1s_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using constant tau_it - based on 1s winds, using fQ
eta_active_fQ_1s_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using TFEM u_th - based on 1s winds, using fQ
eta_active_fQ_LP_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using constant tau_it - based on LP winds, using fQ
eta_active_fQ_LP_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using TFEM u_th - based on LP winds, using fQ
eta_inactive_fU_1s_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, constant tau_it - based on 1s winds, using fU
eta_inactive_fU_1s_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on 1s winds, using fU
eta_inactive_fU_LP_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, constant tau_it - based on LP winds, using fU
eta_inactive_fU_LP_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on LP winds, using fU
eta_active_fU_1s_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using constant tau_it - based on 1s winds, using fU
eta_active_fU_1s_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using TFEM u_th - based on 1s winds, using fU
eta_active_fU_LP_constthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using constant tau_it - based on LP winds, using fU
eta_active_fU_LP_TFEMthr_all = cell(N_Sites,1); %mean stress partition based on mean stress during active transport, using TFEM u_th - based on LP winds, using fU

%inialize arrays of low-pass timeseries
q_LP_all = cell(N_Sites,1); %low-pass flux timeseries
u_LP_all = cell(N_Sites,1); %low-pass wind timeseries
t_q_LP_all = cell(N_Sites,1); %times corresponding to low-pass flux
t_u_LP_all = cell(N_Sites,1); %times corresponding to low-pass wind

%initialize arrays of window average timeseries
Qbinary_1s_all = cell(N_Sites,1); %1 second wind average flux occurrence timeseries
u_1s_all = cell(N_Sites,1); %1 second window average wind velocity timeseries
t_1s_all = cell(N_Sites,1); %times for 1 second window average timeseries

% %initialize lists of profile values
% u_profile_all = cell(N_Sites,1); %calibrated velocity profile
% z_profile_all = cell(N_Sites,1); %heights for profile
% Anemometer_profile_all = cell(N_Sites,1); %names of anemometers for profile
% u_lowest_all = cell(N_Sites,1); %calibrated velocity profile - lowest anemometers
% z_lowest_all = cell(N_Sites,1); %heights of lowest anemometers
% Anemometer_lowest_all = cell(N_Sites,1); %names of lowest anemometers for profile
% ustLog_all = cell(N_Sites,1); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
% z0_profile_all = cell(N_Sites,1); %z0 calculated from mean of calibrated velocities for lowest anemometers in velocity profile

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% choose anemometer type based on site of interest
    if strcmp(Sites{i},'Oceano')
        AnemometerType = 'Sonic';
        BaseAnemometer = 'S1';
        dt_u_s = 0.02; %time interval of sonic (s)
    elseif strcmp(Sites{i},'RanchoGuadalupe')
        AnemometerType = 'Ultrasonic';
        BaseAnemometer = 'U1';
        dt_u_s = 0.04; %time interval of sonic (s)
    elseif strcmp(Sites{i},'Jericoacoara')
        AnemometerType = 'Ultrasonic';
        BaseAnemometer = 'U1';
        dt_u_s = 0.04; %time interval of sonic (s)
    end
    
    
    %% extract wind and flux data from overall processed data file
    WindDataAll = Data{i}.ProcessedData.(AnemometerType); %all wind data
    Anemometer_profile = fieldnames(WindDataAll); %get names of anemometers in profile
    N_Anemometers = length(Anemometer_profile); %get number of anemometers in profile
    WindDataBase = WindDataAll.(BaseAnemometer); %data only for base anemometer
    FluxData = Data{i}.ProcessedData.FluxWenglor; %Wenglor flux data
    
    
    %% get start times and end times for flux and wind observations, using times for base anemometer
    WindStartTimes = [WindDataBase.StartTime]';
    WindEndTimes = [WindDataBase.EndTime]';
    FluxStartTimes = [FluxData.StartTime]';
    FluxEndTimes = [FluxData.EndTime]';
    
    %get start and end times for intersecting flux and wind intervals
    [StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(WindStartTimes,WindEndTimes,FluxStartTimes,FluxEndTimes);
    
    %create time blocks based on intersecting time intervals
    [BlockStartTimes, BlockEndTimes] = ...
            CreateTimeBlocks(StartTimesIntersecting, EndTimesIntersecting, WindowTimeInterval);
    
    %Add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [BlockStartTimes_j, BlockEndTimes_j] = ...
            CreateTimeBlocks(StartTimesIntersecting+(j*RunningTimeInterval), EndTimesIntersecting, WindowTimeInterval);
        BlockStartTimes = [BlockStartTimes; BlockStartTimes_j];
        BlockEndTimes = [BlockEndTimes; BlockEndTimes_j];
    end
    BlockStartTimes = sort(BlockStartTimes);
    BlockEndTimes = sort(BlockEndTimes);
    N_Blocks = length(BlockStartTimes);
    
    
    %% initialize lists of values

    %lists of dates, times, and grain size values
    d50_all{i} = zeros(N_Blocks,1); %lists of surface grain sizes corresponding to calculated values
    d10_all{i} = zeros(N_Blocks,1);
    d90_all{i} = zeros(N_Blocks,1);
    date_all{i} = datetime(BlockStartTimes.Year, BlockStartTimes.Month, BlockStartTimes.Day); %lists of dates corresponding to calculated values
    StartTime_all{i} = BlockStartTimes;
    EndTime_all{i} = BlockEndTimes;
    
    %flux values
    Q_all{i} = zeros(N_Blocks,1)*NaN; %total flux
    sigma_Q_all{i} = zeros(N_Blocks,1)*NaN; %total flux uncertainty
    zq_all{i} = zeros(N_Blocks,1)*NaN; %characteristic flux height
    sigma_zq_all{i} = zeros(N_Blocks,1)*NaN; %characteristic flux height uncertainty
    fD_all{i} = zeros(N_Blocks,N_dt_fQ)*NaN; %Wenglor detection frequency
    fD1_all{i} = zeros(N_Blocks,1)*NaN; %Wenglor detection frequency - 1 s
    fQ_all{i} = zeros(N_Blocks,N_dt_fQ)*NaN; %Wenglor total transport frequency
    fQ1_all{i} = zeros(N_Blocks,1)*NaN; %Wenglor total transport frequency - 1 s
    fQalt_all{i} = zeros(N_Blocks,N_dt_fQ)*NaN; %Wenglor total transport frequency, alternate calc
    q_all{i} = cell(N_Blocks,1); %partial flux
    sigma_q_all{i} = cell(N_Blocks,1); %partial flux uncertainty
    zW_all{i} = cell(N_Blocks,1); %Wenglor height
    sigma_zW_all{i} = cell(N_Blocks,1); %Wenglor height uncertainty
    qcal_all{i} = cell(N_Blocks,1); %calibration factor
    
    %initialize lists of values for lowest anemometer
    zU_all{i} = zeros(N_Blocks,1)*NaN; %height of anemometer
    u_bar_all{i} = zeros(N_Blocks,1)*NaN; %raw mean u values
    u_std_all{i} = zeros(N_Blocks,1)*NaN; %standard deviation u values
    u2_all{i} = zeros(N_Blocks,1)*NaN; %mean u^2
    theta_all{i} = zeros(N_Blocks,1)*NaN; %mean theta
    zL_all{i} = zeros(N_Blocks,1)*NaN; %stability parameter
    ustRe_all{i} = zeros(N_Blocks,1)*NaN; %calibrated u* for Reynolds stress
    tauRe_all{i} = zeros(N_Blocks,1)*NaN; %calibrated tau for Reynolds stress
    f_uaboveft_all{i} = zeros(N_Blocks,1)*NaN; %frequency of u above fluid threshold
    f_ubelowit_all{i} = zeros(N_Blocks,1)*NaN; %frequency of u below impact threshold
    f_ubetween_all{i} = zeros(N_Blocks,1)*NaN; %frequency of u between thresholds
    f_ubetween_frombelow_all{i} = zeros(N_Blocks,1)*NaN; %frequency of u between thresholds approaching from below it
    f_ubetween_fromabove_all{i} = zeros(N_Blocks,1)*NaN; %frequency of u between thresholds approaching from above ft
    uth_TFEM_all{i} = zeros(N_Blocks,1)*NaN; %1 second TFEM calc of threshold
    tauth_TFEM_all{i} = zeros(N_Blocks,1)*NaN; %1 second TFEM calc of threshold stress
    zs_all{i} = zeros(N_Blocks,1)*NaN; %observed roughness height from lowest anemometer
    
    %stress partition values
    eta_inactive_1s_Q0_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, times of zero flux - based on 1s winds, using fQ
    eta_inactive_fQ_1s_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, constant tau_it - based on 1s winds, using fQ
    eta_inactive_fQ_1s_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on 1s winds, using fQ
    eta_inactive_fQ_LP_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, constant tau_it - based on LP winds, using fQ
    eta_inactive_fQ_LP_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on LP winds, using fQ
    eta_active_fQ_1s_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using constant tau_it - based on 1s winds, using fQ
    eta_active_fQ_1s_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using TFEM u_th - based on 1s winds, using fQ
    eta_active_fQ_LP_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using constant tau_it - based on LP winds, using fQ
    eta_active_fQ_LP_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using TFEM u_th - based on LP winds, using fQ
    
    eta_inactive_fU_1s_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, constant tau_it - based on 1s winds, using fU
    eta_inactive_fU_1s_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on 1s winds, using fU
    eta_inactive_fU_LP_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, constant tau_it - based on LP winds, using fU
    eta_inactive_fU_LP_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during inactive transport, TFEM u_th - based on LP winds, using fU
    eta_active_fU_1s_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using constant tau_it - based on 1s winds, using fU
    eta_active_fU_1s_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using TFEM u_th - based on 1s winds, using fU
    eta_active_fU_LP_constthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using constant tau_it - based on LP winds, using fU
    eta_active_fU_LP_TFEMthr_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on mean stress during active transport, using TFEM u_th - based on LP winds, using fU
     
    %inialize arrays of low-pass timeseries
    q_LP_all{i} = cell(N_Blocks,1); %low-pass flux timeseries
    u_LP_all{i} = cell(N_Blocks,1); %low-pass wind timeseries
    t_q_LP_all{i} = (0:dt_q_s:seconds(WindowTimeInterval))';
    t_u_LP_all{i} = (0:dt_u_s:seconds(WindowTimeInterval))';
    
    %initialize arrays of window average timeseries
    Qbinary_1s_all{i} = cell(N_Blocks,1); %1 second wind average flux occurrence timeseries
    u_1s_all{i} = cell(N_Blocks,1); %1 second window average wind velocity timeseries
    t_1s_all{i} = (0:seconds(WindowTimeInterval))'; %times for 1 second window average timeseries
    
%     %initialize lists of profile values
%     u_profile_all{i} = cell(N_Blocks,1); %calibrated velocity profile
%     z_profile_all{i} = cell(N_Blocks,1); %heights for profile
%     Anemometer_profile_all{i} = cell(N_Blocks,1); %names of anemometers for profile
%     u_lowest_all{i} = cell(N_Blocks,1); %calibrated velocity profile - lowest anemometers
%     z_lowest_all{i} = cell(N_Blocks,1); %heights of lowest anemometers
%     Anemometer_lowest_all{i} = cell(N_Blocks,1); %names of lowest anemometers for profile
%     ustLog_all{i} = zeros(N_Blocks,1); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
%     z0_profile_all{i} = zeros(N_Blocks,1); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile
    
    
    %% go through time blocks
    for j = 1:N_Blocks

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Blocks),', ',datestr(now)]
        
        %get specific start and end time
        StartTime = BlockStartTimes(j);
        EndTime = BlockEndTimes(j);

        
        %% GRAIN SIZE INFO
        ind_date = find([Metadata{i}.GrainSize_Surface.Date]==datetime(StartTime.Year, StartTime.Month, StartTime.Day));
        d50_all{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_50_mm]);
        d10_all{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_10_mm]);
        d90_all{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_90_mm]);
        
        
        %% FLUX CALCULATIONS FOR INTERVAL
        
        %extract time interval - get interval number and indices within interval for analysis
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(FluxData,StartTime,EndTime,'t','t','t');
        
        %determine if values span entire analysis interval
        if ~isempty(IntervalN)
            %use only first interval
            IntervalN = IntervalN(1);
            IntervalInd = IntervalInd{1};
            StartTime_Extraction = FluxData(IntervalN).t.t(IntervalInd(1));
            EndTime_Extraction = FluxData(IntervalN).t.t(IntervalInd(end));
                        
            %further reduce IntervalInd based on eliminating error times
            [~, ErrInd, ~] = intersect(FluxData(IntervalN).t.t,FluxData(IntervalN).t.err);
            IntervalInd = setdiff(IntervalInd,ErrInd);
            t_Interval = FluxData(IntervalN).t.t(IntervalInd);
            
        else %generate arbitrary values if there are no data
            StartTime_Extraction = datetime(0,0,0);
            EndTime_Extraction = datetime(0,0,0);
        end
        
        %restrict to time intervals with data spanning entire interval
        if (StartTime_Extraction==StartTime)&&(EndTime_Extraction==EndTime)
            
            %get raw qz and qcal values, determine n from this
            qz = FluxData(IntervalN).qz.qz(IntervalInd,:);
            qcal = FluxData(IntervalN).qz.qzPerCount(IntervalInd,:);
            sigma_qcal = FluxData(IntervalN).qz.sigma_qzPerCount(IntervalInd,:);
            n = FluxData(IntervalN).qz.n(IntervalInd,:); %particle arrival rates
            N = sum(sum(n)); %total detected particles
            
            %get profile information for fitting
            q_profile = mean(qz); %compute mean qz profile
            sigma_q_profile = mean(sigma_qcal.*n); %compute uncertainty in qz profile
            zW_profile = mean(FluxData(IntervalN).qz.z(IntervalInd,:)); %compute mean Wenglor heights
            sigma_zW_profile = mean(FluxData(IntervalN).qz.sigma_z(IntervalInd,:)); %compute uncertainty in Wenglor heights

            %Perform profile fit to get q0, zq, and Q
            [q0,zq,sigma_q0,sigma_zq] = qz_profilefit(q_profile,zW_profile,sigma_q_profile,sigma_zW_profile);
            Q = q0*zq; %get total flux [g/m/s]
            sigma_Q = sqrt((sigma_q0*zq)^2+(sigma_zq*q0)^2); %estimate uncertainty in total flux
                
            %convert to 0 if Q<0.05 or if Q=NaN AND expected Q<0.05 g/m/s based on lowest Wenglor and zq = 5 cm
            if Q<0.05||(isnan(Q)&&((0.1*q_profile(1))/exp(-zW_profile(1)/0.1)<0.05));
                Q=0;
                zq=0;
                sigma_q0=0;
                sigma_zq=0;
            end
            
            %Determine flux frequencies for different window averaging times
            qsum = sum(qz')'; %create vector qsum, which is sum of all q's for each time increment
            if Q==0 %if no flux, set frequencies to zero
                fD_all{i}(j,:) = 0;
                fQ_all{i}(j,:) = 0;
                fQalt_all{i}(j,:) = 0;
            else %otherwise, do calculations
                for k = 1:N_dt_fQ; %go through all window averaging times
                    qsum_avg = window_average(qsum, t_Interval, dt_fQ_windowaverage(k)); %compute window average
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

                    %alternate estimate of fQ
                    fQalt = fD/(1-exp(-N/T)); %calculate new fQ

                    %remove points with fQ>1
                    if fQ>1
                        fQ = NaN;
                    end
                    if fQalt>1
                        fQalt = NaN;
                    end 
                    
                    %add to lists of all fD and fQ
                    fD_all{i}(j,k) = fD;
                    fQ_all{i}(j,k) = fQ;
                    fQalt_all{i}(j,k) = fQalt;
                end
            end
            
            %Timesteps with flux for 1 second wind average
            qsum1_avg = window_average(qsum, t_Interval, duration(0,0,1));
            Qbinary_1s_all{i}{j} = qsum1_avg~=0;
            
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
            fD1_all{i}(j) = fD_all{i}(j,ind_dt_fQ_1s); %1s detection frequency
            fQ1_all{i}(j) = fQ_all{i}(j,ind_dt_fQ_1s); %1s transport frequency
            
            %generate low-pass flux timeseries and add to array
            q_LP = zeros(size(qz));
            for k = 1:length(q_profile)
                q_LP(:,k) = LowPassFilter(qz(:,k), dt_q_s, 0.04);
            end
            q_LP_all{i}{j} = q_LP;
        end

        
        %% WIND CALCULATIONS FOR INTERVAL - BASE ANEMOMETER
        %extract time interval
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindDataBase,StartTime,EndTime,'u','int','int');
        
        %determine if values span entire analysis interval
        if ~isempty(IntervalN)
            %use only first interval
            IntervalN = IntervalN(1);
            IntervalInd = IntervalInd{1};
            StartTime_Extraction = WindDataBase(IntervalN).t.int(IntervalInd(1));
            EndTime_Extraction = WindDataBase(IntervalN).t.int(IntervalInd(end));
            
            %further reduce IntervalInd based on eliminating error times
            [~, ErrInd, ~] = intersect(WindDataBase(IntervalN).t.int,WindDataBase(IntervalN).t.err);
            IntervalInd_NoErr = setdiff(IntervalInd,ErrInd);
            
        else %generate arbitrary values if there are no data
            StartTime_Extraction = datetime(0,0,0);
            EndTime_Extraction = datetime(0,0,0);
        end
        
        %restrict to time intervals with data spanning entire interval
        if (StartTime_Extraction==StartTime)&&(EndTime_Extraction==EndTime)
            
            %get anemometer height
            zU = WindDataBase(IntervalN).z.z;
            
            %get velocity values using no error points
            u = WindDataBase(IntervalN).u.int(IntervalInd_NoErr);
            v = WindDataBase(IntervalN).v.int(IntervalInd_NoErr);
            w = WindDataBase(IntervalN).w.int(IntervalInd_NoErr);
            t = WindDataBase(IntervalN).t.int(IntervalInd_NoErr);
            T_raw = WindDataBase(IntervalN).T.int(IntervalInd_NoErr);

            %compute wind angle
            theta_all{i}(j) = atan(mean(v./u))*180/pi;
            
            %rotate instrument, call these 'raw' values
            [u_raw, ~, w_raw] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument
            
            %get calibration factors
            CalibrationFactor_u = WindDataBase(IntervalN).u.CalibrationFactor;
            CalibrationFactor_w = WindDataBase(IntervalN).w.CalibrationFactor;
            CalibrationFactor_T = WindDataBase(IntervalN).T.CalibrationFactor;
            CalibrationIntercept_u = WindDataBase(IntervalN).u.CalibrationIntercept;
            CalibrationIntercept_w = WindDataBase(IntervalN).w.CalibrationIntercept;
            CalibrationIntercept_T = WindDataBase(IntervalN).T.CalibrationIntercept;
            
            %apply calibration, call these 'cal' values
            u_cal = (u_raw-CalibrationIntercept_u)/CalibrationFactor_u;
            w_cal = (w_raw-CalibrationIntercept_w)/CalibrationFactor_w;
            T_cal = (T_raw-CalibrationIntercept_T)/CalibrationFactor_T;

            %make computations - calibrated values
            u_bar_cal = mean(u_cal);
            w_bar_cal = mean(w_cal);
            T_bar_cal = mean(T_cal);
            ustRe_kernal = mean((u_cal-u_bar_cal).*(w_cal-w_bar_cal));
            if ustRe_kernal<=0
                ustRe_cal = sqrt(-ustRe_kernal);
            else
                ustRe_cal = NaN;
            end
            heatflux_cal = mean((T_cal-T_bar_cal).*(w_cal-w_bar_cal));
            zL = (-(g./(T_bar_cal+273.15)).*heatflux_cal)./(ustRe_cal.^3./(kappa*zU));
            
            %add velocity values to list
            zU_all{i}(j) = zU; %height of anemometer
            u_bar_all{i}(j) = u_bar_cal; %raw mean u values
            u_std_all{i}(j) = std(u_cal); %standard deviation u values
            u2_all{i}(j) = mean(u_cal.^2); %mean u^2
            
            %add stress values to list
            ustRe_all{i}(j) = ustRe_cal;
            tauRe_all{i}(j) = rho_a*ustRe_cal.^2;
            
            %add stability value to list
            zL_all{i}(j) = zL;
            
            %add zs to list
            zs_all{i}(j) = zU*exp(-kappa*u_bar_cal/ustRe_cal); %observed roughness height from lowest anemometer
            
            
            %% create averaged / low-pass wind timeseries, get velocity values including error points
            u_witherror = WindDataBase(IntervalN).u.int(IntervalInd);
            v_witherror = WindDataBase(IntervalN).v.int(IntervalInd);
            w_witherror = WindDataBase(IntervalN).w.int(IntervalInd);
            t_witherror = WindDataBase(IntervalN).t.int(IntervalInd);
            [u_raw_witherror, ~, ~] = reorient_anemometers_vanboxel2004(u_witherror, v_witherror, w_witherror); %rotate instrument
            u_cal_witherror = (u_raw_witherror-CalibrationIntercept_u)/CalibrationFactor_u;
            u_1s = window_average(u_cal_witherror, t, duration(0,0,1)); %compute 1 second window average wind timeseries
            u_LP = LowPassFilter(u_cal_witherror, dt_u_s, 0.04); %compute 0.04 low-pass filtered wind timeseries
            u_LP_all{i}{j} = u_LP; %low-pass wind timeseries - add to array
            u_1s_all{i}{j} = u_1s; %1-second window average timeseries - add to array
            
            
            %% TFEM threshold calcs
            u_it = ust_it(i)/kappa*log(zU/z0(i)); %compute expected u impact threshold based on log law
            u_ft = ust_ft(i)/kappa*log(zU/z0(i)); %compute expected u fluid threshold based on log low
            if fQ_all{i}(j)~=1 %if 1, just use u_it and tau_it
                uth_TFEM = prctile(u_1s,100*(1-fQ1_all{i}(j))); %generate uth from TFEM method based on fQ and 1 s window averaged u
                tauth_TFEM = rho_a*kappa^2*uth_TFEM^2/log(zU/z0(i))^2; %generate tauth from TFEM method based on uth_TFEM 
            else
                uth_TFEM = u_it; %if fQ=1, just use assumed value for uth
                tauth_TFEM = tau_it(i); %if fQ=1, just use assumed value for tauth
            end
            ind_ubelowit = find(u_1s<u_it);
            ind_uaboveit = find(u_1s>=u_it);
            ind_ubelowft = find(u_1s<u_ft);
            ind_uaboveft = find(u_1s>=u_ft);
            f_uaboveft = length(ind_uaboveft)/length(u_1s);
            f_ubelowit = length(ind_ubelowit)/length(u_1s);
            ind_ubetween = intersect(ind_uaboveit,ind_ubelowft);
            N_ubetween = length(ind_ubetween);
            f_ubetween = N_ubetween/length(u_1s);
            N_ubetween_frombelow = 0; %initialize u's in between from below
            N_ubetween_fromabove = 0; %initialize u's in between from above
            N_ubetween_isempty = 0; %initialize u's in between without history available
            for k = 1:N_ubetween;
                ind_last_ubelowit = find(ind_ubelowit<ind_ubetween(k), 1, 'last' );
                ind_last_uaboveft = find(ind_uaboveft<ind_ubetween(k), 1, 'last' );
                if(isempty(ind_last_uaboveft))
                    if(isempty(ind_last_ubelowit))
                        N_ubetween_isempty = N_ubetween_isempty+1;
                    else
                        N_ubetween_frombelow = N_ubetween_frombelow+1;
                    end
                elseif(isempty(ind_last_ubelowit))
                    if(isempty(ind_last_uaboveft))
                        N_ubetween_isempty = N_ubetween_isempty+1;
                    else
                        N_ubetween_fromabove = N_ubetween_fromabove+1;
                    end
                else
                    if(ind_last_ubelowit>ind_last_uaboveft)
                        N_ubetween_frombelow = N_ubetween_frombelow+1;
                    elseif(ind_last_uaboveft>ind_last_ubelowit)
                        N_ubetween_fromabove = N_ubetween_fromabove+1;
                    end
                end
            end
            f_ubetween_frombelow = f_ubetween*N_ubetween_frombelow/(N_ubetween_frombelow+N_ubetween_fromabove); %frequency of u between thresholds approaching from below it
            f_ubetween_fromabove = f_ubetween*N_ubetween_fromabove/(N_ubetween_frombelow+N_ubetween_fromabove); %frequency of u between thresholds approaching from above ft
                    
            %add thresholds to lists
            uth_TFEM_all{i}(j) = uth_TFEM; %uth from TFEM
            tauth_TFEM_all{i}(j) = tauth_TFEM; %tauth from TFEM
            f_uaboveft_all{i}(j) = f_uaboveft; %frequency of u above fluid threshold
            f_ubelowit_all{i}(j) = f_ubelowit; %frequency of u below impact threshold
            f_ubetween_all{i}(j) = f_ubetween; %frequency of u between thresholds
            f_ubetween_frombelow_all{i}(j) = f_ubetween_frombelow; %frequency of u between thresholds approaching from below it
            f_ubetween_fromabove_all{i}(j) = f_ubetween_fromabove; %frequency of u between thresholds approaching from above ft
            
            %% stress partition calculations - based on mean stress during inactive transport
            % Based on mean stress of times with no transport detection - u_1s
            u2bar_inactive_Q0 = mean(u_1s(Qbinary_1s_all{i}{j}==0).^2); %get mean u^2 during no transport based on 1 second flux/wind timeseries
            tau_inactive_Q0 = rho_a*kappa^2*u2bar_inactive_Q0/log(zU/z0(i))^2; %get mean stress during no transport based on 1 second flux timeseries
            tau0_bar_Q0 = fQ1_all{i}(j)*tau_it(i)+(1-fQ1_all{i}(j))*tau_inactive_Q0; %get mean stress dissipation
            eta_inactive_1s_Q0_all{i}(j) = max([1-tau0_bar_Q0/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, times of zero flux
            
            % Based on mean stress when wind is below impact threshold - u_1s
            u2bar_inactive = mean(u_1s(u_1s<u_it).^2); %get mean u^2 during no transport based on u_1s winds less than uth for fixed tauit
            tau_inactive = rho_a*kappa^2*u2bar_inactive/log(zU/z0(i))^2; %get mean stress during no transport based on 1 second winds less than uth for fixed tauit
            tau0_bar = fQ1_all{i}(j)*tau_it(i)+(1-fQ1_all{i}(j))*tau_inactive; %get mean stress dissipation
            eta_inactive_fQ_1s_constthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            % redo calculation using fU instead of fQ
            fU = length(find(u_1s>=u_it))/length(u_1s);
            tau0_bar = fU*tau_it(i)+(1-fU)*tau_inactive; %get mean stress dissipation
            eta_inactive_fU_1s_constthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            
            % Based on mean stress when wind is below TFEM threshold - u_1s
            u2bar_inactive = mean(u_1s(u_1s<uth_TFEM).^2); %get mean u^2 during no transport based on 1 second winds less than uth from TFEM
            tau_inactive = rho_a*kappa^2*u2bar_inactive/log(zU/z0(i))^2; %get mean stress during no transport based on 1 second winds less than uth from TFEM
            tau0_bar = fQ1_all{i}(j)*tau_it(i)+(1-fQ1_all{i}(j))*tau_inactive; %get mean stress dissipation
            eta_inactive_fQ_1s_TFEMthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            % redo calculation using fU instead of fQ
            fU = length(find(u_1s>=uth_TFEM))/length(u_1s);
            tau0_bar = fU*tau_it(i)+(1-fU)*tau_inactive; %get mean stress dissipation
            eta_inactive_fU_1s_TFEMthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            
            % Based on mean stress when wind is below impact threshold - u_LP 
            u2bar_inactive = mean(u_LP(u_LP<u_it).^2); %get mean u^2 during no transport based on u_LP winds less than uth for fixed tauit
            tau_inactive = rho_a*kappa^2*u2bar_inactive/log(zU/z0(i))^2; %get mean stress during no transport based on 1 second winds less than uth for fixed tauit
            tau0_bar = fQ1_all{i}(j)*tau_it(i)+(1-fQ1_all{i}(j))*tau_inactive; %get mean stress dissipation
            eta_inactive_fQ_LP_constthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            % redo calculation using fU instead of fQ
            fU = length(find(u_LP>=u_it))/length(u_LP);
            tau0_bar = fU*tau_it(i)+(1-fU)*tau_inactive; %get mean stress dissipation
            eta_inactive_fU_LP_constthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            
            % Based on mean stress when wind is below TFEM threshold - u_LP 
            u2bar_inactive = mean(u_LP(u_LP<uth_TFEM).^2); %get mean u^2 during no transport based on 1 second winds less than uth from TFEM
            tau_inactive = rho_a*kappa^2*u2bar_inactive/log(zU/z0(i))^2; %get mean stress during no transport based on 1 second winds less than uth from TFEM
            tau0_bar = fQ1_all{i}(j)*tau_it(i)+(1-fQ1_all{i}(j))*tau_inactive; %get mean stress dissipation
            eta_inactive_fQ_LP_TFEMthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            % redo calculation using fU instead of fQ
            fU = length(find(u_LP>=u_it))/length(u_LP);
            tau0_bar = fU*tau_it(i)+(1-fU)*tau_inactive; %get mean stress dissipation
            eta_inactive_fU_LP_constthr_all{i}(j) = max([1-tau0_bar/tauRe_all{i}(j), 0]); %mean stress partition based on flux frequency method, constant tau_it
            
            %% Stress partition calculations -- based on mean stress during active transport
            %generate zs(u) rating curve, compute zs for each timestep
            u_ratingcurve = linspace(min(u_1s),max(u_1s),N_zs_ratingcurve)'; %initialize list of u values for rating curve
            zs_ratingcurve = linspace(min(u_1s),max(u_1s),N_zs_ratingcurve)'; %initialize list of corresponding zs(u) values for rating curve
            for k = 1:N_zs_ratingcurve 
                if u_ratingcurve(k)<u_it %if below threshold, zs=z0
                    zs_ratingcurve(k) = z0(i);
                else %if above threshold, iteratively determine zs and tau
                    zs_old = z0(i); %initialize zs guess based on z0
                    tau_old = rho_a*kappa^2*u_ratingcurve(k).^2/log(zU/z0(i))^2; %initialize tau guess based on z0
                    zs_new = z0(i)+k_zs*(tau_old-tau_it(i)); %revise zs guess based on tau
                    tau_new = rho_a*kappa^2*u_ratingcurve(k).^2/log(zU/zs_new)^2; %compute new tau based on zs
                    while abs(tau_new-tau_old)>0.0001 %if new and old taus are too far apart
                        zs_old = zs_new; %set the "new" guesses to the "old" guesses
                        tau_old = tau_new; %set the "new" guesses to the "old" guesses
                        zs_new = z0(i)+k_zs*(tau_old-tau_it(i)); %recompute zs
                        tau_new = rho_a*kappa^2*u_ratingcurve(k).^2/log(zU/zs_new)^2; %recompute tau
                    end
                    zs_ratingcurve(k) = zs_new; %add roughness length to rating curve list
                end
            end
            
            %get zs values corresponding to u timeseries
            zs_1s = zeros(size(u_1s));
            for k = 1:length(u_1s)
                zs_1s(k) = zs_ratingcurve(abs(u_ratingcurve-u_1s(k))==min(abs(u_ratingcurve-u_1s(k))));
            end
            zs_LP = zeros(size(u_LP));
            for k = 1:length(u_LP)
                zs_LP(k) = zs_ratingcurve(abs(u_ratingcurve-u_LP(k))==min(abs(u_ratingcurve-u_LP(k))));
            end
            
            %compute eta value on u and zs, for u_1s>u_it
            tau_inst = rho_a*kappa^2*u_1s.^2./log(zU./zs_1s).^2;
            eta_inst = zeros(size(u_1s));
            ind_abovethr = find(u_1s>=u_it);
            eta_inst(ind_abovethr)=1-(tau_it(i)./tau_inst(ind_abovethr));
            eta_active_fQ_1s_constthr_all{i}(j) = fQ1_all{i}(j)*mean(eta_inst(ind_abovethr)); %mean stress partition for active transport, constant tau_it
            eta_active_fU_1s_constthr_all{i}(j) = mean(eta_inst); %same as before, but now use the actual frequency of wind above threshold
            
            %compute eta value on u and zs, for u_1s>u_th (TFEM threshold)
            tau_inst = rho_a*kappa^2*u_1s.^2./log(zU./zs_1s).^2;
            eta_inst = zeros(size(u_1s));
            ind_abovethr = find(u_1s>=uth_TFEM);
            eta_inst(ind_abovethr)=1-(tau_it(i)./tau_inst(ind_abovethr));
            eta_active_fQ_1s_TFEMthr_all{i}(j) = fQ1_all{i}(j)*mean(eta_inst(ind_abovethr)); %mean stress partition for active transport, TFEM threshold
            eta_active_fU_1s_TFEMthr_all{i}(j) = mean(eta_inst); %same as before, but now use the actual frequency of wind above threshold
            
            %compute eta value on u and zs, for u_LP>u_it
            tau_inst = rho_a*kappa^2*u_LP.^2./log(zU./zs_LP).^2;
            eta_inst = zeros(size(u_LP));
            ind_abovethr = find(u_LP>=u_it);
            eta_inst(ind_abovethr)=1-(tau_it(i)./tau_inst(ind_abovethr));
            eta_active_fQ_LP_constthr_all{i}(j) = fQ1_all{i}(j)*mean(eta_inst(ind_abovethr)); %mean stress partition for active transport, constant tau_it
            eta_active_fU_LP_constthr_all{i}(j) = mean(eta_inst); %same as before, but now use the actual frequency of wind above threshold
            
            %compute eta value on u and zs, for u_LP>u_th (TFEM threshold)
            tau_inst = rho_a*kappa^2*u_LP.^2./log(zU./zs_LP).^2;
            eta_inst = zeros(size(u_LP));
            ind_abovethr = find(u_LP>=uth_TFEM);
            eta_inst(ind_abovethr)=1-(tau_it(i)./tau_inst(ind_abovethr));
            eta_active_fQ_LP_TFEMthr_all{i}(j) = fQ1_all{i}(j)*mean(eta_inst(ind_abovethr)); %mean stress partition for active transport, TFEM threshold
            eta_active_fU_LP_TFEMthr_all{i}(j) = mean(eta_inst); %same as before, but now use the actual frequency of wind above threshold
        end
%          
%         %% WIND CALCULATIONS FOR INTERVAL - FULL PROFILE
%         u_cal_profile = zeros(N_Anemometers,1)*NaN;
%         z_profile = zeros(N_Anemometers,1)*NaN;
%         
%         for k = 1:N_Anemometers;
%             %extract time interval
%             [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindDataAll.(Anemometer_profile{k}),StartTime,EndTime,'u','int','int');
% 
%             %determine if values span entire analysis interval
%             if ~isempty(IntervalN)
%                 %use only first interval
%                 IntervalN = IntervalN(1);
%                 IntervalInd = IntervalInd{1};
%                 StartTime_Extraction = WindDataAll.(Anemometer_profile{k})(IntervalN).t.int(IntervalInd(1));
%                 EndTime_Extraction = WindDataAll.(Anemometer_profile{k})(IntervalN).t.int(IntervalInd(end));
% 
%                 %further reduce IntervalInd based on eliminating error times
%                 [~, ErrInd, ~] = intersect(WindDataAll.(Anemometer_profile{k})(IntervalN).t.int,WindDataAll.(Anemometer_profile{k})(IntervalN).t.err);
%                 IntervalInd = setdiff(IntervalInd,ErrInd);
% 
%             else %generate arbitrary values if there are no data
%                 StartTime_Extraction = datetime(0,0,0);
%                 EndTime_Extraction = datetime(0,0,0);
%             end
% 
%             %restrict to time intervals with data spanning entire interval
%             if (StartTime_Extraction==StartTime)&&(EndTime_Extraction==EndTime)
%                 
%                 %get velocity values
%                 u = WindDataAll.(Anemometer_profile{k})(IntervalN).u.int(IntervalInd);
%                 v = WindDataAll.(Anemometer_profile{k})(IntervalN).v.int(IntervalInd);
%                 w = WindDataAll.(Anemometer_profile{k})(IntervalN).w.int(IntervalInd);
% 
%                 %rotate instrument, call these 'raw' values
%                 [u_raw, ~, ~] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument
% 
%                 %get calibration factors
%                 CalibrationFactor_u = WindDataAll.(Anemometer_profile{k})(IntervalN).u.CalibrationFactor;
%                 CalibrationIntercept_u = WindDataAll.(Anemometer_profile{k})(IntervalN).u.CalibrationIntercept;
% 
%                 %apply calibration, call these 'cal' values
%                 u_cal = (u_raw-CalibrationIntercept_u)/CalibrationFactor_u;
% 
%                 %add mean velocities to profile
%                 u_cal_profile(k) = mean(u_cal);
% 
%                 %add height to profile
%                 z_profile(k) = WindDataAll.(Anemometer_profile{k})(IntervalN).z.z;
%             end
%         end
%         
%         %get indices of NaN entries in profile
%         ind_notnan = find(~isnan(u_cal_profile));
%         
%         %if some profile exists...
%         if ~isempty(ind_notnan)
%             
%             %remove NaN's from profile, add to 'profile_all' lists
%             u_cal_profile = u_cal_profile(ind_notnan);
%             z_profile = z_profile(ind_notnan);
%             
%             %add to lists
%             u_profile_all{i}{j} = u_cal_profile;
%             z_profile_all{i}{j} = z_profile;
%             Anemometer_profile_all{i}{j} = Anemometer_profile{ind_notnan};
%             
%             %get profile for three lowest anemometers
%             [~, ind_sort] = sort(z_profile);
%         
%             %check to make sure there are enough values for fitting log
%             if length(ind_sort)>=3
%                 ind_lowest = ind_sort(1:3);
%                 u_cal_lowest = u_cal_profile(ind_lowest);
%                 z_lowest = z_profile(ind_lowest);
% 
%                 %compute u* and z0, add to list
%                 P_cal = polyfit(log(z_lowest),u_cal_lowest,1);
%                 ustLog_cal = kappa*P_cal(1);
%                 z0_cal = exp(-P_cal(2)/P_cal(1));
% 
%                 %deal with poorly fit profiles
%                 if(ustLog_cal<0||isnan(ustLog_cal))
%                     ustLog_cal = NaN;
%                     z0_cal = NaN;
%                 end
%             
%                 %add to lists
%                 u_lowest_all{i}{j} = u_cal_lowest;
%                 z_lowest_all{i}{j} = z_lowest;
%                 Anemometer_lowest_all{i}{j} = Anemometer_profile{ind_lowest};
%                 ustLog_all{i}(j) = ustLog_cal;
%                 z0_profile_all{i}(j) = z0_cal;
%             end
%         end 
    end
end

% SAVE DATA
save(SaveFullData_Path,'Sites','*all','-v7.3'); %save full size file including all sub-timeseries
clear('*LP_all','*1s_all'); %remove timeseries to make file smaller
save(SaveData_Path,'Sites','*all'); %save reduced file in GitHub folder