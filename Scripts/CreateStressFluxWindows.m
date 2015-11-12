%% SCRIPT TO GENERATE WINDOWS OF SALTATION FLUX AND STRESS VALUES FOR ANALYSIS
% SCRIPT DEPENDENCIES:
% CreateStressFluxWindows.m
% CreateTimeBlocks.m
% ExtractVariableTimeInterval.m
% IntersectingTimeIntervals.m
% qz_profilefit
% reorient_anemometers_vanboxel2004
% window_average

%% initialize
clearvars;

%% parameter values
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
g = 9.8; %gravity m/s^2
kappa = 0.39; %von Karman parameter

%% information about sites for analysis
%Sites = {'Jericoacoara'};
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%threshold information by site
ust_it = [0.35, 0.28, 0.28]; %shear velocity (m/s) threshold
z0f = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m) at threshold
k_zs = 0.004; %proportionality constant for zs (m^2 s^2 kg^-1)

%z0f = [9e-5, 1e-4, 2e-4]; %aerodynamic roughness length (m) at threshold
%z0f = [1e-5, 1e-5, 1e-5]; %aerodynamic roughness length (m) at threshold
tau_it = rho_a*ust_it.^2; %shear stress (Pa) threshold

%% set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
N_RunningPerProfile = floor(ProfileTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times

%% set durations of window average for flux frequencies
dt_min_s = 1;
dt_max_s = 1;
N_dt = 1;
% dt_min_s = 0.04; %mininum window average time
% dt_max_s = 600; %maximum window average time
% N_dt = 30; %number of durations of window average
dt_s = unique([round(logspace(0,log10(dt_max_s/dt_min_s),N_dt))*dt_min_s, 1]); %create window average dt's, make sure to include 1 second in list
N_dt = length(dt_s); %recalculate number of window average dt's after removing repeats
dt_windowaverage = duration(0,0,dt_s);
ind_dt_1s = find(dt_s==1); %index in list of dt's for 1 second

%set number of points to use for rating curve of Cf versus u
N_Cf_ratingcurve = 1000;

%% information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../AnalysisData/'; %folder for storing outputs of this analysis
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

%initialize lists of flux values
Q_all = cell(N_Sites,1); %total flux
zq_all = cell(N_Sites,1); %characteristic flux height
fD_all = cell(N_Sites,1); %Wenglor detection frequency matrix
fD1_all = cell(N_Sites,1); %Wenglor detection frequency list - 1 second
fQ_all = cell(N_Sites,1); %Wenglor transport frequency matrix
fQ1_all = cell(N_Sites,1); %Wenglor transport frequency matrix - 1 second
fQalt_all = cell(N_Sites,1); %Wenglor transport frequency matrix - alternate (non-iterative) calc
Q1binary_all = cell(N_Sites,1); %Timesteps with flux for 1 second wind average
q_all = cell(N_Sites,1); %partial flux
zW_all = cell(N_Sites,1); %Wenglor flux height list
qcal_all = cell(N_Sites,1); %calibration factor list

%initialize lists of wind values
u_bar_raw_all = cell(N_Sites,1); %raw mean u values
u_std_raw_all = cell(N_Sites,1); %standard deviation u values
u2_raw_all = cell(N_Sites,1); %mean u^2
u_bar_cal_all = cell(N_Sites,1); %raw mean u values
u_std_cal_all = cell(N_Sites,1); %standard deviation u values
u2_cal_all = cell(N_Sites,1); %mean u^2

%initialize lists of stress values
ustRe_raw_all = cell(N_Sites,1); %raw u* for Reynolds stress
tauRe_raw_all = cell(N_Sites,1); %raw tau for Reynolds stress
ustRe_cal_all = cell(N_Sites,1); %calibrated u* for Reynolds stress
tauRe_cal_all = cell(N_Sites,1); %calibrated tau for Reynolds stress

%initialize lists of eta values
eta_fQ_Q0_all = cell(N_Sites,1); %mean stress partition based on flux frequency method, times of zero flux
eta_fQ_tauit_all = cell(N_Sites,1); %mean stress partition based on flux frequency method, constant tau_it
eta_fQ_TFEM_all = cell(N_Sites,1); %mean stress partition based on flux frequency method, TFEM u_th
eta_zs_all = cell(N_Sites,1); %mean stress partition based on friction coefficient estimation from z_s

%initialize lists of threshold values
uth_TFEM_all = cell(N_Sites,1); %1-second TFEM estimate of uth

%initialize lists of grain sizes, dates, start times, and end times
d50_all = cell(N_Sites,1); %lists of surface grain sizes corresponding to calculated values
d10_all = cell(N_Sites,1);
d90_all = cell(N_Sites,1);
date_all = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_all = cell(N_Sites,1); %lists of block start times
EndTime_all = cell(N_Sites,1); %lists of block end times

% %initialize lists of profile values
% u_raw_profile_all = cell(N_Sites,1); %raw velocity profile
% u_cal_profile_all = cell(N_Sites,1); %calibrated velocity profile
% z_profile_all = cell(N_Sites,1); %heights for profile
% Anemometer_profile_all = cell(N_Sites,1); %names of anemometers for profile
% u_raw_lowest_all = cell(N_Sites,1); %raw velocity profile - lowest anemometers
% u_cal_lowest_all = cell(N_Sites,1); %calibrated velocity profile - lowest anemometers
% z_lowest_all = cell(N_Sites,1); %heights of lowest anemometers
% Anemometer_lowest_all = cell(N_Sites,1); %names of lowest anemometers for profile
% ustLog_raw_all = cell(N_Sites,1); %u* calculated from mean of raw velocities for lowest anemometers in profile
% z0_raw_all = cell(N_Sites,1); %z0 calculated from mean of raw velocities for lowest anemometers in profile
% ustLog_cal_all = cell(N_Sites,1); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
% z0_cal_all = cell(N_Sites,1); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
    
    %% choose anemometer type based on site of interest
    if strcmp(Sites{i},'Oceano')
        AnemometerType = 'Sonic';
        BaseAnemometer = 'S1';
    elseif strcmp(Sites{i},'RanchoGuadalupe')
        AnemometerType = 'Ultrasonic';
        BaseAnemometer = 'U1';
    elseif strcmp(Sites{i},'Jericoacoara')
        AnemometerType = 'Ultrasonic';
        BaseAnemometer = 'U1';
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
            CreateTimeBlocks(StartTimesIntersecting, EndTimesIntersecting, ProfileTimeInterval);
    
    %Add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [BlockStartTimes_j, BlockEndTimes_j] = ...
            CreateTimeBlocks(StartTimesIntersecting+(j*RunningTimeInterval), EndTimesIntersecting, ProfileTimeInterval);
        BlockStartTimes = [BlockStartTimes; BlockStartTimes_j];
        BlockEndTimes = [BlockEndTimes; BlockEndTimes_j];
    end
    BlockStartTimes = sort(BlockStartTimes);
    BlockEndTimes = sort(BlockEndTimes);
    N_Blocks = length(BlockStartTimes);
    
    %% initialize lists of flux values
    Q_all{i} = zeros(N_Blocks,1)*NaN; %total flux
    zq_all{i} = zeros(N_Blocks,1)*NaN; %characteristic flux height
    fD_all{i} = zeros(N_Blocks,N_dt)*NaN; %Wenglor detection frequency
    fD1_all{i} = zeros(N_Blocks,1)*NaN; %Wenglor detection frequency - 1 s
    fQ_all{i} = zeros(N_Blocks,N_dt)*NaN; %Wenglor total transport frequency
    fQ1_all{i} = zeros(N_Blocks,1)*NaN; %Wenglor total transport frequency - 1 s
    fQalt_all{i} = zeros(N_Blocks,N_dt)*NaN; %Wenglor total transport frequency, alternate calc
    Q1binary_all{i} = cell(N_Blocks,1); %timesteps with flux for 1 second wind average
    q_all{i} = cell(N_Blocks,1); %partial flux
    zW_all{i} = cell(N_Blocks,1); %Wenglor height
    qcal_all{i} = cell(N_Blocks,1); %calibration factor
    
    %initialize lists of wind values
    u_bar_raw_all{i} = zeros(N_Blocks,1)*NaN; %raw mean u values
    u_std_raw_all{i} = zeros(N_Blocks,1)*NaN; %standard deviation u values
    u2_raw_all{i} = zeros(N_Blocks,1)*NaN; %mean u^2
    u_bar_cal_all{i} = zeros(N_Blocks,1)*NaN; %raw mean u values
    u_std_cal_all{i} = zeros(N_Blocks,1)*NaN; %standard deviation u values
    u2_cal_all{i} = zeros(N_Blocks,1)*NaN; %mean u^2
    
    %initialize lists of thresholds
    uth_TFEM_all{i} = zeros(N_Blocks,1)*NaN; %1 second TFEM calc of threshold
    
    %initialize lists of stress values
    ustRe_raw_all{i} = zeros(N_Blocks,1)*NaN; %raw u* for Reynolds stress
    tauRe_raw_all{i} = zeros(N_Blocks,1)*NaN; %raw tau for Reynolds stress
    ustRe_cal_all{i} = zeros(N_Blocks,1)*NaN; %calibrated u* for Reynolds stress
    tauRe_cal_all{i} = zeros(N_Blocks,1)*NaN; %calibrated tau for Reynolds stress
    
    %initialize lists of eta values
    eta_fQ_Q0_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on flux frequency method, times of zero flux
    eta_fQ_tauit_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on flux frequency method, constant tau_it
    eta_fQ_TFEM_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on flux frequency method, TFEM u_th
    eta_zs_all{i} = zeros(N_Blocks,1)*NaN; %mean stress partition based on friction coefficient estimation from z_s, constant tau_it

    %create list of dates
    date_all{i} = datetime(BlockStartTimes.Year, BlockStartTimes.Month, BlockStartTimes.Day); %lists of dates corresponding to calculated values
    
    %initialize lists of grain sizes
    d50_all{i} = zeros(N_Blocks,1); %lists of surface grain sizes corresponding to calculated values
    d10_all{i} = zeros(N_Blocks,1);
    d90_all{i} = zeros(N_Blocks,1);
    
%     %initialize lists of profile values
%     u_raw_profile_all{i} = cell(N_Blocks,1); %raw velocity profile
%     u_cal_profile_all{i} = cell(N_Blocks,1); %calibrated velocity profile
%     z_profile_all{i} = cell(N_Blocks,1); %heights for profile
%     Anemometer_profile_all{i} = cell(N_Blocks,1); %names of anemometers for profile
%     u_raw_lowest_all{i} = cell(N_Blocks,1); %raw velocity profile - lowest anemometers
%     u_cal_lowest_all{i} = cell(N_Blocks,1); %calibrated velocity profile - lowest anemometers
%     z_lowest_all{i} = cell(N_Blocks,1); %heights of lowest anemometers
%     Anemometer_lowest_all{i} = cell(N_Blocks,1); %names of lowest anemometers for profile
%     ustLog_raw_all{i} = zeros(N_Blocks,1); %u* calculated from mean of raw velocities for lowest anemometers in profile
%     z0_raw_all{i} = zeros(N_Blocks,1); %z0 calculated from mean of raw velocities for lowest anemometers in profile
%     ustLog_cal_all{i} = zeros(N_Blocks,1); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
%     z0_cal_all{i} = zeros(N_Blocks,1); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile
    
    %% go through time blocks
    for j = 1:N_Blocks

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Blocks),', ',datestr(now)]
        
        %get specific start and end time
        StartTime = BlockStartTimes(j);
        EndTime = BlockEndTimes(j);
        
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
            n = qz./qcal; %particle arrival rates
            N = sum(sum(n)); %total detected particles
            
            %compute mean qz profile
            q_profile = mean(qz);
            zW_profile = FluxData(IntervalN).z.z;

            %Perform profile fit to get q0, zq, and Q
            [q0,zq] = qz_profilefit(q_profile,zW_profile,0,0); %For now, ignore uncertainty in profile values
            Q = q0*zq; %get total flux [g/m/s]
            
            %convert to 0 if NaN and expected Q<0.05 g/m/s based on lowest Wenglor and ze = 5 cm
            if isnan(Q)&&((0.1*q_profile(1))/exp(-zW_profile(1)/0.1)<0.05);
                Q=0;
                zq=0;
            end
            
            %Determine flux frequencies for different window averaging times
            qsum = sum(qz')'; %create vector qsum, which is sum of all q's for each time increment
            if Q==0 %if no flux, set frequencies to zero
                fD_all{i}(j,:) = 0;
                fQ_all{i}(j,:) = 0;
                fQalt_all{i}(j,:) = 0;
            else %otherwise, do calculations
                for k = 1:N_dt; %go through all window averaging times
                    qsum_avg = window_average(qsum, t_Interval, dt_windowaverage(k)); %compute window average
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
            Q1binary_all{i}{j} = qsum1_avg~=0;
            
            %get calibration values
            qcal_profile = mean(qcal);
            
            %add to list
            Q_all{i}(j) = Q; %total flux
            zq_all{i}(j) = zq; %characteristic flux height
            q_all{i}{j} = q_profile; %partial flux
            zW_all{i}{j} = zW_profile; %Wenglor flux height
            qcal_all{i}{j} = qcal_profile; %calibration factor
            fD1_all{i}(j) = fD_all{i}(j,ind_dt_1s); %1s detection frequency
            fQ1_all{i}(j) = fQ_all{i}(j,ind_dt_1s); %1s transport frequency
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
            
            %rotate instrument, call these 'raw' values
            [u_raw, ~, w_raw] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

            %get calibration factors
            CalibrationFactor_u = WindDataBase(IntervalN).u.CalibrationFactor;
            CalibrationFactor_w = WindDataBase(IntervalN).w.CalibrationFactor;
            CalibrationIntercept_u = WindDataBase(IntervalN).u.CalibrationIntercept;
            CalibrationIntercept_w = WindDataBase(IntervalN).w.CalibrationIntercept;

            %apply calibration, call these 'cal' values
            u_cal = (u_raw-CalibrationIntercept_u)/CalibrationFactor_u;
            w_cal = (w_raw-CalibrationIntercept_w)/CalibrationFactor_w;

            %make computations - raw values
            u_bar_raw = mean(u_raw);
            w_bar_raw = mean(w_raw);
            ustRe_kernal = mean((u_raw-u_bar_raw).*(w_raw-w_bar_raw));
            if ustRe_kernal<=0
                ustRe_raw = sqrt(-ustRe_kernal);
            else
                ustRe_raw = NaN;
            end

            %make computations - calibrated values
            u_bar_cal = mean(u_cal);
            w_bar_cal = mean(w_cal);
            ustRe_kernal = mean((u_cal-u_bar_cal).*(w_cal-w_bar_cal));
            if ustRe_kernal<=0
                ustRe_cal = sqrt(-ustRe_kernal);
            else
                ustRe_cal = NaN;
            end
            
            %add velocity values to list
            u_bar_raw_all{i}(j) = u_bar_raw; %raw mean u values
            u_std_raw_all{i}(j) = std(u_raw); %standard deviation u values
            u2_raw_all{i}(j) = mean(u_raw.^2); %mean u^2
            u_bar_cal_all{i}(j) = u_bar_cal; %raw mean u values
            u_std_cal_all{i}(j) = std(u_cal); %standard deviation u values
            u2_cal_all{i}(j) = mean(u_cal.^2); %mean u^2
            
            %add stress values
            ustRe_raw_all{i}(j) = ustRe_raw;
            tauRe_raw_all{i}(j) = rho_a*ustRe_raw.^2;
            ustRe_cal_all{i}(j) = ustRe_cal;
            tauRe_cal_all{i}(j) = rho_a*ustRe_cal.^2;
            
            %get velocity values including error points
            u_witherror = WindDataBase(IntervalN).u.int(IntervalInd);
            v_witherror = WindDataBase(IntervalN).v.int(IntervalInd);
            w_witherror = WindDataBase(IntervalN).w.int(IntervalInd);
            t_witherror = WindDataBase(IntervalN).t.int(IntervalInd);
            [u_raw_witherror, ~, ~] = reorient_anemometers_vanboxel2004(u_witherror, v_witherror, w_witherror); %rotate instrument
            u_cal_witherror = (u_raw_witherror-CalibrationIntercept_u)/CalibrationFactor_u;
            
            %stress partition calculations
            u_1s = window_average(u_cal_witherror, t, duration(0,0,1)); %compute 1 second window average wind timeseries
            u_LP = LowPassFilter(u_cal_witherror, seconds(mode(diff(t))), 0.04); %compute 0.04 low-pass filtered wind timeseries
            uth_tauit = ust_it(i)/kappa*log(zU/z0f(i)); %compute expected threshold based on log law
            uth_TFEM = prctile(u_1s,100*(1-fQ1_all{i}(j))); %generate uth from TFEM method based on fQ and 1-s window averaged u
            tauth_TFEM = rho_a*kappa^2*uth_TFEM^2/log(zU/z0f(i))^2; %generate tauth from TFEM method based on uth_TFEM 
            ubarQ0_Q0 = mean(u_1s(Q1binary_all{i}{j}==0)); %get mean u during no transport based on 1 second flux/wind timeseries
            tauQ0_Q0 = rho_a*kappa^2*ubarQ0_Q0^2/log(zU/z0f(i))^2; %get mean stress during no transport based on 1 second flux timeseries
            ubarQ0_tauit = mean(u_1s(u_1s<uth_tauit)); %get mean u during no transport based on 1 second winds less than uth for fixed tauit
            tauQ0_tauit = rho_a*kappa^2*ubarQ0_tauit^2/log(zU/z0f(i))^2; %get mean stress during no transport based on 1 second winds less than uth for fixed tauit
            ubarQ0_TFEM = mean(u_1s(u_1s<uth_TFEM)); %get mean u during no transport based on 1 second winds less than uth from TFEM
            tauQ0_TFEM = rho_a*kappa^2*ubarQ0_TFEM^2/log(zU/z0f(i))^2; %get mean stress during no transport based on 1 second winds less than uth from TFEM

            %generate Cf(u) rating curve, compute Cf for each timestep
            u_ratingcurve = linspace(min(u_LP),max(u_LP),N_Cf_ratingcurve)'; %generate list of u values for rating curve
            Cf_ratingcurve = zeros(N_Cf_ratingcurve,1); %initialize list of associated Cf values for rating curve
            for k = 1:N_Cf_ratingcurve 
                if u_ratingcurve(k)<uth_tauit %if below threshold, Cf = 1
                    Cf_ratingcurve(k) = 1;
                else %if above threshold, iteratively determine zs and tau
                    zs_old = z0f(i); %initialize zs guess based on z0
                    tau_old = rho_a*kappa^2*u_ratingcurve(k).^2/log(zU/z0f(i))^2; %initialize tau guess based on z0
                    zs_new = z0f(i)+k_zs*(tau_old-tau_it(i)); %revise zs guess based on tau
                    tau_new = rho_a*kappa^2*u_ratingcurve(k).^2/log(zU/zs_new)^2; %compute new tau based on zs
                    while abs(tau_new-tau_old)>0.0001 %if new and old taus are too far apart
                        zs_old = zs_new; %set the "new" guesses to the "old" guesses
                        tau_old = tau_new; %set the "new" guesses to the "old" guesses
                        zs_new = z0f(i)+k_zs*(tau_old-tau_it(i)); %recompute zs
                        tau_new = rho_a*kappa^2*u_ratingcurve(k).^2/log(zU/zs_new)^2; %recompute tau
                    end
                    Cf_ratingcurve(k) = log(zU/zs_new)^2/log(zU/z0f(i))^2; %compute friction coefficient
                end
            end
            Cf = zeros(length(u_LP),1);
            for k = 1:length(u_LP)
                Cf(k) = Cf_ratingcurve(abs(u_ratingcurve-u_LP(k))==min(abs(u_ratingcurve-u_LP(k))));
            end
            %compute instantaneous eta value based on 1-Cf*(uth^2/uLP^2), for uLP>uth
            eta_inst = zeros(size(u_LP));
            ind_abovethr = find(u_LP>uth_tauit);
            eta_inst(u_LP>uth_tauit)=1-Cf(ind_abovethr).*(uth_tauit^2./u_LP(ind_abovethr).^2);
            
            %calculate various versions of eta
            eta_fQ_Q0_all{i}(j) = (tauRe_cal_all{i}(j)-fQ1_all{i}(j)*tau_it(i)-(1-fQ1_all{i}(j))*tauQ0_Q0)/tauRe_cal_all{i}(j); %mean stress partition based on flux frequency method, times of zero flux
            eta_fQ_tauit_all{i}(j) = (tauRe_cal_all{i}(j)-fQ1_all{i}(j)*tau_it(i)-(1-fQ1_all{i}(j))*tauQ0_tauit)/tauRe_cal_all{i}(j); %mean stress partition based on flux frequency method, constant tau_it
            eta_fQ_TFEM_all{i}(j) = (tauRe_cal_all{i}(j)-fQ1_all{i}(j)*tau_it(i)-(1-fQ1_all{i}(j))*tauQ0_TFEM)/tauRe_cal_all{i}(j); %mean stress partition based on flux frequency method, TFEM u_th
            eta_zs_all{i}(j) = mean(eta_inst); %mean stress partition based on friction coefficient estimation from z_s, constant tau_it
            
            %add thresholds to lists
            uth_TFEM_all{i}(j) = uth_TFEM; %uth from TFEM
        end
        
        %GRAIN SIZE INFO
        ind_date = find([Metadata{i}.GrainSize_Surface.Date]==WindDataBase(IntervalN).Date);
        d50_all{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_50_mm]);
        d10_all{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_10_mm]);
        d90_all{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_90_mm]);

%         %% WIND CALCULATIONS FOR INTERVAL - FULL PROFILE
%         u_raw_profile = zeros(N_Anemometers,1)*NaN;
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
%                 u_raw_profile(k) = mean(u_raw);
%                 u_cal_profile(k) = mean(u_cal);
% 
%                 %add height to profile
%                 z_profile(k) = WindDataAll.(Anemometer_profile{k})(IntervalN).z.z;
%             end
%         end
%         
%         %get indices of NaN entries in profile
%         ind_notnan = find(~isnan(u_raw_profile));
%         
%         %if some profile exists...
%         if ~isempty(ind_notnan)
%             
%             %remove NaN's from profile, add to 'profile_all' lists
%             u_raw_profile = u_raw_profile(ind_notnan);
%             u_cal_profile = u_cal_profile(ind_notnan);
%             z_profile = z_profile(ind_notnan);
%             
%             %add to lists
%             u_raw_profile_all{i}{j} = u_raw_profile;
%             u_cal_profile_all{i}{j} = u_cal_profile;
%             z_profile_all{i}{j} = z_profile;
%             Anemometer_profile_all{i}{j} = Anemometer_profile{ind_notnan};
%             
%             %get profile for three lowest anemometers
%             [~, ind_sort] = sort(z_profile);
%         
%             %check to make sure there are enough values for fitting log
%             if length(ind_sort)>=3
%                 ind_lowest = ind_sort(1:3);
%                 u_raw_lowest = u_raw_profile(ind_lowest);
%                 u_cal_lowest = u_cal_profile(ind_lowest);
%                 z_lowest = z_profile(ind_lowest);
% 
%                 %compute u* and z0, add to list
%                 P_raw = polyfit(log(z_lowest),u_raw_lowest,1);
%                 ustLog_raw = kappa*P_raw(1);
%                 z0_raw = exp(-P_raw(2)/P_raw(1));
%                 P_cal = polyfit(log(z_lowest),u_cal_lowest,1);
%                 ustLog_cal = kappa*P_cal(1);
%                 z0_cal = exp(-P_cal(2)/P_cal(1));
% 
%                 %deal with poorly fit profiles
%                 if(ustLog_raw<0||isnan(ustLog_raw))
%                     ustLog_raw = NaN;
%                     z0_raw = NaN;
%                 end
%                 if(ustLog_cal<0||isnan(ustLog_cal))
%                     ustLog_cal = NaN;
%                     z0_cal = NaN;
%                 end
%             
%                 %add to lists
%                 u_raw_lowest_all{i}{j} = u_raw_lowest;
%                 u_cal_lowest_all{i}{j} = u_cal_lowest;
%                 z_lowest_all{i}{j} = z_lowest;
%                 Anemometer_lowest_all{i}{j} = Anemometer_profile{ind_lowest};
%                 ustLog_raw_all{i}(j) = ustLog_raw;
%                 z0_raw_all{i}(j) = z0_raw;
%                 ustLog_cal_all{i}(j) = ustLog_cal;
%                 z0_cal_all{i}(j) = z0_cal;
%             end
%         end 
    end
    
    %KEEP ONLY INTERVALS WHERE FLUX AND STRESS ARE WELL DEFINED
    ind_good = intersect(find(~isnan(Q_all{i})),find(~isnan(ustRe_raw_all{i})));

    %flux values
    Q_all{i} = Q_all{i}(ind_good); %total flux
    zq_all{i} = zq_all{i}(ind_good); %characteristic flux height
    fD_all{i} = fD_all{i}(ind_good,:); %Wenglor detection frequency
    fD1_all{i} = fD1_all{i}(ind_good); %Wenglor 1 s detection frequency
    fQ_all{i} = fQ_all{i}(ind_good,:); %Wenglor total transport frequency
    fQ1_all{i} = fQ1_all{i}(ind_good); %Wenglor 1 s transport frequency
    fQalt_all{i} = fQalt_all{i}(ind_good); %Wenglor total transport frequency, normalized by Poisson expectation
    Q1binary_all{i} = Q1binary_all{i}(ind_good); %timesteps with flux for 1 second wind average
    q_all{i} = q_all{i}(ind_good); %partial flux
    zW_all{i} = zW_all{i}(ind_good); %Wenglor flux height
    qcal_all{i} = qcal_all{i}(ind_good); %calibration factor
    
    %wind values
    u_bar_raw_all{i} = u_bar_raw_all{i}(ind_good); %raw mean u values
    u_std_raw_all{i} = u_std_raw_all{i}(ind_good); %standard deviation u values
    u2_raw_all{i} = u2_raw_all{i}(ind_good); %mean u^2
    u_bar_cal_all{i} = u_bar_cal_all{i}(ind_good); %raw mean u values
    u_std_cal_all{i} = u_std_cal_all{i}(ind_good); %standard deviation u values
    u2_cal_all{i} = u2_cal_all{i}(ind_good); %mean u^2
    uth_TFEM_all{i} = uth_TFEM_all{i}(ind_good); %uth for 1s TFEM
    
    %stress values
    ustRe_raw_all{i} = ustRe_raw_all{i}(ind_good); %raw u* for Reynolds stress
    tauRe_raw_all{i} = tauRe_raw_all{i}(ind_good); %raw tau for Reynolds stress
    ustRe_cal_all{i} = ustRe_cal_all{i}(ind_good); %calibrated u* for Reynolds stress
    tauRe_cal_all{i} = tauRe_cal_all{i}(ind_good); %calibrated tau for Reynolds stress
    
    %stress partition values
    eta_fQ_Q0_all{i} = eta_fQ_Q0_all{i}(ind_good); %mean stress partition based on flux frequency method, times of zero flux
    eta_fQ_tauit_all{i} = eta_fQ_tauit_all{i}(ind_good); %mean stress partition based on flux frequency method, constant tau_it
    eta_fQ_TFEM_all{i} = eta_fQ_TFEM_all{i}(ind_good); %mean stress partition based on flux frequency method, TFEM u_th
    eta_zs_all{i} = eta_zs_all{i}(ind_good); %mean stress partition based on friction coefficient estimation from z_s, constant tau_it
    
    %grain size values
    d50_all{i} = d50_all{i}(ind_good); %d50 of surface sand
    d10_all{i} = d10_all{i}(ind_good); %d10 of surface sand
    d90_all{i} = d90_all{i}(ind_good); %d90 of surface sand
    
    %date/time values
    date_all{i} = date_all{i}(ind_good); %date of observation
    StartTime_all{i} = BlockStartTimes(ind_good); %start time of block
    EndTime_all{i} = BlockEndTimes(ind_good); %end time of block
    
%     %profile values
%     u_raw_profile_all{i} = u_raw_profile_all{i}(ind_good); %raw velocity profile
%     u_cal_profile_all{i} = u_cal_profile_all{i}(ind_good); %calibrated velocity profile
%     z_profile_all{i} = z_profile_all{i}(ind_good); %heights for profile
%     Anemometer_profile_all{i} = Anemometer_profile_all{i}(ind_good); %names of anemometers for profile
%     u_raw_lowest_all{i} = u_raw_lowest_all{i}(ind_good); %raw velocity profile - lowest anemometers
%     u_cal_lowest_all{i} = u_cal_lowest_all{i}(ind_good); %calibrated velocity profile - lowest anemometers
%     z_lowest_all{i} = z_lowest_all{i}(ind_good); %heights of lowest anemometers
%     Anemometer_lowest_all{i} = Anemometer_lowest_all{i}(ind_good); %names of lowest anemometers for profile
%     ustLog_raw_all{i} = ustLog_raw_all{i}(ind_good); %u* calculated from mean of raw velocities for lowest anemometers in profile
%     z0_raw_all{i} = z0_raw_all{i}(ind_good); %z0 calculated from mean of raw velocities for lowest anemometers in profile
%     ustLog_cal_all{i} = ustLog_cal_all{i}(ind_good); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
%     z0_cal_all{i} = z0_cal_all{i}(ind_good); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile
end

% SAVE DATA
save(SaveData_Path,'Sites','*all');