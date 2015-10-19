%% SCRIPT TO GENERATE WINDOWS OF SALTATION FLUX AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;

%% parameter values
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

%% information about sites for analysis
%Sites = {'Jericoacoara'};
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
N_RunningPerProfile = floor(ProfileTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times

%% set duration of window average for flux frequencies
T_windowaverage = duration(0,0,1);

%% information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
SaveData_Path = strcat(folder_ProcessedData,'StressFluxWindows');

%% load processed and metadata for each site, add to structured arrays of all data and metadata
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Sites{i});
    Data{i} = load(ProcessedData_Path); %load processed data
    Metadata{i} = load(Metadata_Path); %load metadata
end

%% initialize variable lists

%initialize lists of flux values
Q_list = cell(N_Sites,1); %total flux
ze_list = cell(N_Sites,1); %flux height
fQ_list = cell(N_Sites,1); %Wenglor total transport frequency list
q_list = cell(N_Sites,1); %partial flux
zq_list = cell(N_Sites,1); %Wenglor flux height list
fq_list = cell(N_Sites,1); %Wenglor transport frequency list
qcal_list = cell(N_Sites,1); %calibration factor list

%initialize lists of wind values
ustRe_raw_list = cell(N_Sites,1); %raw u* for Reynolds stress
tauRe_raw_list = cell(N_Sites,1); %raw tau for Reynolds stress
ustRe_cal_list = cell(N_Sites,1); %calibrated u* for Reynolds stress
tauRe_cal_list = cell(N_Sites,1); %calibrated tau for Reynolds stress

%initialize lists of grain sizes, dates, start times, and end times
d50_list = cell(N_Sites,1); %lists of surface grain sizes corresponding to calculated values
d10_list = cell(N_Sites,1);
d90_list = cell(N_Sites,1);
date_list = cell(N_Sites,1); %lists of dates corresponding to calculated values
StartTime_list = cell(N_Sites,1); %lists of block start times
EndTime_list = cell(N_Sites,1); %lists of block end times

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
    
    %% choose anemometer type based on site of interest
    if strcmp(Sites{i},'Oceano')
        AnemometerType = 'Sonic';
        Anemometer = 'S1';
    elseif strcmp(Sites{i},'RanchoGuadalupe')
        AnemometerType = 'Ultrasonic';
        Anemometer = 'U1';
    elseif strcmp(Sites{i},'Jericoacoara')
        AnemometerType = 'Ultrasonic';
        Anemometer = 'U1';
    end
    
    %% extract wind and flux data from overall processed data file
    WindData = Data{i}.ProcessedData.(AnemometerType).(Anemometer);
    FluxData = Data{i}.ProcessedData.FluxWenglor;
    
    %% get start times and end times for flux and wind observations
    WindStartTimes = [WindData.StartTime]';
    WindEndTimes = [WindData.EndTime]';
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
    Q_list{i} = zeros(N_Blocks,1)*NaN; %total flux
    ze_list{i} = zeros(N_Blocks,1)*NaN; %flux height
    fQ_list{i} = zeros(N_Blocks,1)*NaN; %Wenglor total transport frequency
    q_list{i} = cell(N_Blocks,1); %partial flux
    zq_list{i} = cell(N_Blocks,1); %Wenglor flux height
    fq_list{i} = cell(N_Blocks,1); %Wenglor transport frequency
    qcal_list{i} = cell(N_Blocks,1); %calibration factor
    
    %initialize lists of wind values
    ustRe_raw_list{i} = zeros(N_Blocks,1)*NaN; %raw u* for Reynolds stress
    tauRe_raw_list{i} = zeros(N_Blocks,1)*NaN; %raw tau for Reynolds stress
    ustRe_cal_list{i} = zeros(N_Blocks,1)*NaN; %calibrated u* for Reynolds stress
    tauRe_cal_list{i} = zeros(N_Blocks,1)*NaN; %calibrated tau for Reynolds stress
    
    %create list of dates
    date_list{i} = datetime(BlockStartTimes.Year, BlockStartTimes.Month, BlockStartTimes.Day); %lists of dates corresponding to calculated values
    
    %initialize lists of grain sizes
    d50_list{i} = zeros(N_Blocks,1); %lists of surface grain sizes corresponding to calculated values
    d10_list{i} = zeros(N_Blocks,1);
    d90_list{i} = zeros(N_Blocks,1);
    
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
            
            %compute mean qz profile
            q_profile = mean(FluxData(IntervalN).qz.qz(IntervalInd,:));
            zq_profile = FluxData(IntervalN).z.z;

            %Perform profile fit to get q0, ze, and Q
            [q0,ze] = qz_profilefit(q_profile,zq_profile,0,0); %For now, ignore uncertainty in profile values
            Q = q0*ze; %get total flux [g/m/s]
            
            %create 1 second window-averaged flux timeseries
            N_zq = length(zq_profile);
            qz_avg = cell(1,N_zq);
            for k = 1:N_zq
                [qz_avg_z, ~] = window_average(FluxData(IntervalN).qz.qz(IntervalInd,k), t_Interval, T_windowaverage);
                qz_avg{k} = qz_avg_z;
            end
            qz_avg = cell2mat(qz_avg);
            
            %get frequency of detection
            fq_profile = sum(qz_avg>0)/length(qz_avg);
            fQ = (sum(sum(qz_avg')'>0))/length(qz_avg);
                  
            %get calibration values
            qcal_profile = mean(FluxData(IntervalN).qz.qzPerCount(IntervalInd,:));
            
            %convert to 0 if NaN
            if isnan(Q)
                Q=0;
            end
            if isnan(ze)
                ze=0;
            end
            
            %add to list
            Q_list{i}(j) = Q; %total flux
            ze_list{i}(j) = ze; %flux height
            fQ_list{i}(j) = fQ; %frquency of transport
            q_list{i}{j} = q_profile; %partial flux
            zq_list{i}{j} = zq_profile; %Wenglor flux height
            fq_list{i}{j} = fq_profile; %Wenglor transport frequency
            qcal_list{i}{j} = qcal_profile; %calibration factor
        end
        
        %% WIND CALCULATIONS FOR INTERVAL
        %extract time interval
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData,StartTime,EndTime,'u','int','int');
        
        %determine if values span entire analysis interval
        if ~isempty(IntervalN)
            %use only first interval
            IntervalN = IntervalN(1);
            IntervalInd = IntervalInd{1};
            StartTime_Extraction = WindData(IntervalN).t.int(IntervalInd(1));
            EndTime_Extraction = WindData(IntervalN).t.int(IntervalInd(end));
            
            %further reduce IntervalInd based on eliminating error times
            [~, ErrInd, ~] = intersect(WindData(IntervalN).t.int,WindData(IntervalN).t.err);
            IntervalInd = setdiff(IntervalInd,ErrInd);
            
        else %generate arbitrary values if there are no data
            StartTime_Extraction = datetime(0,0,0);
            EndTime_Extraction = datetime(0,0,0);
        end
        
        %restrict to time intervals with data spanning entire interval
        if (StartTime_Extraction==StartTime)&&(EndTime_Extraction==EndTime)
            
            %get velocity values
            u = WindData(IntervalN).u.int(IntervalInd);
            v = WindData(IntervalN).v.int(IntervalInd);
            w = WindData(IntervalN).w.int(IntervalInd);
            
            %rotate instrument, call these 'raw' values
            [u_raw, ~, w_raw] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

            %get calibration factors
            CalibrationFactor_u = WindData(IntervalN).u.CalibrationFactor;
            CalibrationFactor_w = WindData(IntervalN).w.CalibrationFactor;
            CalibrationIntercept_u = WindData(IntervalN).u.CalibrationIntercept;
            CalibrationIntercept_w = WindData(IntervalN).w.CalibrationIntercept;

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

            %add to list
            ustRe_raw_list{i}(j) = ustRe_raw;
            tauRe_raw_list{i}(j) = rho_a*ustRe_raw.^2;
            ustRe_cal_list{i}(j) = ustRe_cal;
            tauRe_cal_list{i}(j) = rho_a*ustRe_cal.^2;
        end
        
        %GRAIN SIZE INFO
        ind_date = find([Metadata{i}.GrainSize_Surface.Date]==WindData(IntervalN).Date);
        d50_list{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_50_mm]);
        d10_list{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_10_mm]);
        d90_list{i}(j) = median([Metadata{i}.GrainSize_Surface(ind_date).d_90_mm]);
    end

    %KEEP ONLY INTERVALS WHERE FLUX AND STRESS ARE WELL DEFINED
    ind_good = intersect(find(~isnan(Q_list{i})),find(~isnan(ustRe_raw_list{i})));

    Q_list{i} = Q_list{i}(ind_good); %total flux
    ze_list{i} = ze_list{i}(ind_good); %flux height
    fQ_list{i} = fQ_list{i}(ind_good); %Wenglor total transport frequency
    q_list{i} = {q_list{i}{ind_good}}'; %partial flux
    zq_list{i} = {zq_list{i}{ind_good}}'; %Wenglor flux height
    fq_list{i} = {fq_list{i}{ind_good}}'; %Wenglor transport frequency
    qcal_list{i} = {qcal_list{i}{ind_good}}'; %calibration factor
    ustRe_raw_list{i} = ustRe_raw_list{i}(ind_good); %raw u* for Reynolds stress
    tauRe_raw_list{i} = tauRe_raw_list{i}(ind_good); %raw tau for Reynolds stress
    ustRe_cal_list{i} = ustRe_cal_list{i}(ind_good); %calibrated u* for Reynolds stress
    tauRe_cal_list{i} = tauRe_cal_list{i}(ind_good); %calibrated tau for Reynolds stress
    d50_list{i} = d50_list{i}(ind_good); %d50 of surface sand
    d10_list{i} = d10_list{i}(ind_good);
    d90_list{i} = d90_list{i}(ind_good);
    date_list{i} = date_list{i}(ind_good); %date of observation
    StartTime_list{i} = BlockStartTimes(ind_good); %start time of block
    EndTime_list{i} = BlockEndTimes(ind_good); %end time of block
end

%% SAVE DATA
save(SaveData_Path,'Sites',...
    'Q_list','ze_list','fQ_list','q_list','zq_list','fq_list','qcal_list',...
    'ustRe_raw_list','tauRe_raw_list','ustRe_cal_list','tauRe_cal_list',...
    'd50_list','d10_list','d90_list','date_list','StartTime_list','EndTime_list');