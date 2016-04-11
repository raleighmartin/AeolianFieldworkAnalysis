%% function to process Wenglor calibration
 
function [ProcessedWenglors, FluxWenglor] = ProcessWenglors(Data, FluxBSNE, InstrumentMetadata)

%% GO THROUGH BSNE PROFILES TO CALIBRATE WENGLOR COUNTS TO HEIGHT-SPECIFIC FLUXES
%initialize ProcessedWenglors from raw (interpolated) Wenglor data
ProcessedWenglors = Data.Wenglor;

%get names of Wenglors
WenglorNames = fieldnames(ProcessedWenglors);
N_Wenglors = length(WenglorNames);

%go through each Wenglor
for i=1:N_Wenglors
    
    %compute number of intervals for Wenglor
    N_Intervals = length(ProcessedWenglors.(WenglorNames{i}));
    
    %go through each interval, initialize flux column
    for j=1:N_Intervals
        
        %get number of time steps for flux calc
        N_t = length(ProcessedWenglors.(WenglorNames{i})(j).t.int);
        
        %assign 'qz' structured array to new column of master stuctured array
        units = struct('qz','g/m^2/s','z','m','qzPerCount','g/m^2'); %create structured array with units
        ProcessedWenglors.(WenglorNames{i})(j).qz = ...
            struct('qz',zeros(N_t,1),'z',zeros(N_t,1),...
            'sigma_z',zeros(N_t,1),'qzPerCount',zeros(N_t,1),...
            'sigma_qzPerCount',zeros(N_t,1),'units',units);
    end
    
end

%% go through BSNE profiles for calibration
for i=1:length(FluxBSNE);

    %get start and end time of BSNE profile
    StartTimeBSNE_Primary = FluxBSNE(i).StartTime %print on command window to track
    EndTimeBSNE_Primary = FluxBSNE(i).EndTime;

    %get BSNE flux information for time interval
    q0_BSNE = FluxBSNE(i).qz.q0; %qz at 0, (g/m^2/s)
    zq_BSNE = FluxBSNE(i).z.zq; %get BSNE e-folding height (m)
    sigma_q0_BSNE = FluxBSNE(i).qz.sigma_q0;
    sigma_zq_BSNE = FluxBSNE(i).z.sigma_zq;
    sigma2_q0zq_BSNE = FluxBSNE(i).Q.sigma2_q0zq;
    
    %create enlarged time intervals for calibrating Wenglor times outside of BSNE profile times
    if i == 1 %earliest start time is beginning of first BSNE day
        StartTimeBSNE_Enlarged = FluxBSNE(i).Date;
        EndTimeBSNE_Enlarged = mean([FluxBSNE(i).EndTime, FluxBSNE(i+1).StartTime]);
    elseif i == length(FluxBSNE) %latest start time is end of last BSNE day
        StartTimeBSNE_Enlarged = mean([FluxBSNE(i-1).EndTime, FluxBSNE(i).StartTime]);
        EndTimeBSNE_Enlarged = FluxBSNE(i).Date + 1;
    else %middle start and end times are spaced halfway between BSNE end time before and start time after
        StartTimeBSNE_Enlarged = mean([FluxBSNE(i-1).EndTime, FluxBSNE(i).StartTime]);
        EndTimeBSNE_Enlarged = mean([FluxBSNE(i).EndTime, FluxBSNE(i+1).StartTime]);
    end

    %go through each Wenglor for calibration
    for j = 1:N_Wenglors
        
        %get Wenglor counts and interval number (there should be only one) for primary interval
        [PrimaryInterval_n, PrimaryInterval_t, PrimaryInterval_IntervalNumber, ~] = ...
            ExtractVariableTimeInterval(ProcessedWenglors.(WenglorNames{j}),StartTimeBSNE_Primary,EndTimeBSNE_Primary,'n','int','int');
        
        %perform calculations assuming single set of Wenglor counts within interval, otherwise ignore
        if length(PrimaryInterval_IntervalNumber)==1
            
            %get Wenglor height and uncertainty
            z_Wenglor = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).z.z;
            sigma_z_Wenglor = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).z.sigma_z;
            
%             %compute expected flux at Wenglor height
%             qz_pred = q0_BSNE.*exp(-z_Wenglor/zq_BSNE); %(g/m^2/s)
%             
%             %compute uncertainty in expected flux at Wenglor height
%             sigma_qz_pred = sqrt((sigma_q0_BSNE*qz_pred/q0_BSNE).^2 + ...
%                 (sigma_zq_BSNE*qz_pred.*z_Wenglor/zq_BSNE^2).^2 + ...
%                 (2*sigma2_q0zq_BSNE*qz_pred.^2.*z_Wenglor/(q0_BSNE*zq_BSNE^2)));
            
            %compute expected flux and uncertainty at Wenglor height
            [qz_pred, sigma_qz_pred] = qz_prediction(z_Wenglor, q0_BSNE, zq_BSNE, sigma_q0_BSNE, sigma_zq_BSNE, sigma2_q0zq_BSNE);
            
            %compute uncertainty due height uncertainty of Wenglor
            sigma_qz_z = abs((qz_pred/zq_BSNE)*sigma_z_Wenglor); %(g/m^2/s)
            
            %compute calibration flux and uncertainty
            qz_cal = qz_pred; %calibration flux equals BSNE expected flux
            sigma_qz_cal = sqrt(sigma_qz_pred.^2+sigma_qz_z.^2); %calibration flux uncertainty equals expected flux uncertainty plus Wenglor height uncertainty
            
%             %get contribution of uncertainty in Wenglor height to uncertainty in flux
%             sigma_qz_z = abs((-q0_BSNE/zq_BSNE)*exp(-z_Wenglor/zq_BSNE)*sigma_z_Wenglor);
%             
%             %get contribution of uncertainty in BSNE q0 to uncertainty in flux
%             sigma_qz_q0 = abs(exp(-z_Wenglor/zq_BSNE)*sigma_q0_BSNE);
%             sigma_qz_zq = abs((-q0_BSNE*z_Wenglor)/(zq_BSNE^2)*exp(-z_Wenglor/zq_BSNE)*sigma_zq_BSNE);
%             
%             %add uncertainties to get total uncertainty in predicted qz
%             sigma_qz_pred = sqrt(sigma_qz_z^2+sigma_qz_q0^2+sigma_qz_zq^2);
%             
% %             %get uncertainty in fitted flux based on relative uncertainty of fit for nearest BSNE in profile
% %             ind_nearestBSNE = find(abs(FluxBSNE(i).z.z-z_Wenglor)==min(abs(FluxBSNE(i).z.z-z_Wenglor)),1);
% %             sigma_qz_fit = FluxBSNE(i).qz.sigma_qz_pred(ind_nearestBSNE);
% %             
% %             %add uncertainties due to z and fit to get total uncertainty in predicted qz
% %             sigma_qz_pred = sqrt(sigma_qz_z^2+sigma_qz_fit^2);
            
            %compute counts per second for Wenglor during BSNE time interval
            T_BSNE = seconds(EndTimeBSNE_Primary - StartTimeBSNE_Primary);
            W_CountsPerSecond = sum(PrimaryInterval_n)/T_BSNE;
%             sigma_W_CountsPerSecond = sqrt(sum(PrimaryInterval_n))/T_BSNE; %counting error
            
            %get time increment for Wenglor observations
            W_dt = seconds(mode(diff(PrimaryInterval_t)));

            %compute conversion factor from Wenglor counts to flux, "qzPerCount"
            qzPerCount = qz_cal/(W_CountsPerSecond*W_dt);
            sigma_qzPerCount = sigma_qz_cal/(W_CountsPerSecond*W_dt); %calibration uncertainty directly from uncertainty in qz pred
            
%             %uncertainty estimation for calibration coefficient
%             sigma_qzPerCount_qz_pred = abs(sigma_qz_pred/(W_CountsPerSecond*W_dt)); %uncertainty due to qz_pred
%             sigma_qzPerCount_qzPerCount = abs(sigma_W_CountsPerSecond*(qz_pred/(W_CountsPerSecond.^2*W_dt))); %uncertainty due to W_CountsPerSecond
%             sigma_qzPerCount = sqrt(sigma_qzPerCount_qz_pred.^2+sigma_qzPerCount_qzPerCount.^2);
            
%             %uncertainty estimation for calibration coefficient
%             sigma_qzPerCount = (1/W_CountsPerSecond)*sqrt(sigma_qz_pred^2+qz_pred^2/(T_BSNE*W_CountsPerSecond));
            
            %get interval numbers for enlarged interval subject to this calibration (there may be more than one)
            [~, ~, EnlargedInterval_IntervalNumbers, EnlargedInterval_IntervalIndices] = ...
                ExtractVariableTimeInterval(ProcessedWenglors.(WenglorNames{j}),StartTimeBSNE_Enlarged,EndTimeBSNE_Enlarged,'n','int','int');

            %go through each of these intervals and perform calibration
            N_IntervalNumbers = length(EnlargedInterval_IntervalNumbers);

            for k = 1:N_IntervalNumbers

                %get list of counts for calibration
                n_list = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).n.int(EnlargedInterval_IntervalIndices{k});

                %convert to flux values
                qz_list = qzPerCount*n_list;
                                
                %Assign fluxes, calibration factors, and heights to structured array
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.qz(EnlargedInterval_IntervalIndices{k}) = qz_list; %height-specific particle fluxes
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.qzPerCount(EnlargedInterval_IntervalIndices{k}) = qzPerCount; %calibration factor
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.sigma_qzPerCount(EnlargedInterval_IntervalIndices{k}) = sigma_qzPerCount; %uncertainty in calibration factor
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.z(EnlargedInterval_IntervalIndices{k}) = z_Wenglor; %Wenglor height
                ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.sigma_z(EnlargedInterval_IntervalIndices{k}) = sigma_z_Wenglor; %uncertainty in Wenglor height
            end
            
            %ensure that values are in column vectors
            qz = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.qz;
            ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.qz = reshape(qz,[length(qz),1]);
            qzPerCount = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.qzPerCount;
            ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.qzPerCount = reshape(qzPerCount,[length(qzPerCount),1]);
            sigma_qzPerCount = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.sigma_qzPerCount;
            ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.sigma_qzPerCount = reshape(sigma_qzPerCount,[length(sigma_qzPerCount),1]);
            z = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.z;
            ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.z = reshape(z,[length(z),1]);
            sigma_z = ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.sigma_z;
            ProcessedWenglors.(WenglorNames{j})(EnlargedInterval_IntervalNumbers(k)).qz.sigma_z = reshape(sigma_z,[length(sigma_z),1]);
        end
    end
end

%% Create 'FluxWenglor' structured array from calibrated Wenglors
FluxWenglor_dt = duration(0,0,0.04); %time step for Wenglor flux matrices

%get indices of valid Wenglor times in metadata table
ind_MetadataWenglors = intersect(find(strcmp(InstrumentMetadata.InstrumentType,'Wenglor')),find(InstrumentMetadata.ErrorCode==0));

%get raw start and end times for total flux calcs (i.e., start and end times of all possible intervals without any sorting)
FluxWenglor_RawStartTimes = InstrumentMetadata.StartTime(ind_MetadataWenglors);
FluxWenglor_RawEndTimes = InstrumentMetadata.EndTime(ind_MetadataWenglors);
FluxWenglor_RawTimes = unique(union(FluxWenglor_RawStartTimes,FluxWenglor_RawEndTimes));

%compute start times of overlapping intervals
[FluxWenglor_OverlappingStartTimes, FluxWenglor_OverlappingEndTimes] = CombineIntervals(FluxWenglor_RawStartTimes, FluxWenglor_RawEndTimes);
N_OverlappingIntervals = length(FluxWenglor_OverlappingStartTimes);

%now, subdivide these into all distinctive intervals for computing total flux
FluxWenglor_StartTimes = []; %initialize list of all distinctive start times
FluxWenglor_EndTimes = []; %initialize list of all distinctive end times
for i = 1:N_OverlappingIntervals
    OverlappingInterval_AllStartTimes = unique(FluxWenglor_RawTimes(...
        FluxWenglor_RawTimes>=FluxWenglor_OverlappingStartTimes(i)&...
        FluxWenglor_RawTimes<FluxWenglor_OverlappingEndTimes(i)));
    FluxWenglor_StartTimes = [FluxWenglor_StartTimes; OverlappingInterval_AllStartTimes];
    OverlappingInterval_AllEndTimes = unique(FluxWenglor_RawTimes(...
        FluxWenglor_RawTimes>FluxWenglor_OverlappingStartTimes(i)&...
        FluxWenglor_RawTimes<=FluxWenglor_OverlappingEndTimes(i)));
    FluxWenglor_EndTimes = [FluxWenglor_EndTimes; OverlappingInterval_AllEndTimes];
end

%extract associated dates for entries
[y,m,d] = ymd(FluxWenglor_StartTimes);
FluxWenglor_Dates = datetime(y,m,d);

%extract associated sites for entries
N_FluxWenglorIntervals = length(FluxWenglor_StartTimes);
FluxWenglor_Sites = cell(N_FluxWenglorIntervals,1);
for i=1:N_FluxWenglorIntervals
    FluxWenglor_Sites{i} = char(unique(InstrumentMetadata.Site(InstrumentMetadata.Date==FluxWenglor_Dates(i))));
end

%initialize structured array for total flux from profile
FluxWenglor = struct(...
    'Site', cellstr(FluxWenglor_Sites),...
    'Date', num2cell(FluxWenglor_Dates),...
    'StartTime', num2cell(FluxWenglor_StartTimes),...
    'EndTime', num2cell(FluxWenglor_EndTimes),...
    't',struct('t',[],'dt',FluxWenglor_dt),...
    'qz',struct('n',[],'qz',[],'qzPerCount',[],'sigma_qzPerCount',[],...
    'z',[],'sigma_z',[],'WenglorNames',[],'WenglorID',[],'Units','g/m^2/s'),...
    'z',struct('z',[],'Units','m'));

%figure out which Wenglors are contained in each interval
for i = 1:N_Wenglors
    ind_MetadataWenglor = intersect(find(strcmp(InstrumentMetadata.Instrument,WenglorNames{i})),find(InstrumentMetadata.ErrorCode==0));
    Wenglor_StartTimes = InstrumentMetadata.StartTime(ind_MetadataWenglor);
    Wenglor_EndTimes = InstrumentMetadata.EndTime(ind_MetadataWenglor);
    for j = 1:N_FluxWenglorIntervals
        if ~isempty(intersect(find(Wenglor_StartTimes<=FluxWenglor(j).StartTime),find(Wenglor_EndTimes>=FluxWenglor(j).EndTime)))
            FluxWenglor(j).qz.WenglorNames = [FluxWenglor(j).qz.WenglorNames, {WenglorNames{i}}];
        end
    end
end

%fill in times and initialize components of structured array
for i = 1:N_FluxWenglorIntervals
    %print data for script tracking purposes
    FluxWenglor(i).Date
    
    %calculate times
    t = (FluxWenglor(i).StartTime:FluxWenglor_dt:FluxWenglor(i).EndTime)';
    N_t = length(t);
    FluxWenglor(i).t.t = t;
    
    %determine error times, create field in FluxWenglor for this
    WenglorNames = FluxWenglor(i).qz.WenglorNames;
    N_Wenglors = length(WenglorNames);
    t_err_list = []; %initialize list of error times
    for j = 1:N_Wenglors %go through each Wenglor in profile to get error times
        [~, ~, t_IntervalNumber, ~] = ExtractVariableTimeInterval(ProcessedWenglors.(WenglorNames{j}),...
            min(t),max(t),'t','int','int'); %get data window containing error times
        t_err_list = [t_err_list; ProcessedWenglors.(WenglorNames{j})(t_IntervalNumber).t.err]; %get error times and add to list
    end
    t_err_list = intersect(unique(t_err_list),t); %get only unique times in window of times for flux profiles
    FluxWenglor(i).t.err = t_err_list; %add these times to structured array
    
    %get information about Wenglors in profile
    N_qz = length(FluxWenglor(i).qz.WenglorNames);
    FluxWenglor(i).z.z = zeros(1,N_qz);

    %initialize partial fluxes
    FluxWenglor(i).qz.qz = zeros(N_t,N_qz)*NaN;
    FluxWenglor(i).qz.qzPerCount = zeros(N_t,N_qz)*NaN;
    
    %initialize lists of qz_StartInd and qz_EndInd
    qz_StartInd = zeros(N_qz,1);
    qz_EndInd = zeros(N_qz,1);
    
    %initialize lists of Wenglor IDs
    FluxWenglor(i).qz.WenglorID = cell(1,N_qz);
    
    %fill in values
    for j = 1:N_qz
        WenglorName = FluxWenglor(i).qz.WenglorNames{j}; %get specific Wenglor name
        
        %extract information abouve fluxes in time interval
        [~, qz_t, qz_IntervalNumber, qz_ind] = ...
            ExtractVariableTimeInterval(ProcessedWenglors.(WenglorName),FluxWenglor(i).StartTime,FluxWenglor(i).EndTime,'qz','qz','int');
        qz_ind = cell2mat(qz_ind); %make qz_ind into a vector
        
        %line things up with times in flux table
        qz_StartInd(j) = find(t==qz_t(1)); %index in table of first qz from Wenglor, keep track of it
        qz_EndInd(j) = find(t==qz_t(end)); %index in table of last qz from Wenglor, keep track of it
        
        %add values to FluxWenglor array
        FluxWenglor(i).qz.n(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).n.int(qz_ind); %submit counts values into flux table       
        FluxWenglor(i).qz.qz(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).qz.qz(qz_ind); %submit qz's into flux table
        FluxWenglor(i).qz.qzPerCount(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).qz.qzPerCount(qz_ind); %submit calibration values into flux table
        FluxWenglor(i).qz.sigma_qzPerCount(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).qz.sigma_qzPerCount(qz_ind); %submit uncertainty on calibration values into flux table
        FluxWenglor(i).qz.z(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).qz.z(qz_ind); %submit height values into flux table
        FluxWenglor(i).qz.sigma_z(qz_StartInd(j):qz_EndInd(j),j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).qz.sigma_z(qz_ind); %submit uncertainty on height values into flux table
                
        %get information about average Wenglor height
        FluxWenglor(i).z.z(j) = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).z.z;
        
        %get information about Wenglor ID
        FluxWenglor(i).qz.WenglorID{j} = ProcessedWenglors.(WenglorName)(qz_IntervalNumber).InstrumentID;
    end
    
    %sort values by z, compute order of data columns from this
    [z_sorted, z_sort_ind] = sort(FluxWenglor(i).z.z);

    %get indices of rows that are contained in data to keep
    row_keep_ind = max(qz_StartInd):min(qz_EndInd);
    
    %modify matrices based on sorting and determination of rows to keep
    FluxWenglor(i).z.z = z_sorted;
    FluxWenglor(i).t.t = FluxWenglor(i).t.t(row_keep_ind);
    FluxWenglor(i).qz.n = FluxWenglor(i).qz.n(row_keep_ind, z_sort_ind);
    FluxWenglor(i).qz.qz = FluxWenglor(i).qz.qz(row_keep_ind, z_sort_ind);
    FluxWenglor(i).qz.qzPerCount = FluxWenglor(i).qz.qzPerCount(row_keep_ind, z_sort_ind);
    FluxWenglor(i).qz.sigma_qzPerCount = FluxWenglor(i).qz.sigma_qzPerCount(row_keep_ind, z_sort_ind);
    FluxWenglor(i).qz.z = FluxWenglor(i).qz.z(row_keep_ind, z_sort_ind);
    FluxWenglor(i).qz.sigma_z = FluxWenglor(i).qz.sigma_z(row_keep_ind, z_sort_ind);
    FluxWenglor(i).qz.WenglorNames = FluxWenglor(i).qz.WenglorNames(z_sort_ind);
    FluxWenglor(i).qz.WenglorID = FluxWenglor(i).qz.WenglorID(z_sort_ind);
end