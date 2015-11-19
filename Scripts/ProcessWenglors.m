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
        ProcessedWenglors.(WenglorNames{i})(j).qz = struct('qz',zeros(N_t,1),'units','g/m^2/s','qzPerCount',zeros(N_t,1));
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

            %extract information about Wenglor and reference distance sensor for primary interval
            HeightRef = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).HeightRef; %get reference distance sensor for Wenglor height calculation
            z_ref_mm = ExtractVariableTimeInterval(Data.Distance.(HeightRef),StartTimeBSNE_Primary,EndTimeBSNE_Primary,'z','int','int'); %get reference elevations for this interval (mm)
            z_rel = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).StartHeight_m; %get relative height of Wenglor (m)
            sigma_z_rel = ProcessedWenglors.(WenglorNames{j})(PrimaryInterval_IntervalNumber).HeightErr_m; %get uncertainty in relative height of Wenglor from reported value (m)
            
            %get mean and uncertainty in reference height
            z_ref = mean(z_ref_mm*1e-3); %get mean reference height (m)
            sigma_z_ref = std(z_ref_mm*1e-3); %get uncertainty in reference height based on standard deviation (m)
            
            %get mean and uncertainty in overall Wenglor height
            z_Wenglor = z_ref+z_rel; %get Wenglor height by combining reference and relative values
            sigma_z_Wenglor = sqrt(sigma_z_ref.^2+sigma_z_rel.^2); %get total error in Wenglor height by adding relative and reference errors in quadrature
            
            %compute expected flux at Wenglor height
            qz_pred = q0_BSNE.*exp(-z_Wenglor/zq_BSNE); %(g/m^2/s)
            
            %convert uncertainty in Wenglor height to uncertainty in flux
            sigma_qz_z = (-q0_BSNE/zq_BSNE)*exp(-z_Wenglor/zq_BSNE)*sigma_z_Wenglor;
            
            %get uncertainty in fitted flux based on uncertainty of fit for nearest BSNE in profile
            ind_nearestBSNE = find(abs(FluxBSNE(i).z.z-z_Wenglor)==min(abs(FluxBSNE(i).z.z-z_Wenglor)),1);
            sigma_qz_fit = FluxBSNE(i).qz.sigma_qz_pred(ind_nearestBSNE);
            
            %add uncertainties due to z and fit to get total uncertainty in predicted qz
            sigma_qz_pred = sqrt(sigma_qz_z^2+sigma_qz_fit^2);
            
            %compute counts per second for Wenglor during BSNE time interval
            T_BSNE = seconds(EndTimeBSNE_Primary - StartTimeBSNE_Primary);
            W_CountsPerSecond = sum(PrimaryInterval_n)/T_BSNE;
            
            %get time increment for Wenglor observations
            W_dt = seconds(mode(diff(PrimaryInterval_t)));

            %compute conversion factor from Wenglor counts to flux, "qzPerCount"
            qzPerCount = qz_pred/(W_CountsPerSecond*W_dt);
            
            %uncertainty estimation for calibration coefficient
            sigma_qzPerCount = (1/W_CountsPerSecond)*sqrt(sigma_qz_pred^2+qz_pred^2/(T_BSNE*W_CountsPerSecond));
            
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
    'z',[],'sigma_z',[],'WenglorNames',[],'Units','g/m^2/s'),...
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
end