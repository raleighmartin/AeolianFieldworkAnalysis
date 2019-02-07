%% SCRIPT TO EXPORT SALTATION AND WIND DATA INTO ASCII FILES FOR ANALYSIS OUTSIDE MATLAB
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, ExtractVariableTimeInterval.m

%% initialize
clearvars;

%% information about where to load/save data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data

%info about which data to use
% use_raw_wind = 1; %use raw (non-interpolated) wind data? 0 for no, 1 for yes
use_calibrated_flux = 0; %use calibrated partial (height-specific) flux values? 0 for no (raw Wenglor counts), 1 for yes
use_time_outside_BSNE = 0; %use times outside BSNE intervals (affects calibration)? 0 for no, 1 for yes
zW_min = 0.02; %minimum Wenglor height (m) for inclusion in calculations

% %% information for extracting all 30-minute windows at Jericoacoara
% Site_Export = 'Jericoacoara';
% folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Export/30 Minute/Jericoacoara/'; %folder for storing outputs of this extraction
% load('../../AnalysisData/Windowing/TimeWindows_30min_Restricted.mat'); %load time windows
% ind_Site = find(strcmp(Sites,Site_Export));
% StartTimesExport = StartTime_window{ind_Site}(hasfluxdata_window{ind_Site}==1);
% EndTimesExport = EndTime_window{ind_Site}(hasfluxdata_window{ind_Site}==1);
% zq_Site = 0.097; %estimated saltation height for Site, in meters
% sigma_zq_Site = 0.005; %estimated uncertainty in saltation height for Site, in meters
% WindInstrumentType = 'Ultrasonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
% use_specific_times = 1; %specific times for export set

%% information for extracting all 30-minute windows at Rancho Guadalupe
Site_Export = 'RanchoGuadalupe';
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Export/30 Minute/Rancho Guadalupe/'; %folder for storing outputs of this extraction
load('../../AnalysisData/Windowing/TimeWindows_30min_Restricted.mat'); %load time windows
ind_Site = find(strcmp(Sites,Site_Export));
StartTimesExport = StartTime_window{ind_Site}(hasfluxdata_window{ind_Site}==1);
EndTimesExport = EndTime_window{ind_Site}(hasfluxdata_window{ind_Site}==1);
zq_Site = 0.107; %estimated saltation height for Site, in meters
sigma_zq_Site = 0.005; %estimated uncertainty in saltation height for Site, in meters
WindInstrumentType = 'Ultrasonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
use_specific_times = 1; %specific times for export set

% %% information for extracting all 30-minute windows at Oceano
% Site_Export = 'Oceano';
% folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Export/30 Minute/Oceano/'; %folder for storing outputs of this extraction
% load('../../AnalysisData/Windowing/TimeWindows_30min_Restricted.mat'); %load time windows
% ind_Site = find(strcmp(Sites,Site_Export));
% StartTimesExport = StartTime_window{ind_Site}(hasfluxdata_window{ind_Site}==1);
% EndTimesExport = EndTime_window{ind_Site}(hasfluxdata_window{ind_Site}==1);
% zq_Site = 0.055; %estimated saltation height for Site, in meters
% sigma_zq_Site = 0.004; %estimated uncertainty in saltation height for Site, in meters
% WindInstrumentType = 'Sonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
% use_specific_times = 1; %specific times for export set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% information for full data extraction - Jericoacoara
% Site_Export = 'Jericoacoara'; %name of field site
% Dates = [datetime(2014,11,13);...
%     datetime(2014,11,14);...
%     datetime(2014,11,18);...
%     datetime(2014,11,20)]; %dates for data extraction
% WindInstrumentType = 'Ultrasonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
% 
% %% information for full data extraction - Rancho Guadalupe
% Site_Export = 'RanchoGuadalupe'; %name of field site
% Dates = [datetime(2015,3,23);...
%     datetime(2015,3,24)]; %dates for data extraction
% WindInstrumentType = 'Ultrasonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
%
% %% information for full data extraction - Oceano
% Site_Export = 'Oceano'; %name of field site
% Dates = [ ... 
%     datetime(2015,5,15);...
%     datetime(2015,5,16);...
%     datetime(2015,5,17);...
%     datetime(2015,5,18);...
%     datetime(2015,5,19);...
%     datetime(2015,5,20);...
%     datetime(2015,5,21);...
%     datetime(2015,5,22);...
%     datetime(2015,5,23);...
%     datetime(2015,5,24);...
%     datetime(2015,5,26);...
%     datetime(2015,5,27);...
%     datetime(2015,5,28);...
%     datetime(2015,5,29);...
%     datetime(2015,5,30);...
%     datetime(2015,5,31);...
%     datetime(2015,6,1);...
%     datetime(2015,6,2)...
%     datetime(2015,6,3);...
%     datetime(2015,6,4);...
%     ]; %dates for data extraction
% WindInstrumentType = 'Sonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
% use_specific_times = 0; %no specific times for export set
%
% %% information for limited data extraction
% Site_Export = 'Oceano';
% StartTimesExport = datetime(2015,6,2,14,28,0);
% EndTimesExport = datetime(2015,6,2,15,28,0);
% zq_Site = 0.055; %estimated saltation height for Site, in meters
% sigma_zq_Site = 0.004; %estimated uncertainty in saltation height for Site, in meters
% WindInstrumentType = 'Sonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)
% use_specific_times = 1; %specific times for export set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% load processed data and extract relevant variables
LoadData_Path = strcat(folder_LoadData,'ProcessedData_',Site_Export);
load(LoadData_Path); %load processed data
WenglorData = ProcessedData.FluxWenglor; %Wenglor flux data
BSNEData = ProcessedData.FluxBSNE; %BSNE data
WindData = ProcessedData.(WindInstrumentType); %data for all anemometers
clear ProcessedData; %remove 'ProcessedData' to clear up memory

%% get start times and end times for Wenglor and BSNE observations
WenglorStartTimes = [WenglorData.StartTime]';
WenglorEndTimes = [WenglorData.EndTime]';
BSNEStartTimes = [BSNEData.StartTime]';
BSNEEndTimes = [BSNEData.EndTime]';

%% get start times and end times for wind observations for profile of anemometers
AnemometerName_profile = fieldnames(WindData); %anemometer names
N_Anemometers_profile = length(AnemometerName_profile); %number of anemometers
WindStartTimes_profile = cell(N_Anemometers_profile,1); %start times for wind data in profile
WindEndTimes_profile = cell(N_Anemometers_profile,1); %end times for wind data in profile
for j = 1:N_Anemometers_profile
    WindStartTimes_profile{j} = [WindData.(AnemometerName_profile{j}).StartTime]';
    WindEndTimes_profile{j} = [WindData.(AnemometerName_profile{j}).EndTime]';
end

%% get intersecting time intervals for flux
if use_time_outside_BSNE==0 %in this case, consider time intervals only when BSNE traps are active for calibration
    [FluxStartTimes, FluxEndTimes] = IntersectingTimeIntervals(BSNEStartTimes,BSNEEndTimes,WenglorStartTimes,WenglorEndTimes);
else %otherwise, use all Wenglor data
    FluxStartTimes = WenglorStartTimes;
    FluxEndTimes = WenglorEndTimes;
end

%% get intersecting time intervals for all anemometers in profile
WindStartTimes = WindStartTimes_profile{1}; %initialize with times for first anemometer
WindEndTimes = WindEndTimes_profile{1}; %initialize with times for first anemometer
for j = 2:N_Anemometers_profile %get intersecting values for subsequent anemometers
    [WindStartTimes, WindEndTimes] = IntersectingTimeIntervals(WindStartTimes,WindEndTimes,WindStartTimes_profile{j},WindEndTimes_profile{j});
end

%% get flux/wind intersecting time intervals
[StartTimesIntersecting, EndTimesIntersecting] = IntersectingTimeIntervals(WindStartTimes,WindEndTimes,FluxStartTimes,FluxEndTimes);
[y,m,d] = ymd(StartTimesIntersecting); %get associated years / months / dates for intersecting time intervals
DatesIntersecting = datetime(y,m,d); %create list of dates for intersecting time intervals

%% get time intervals for export based on dates of interest
if use_specific_times~=1
    ind_Export = find(ismember(DatesIntersecting,Dates)); %get indices of time intervals for export
    StartTimesExport = StartTimesIntersecting(ind_Export); %get list of start times for export
    EndTimesExport = EndTimesIntersecting(ind_Export); %get list of end times for export
end
    
%% specify arguments for export depending on user choices

% %set argument for wind variables
% if use_raw_wind==1 %if using raw data
%     wind_arg = 'raw';
% else %if using interpolated data
%     wind_arg = 'int';
% end

%set arguments and labels for flux variables
if use_calibrated_flux==1 %if using calibrated fluxes
    flux_arg = 'qz';
    flux_label = 'q';
    flux_units = 'g/m2/s';
else %if using particle counts
    flux_arg = 'n';
    flux_label = 'n';
    flux_units = '(counts/interval)';
end

% go through time intervals to extract and export data
for i = 1:length(StartTimesExport)
    
    % set header name for exported data files
    name_header = [Site_Export,'_',datestr(StartTimesExport(i),'yyyymmdd_HHMM'),'_',datestr(EndTimesExport(i),'HHMM')];
    
    %% export wind variables
    for use_raw_wind = [0, 1]

        %go through each anemometer height
        for j = 1:N_Anemometers_profile

            %set argument for wind variables
            if use_raw_wind==1 %if using raw data
                wind_arg = 'raw';
            else %if using interpolated data
                wind_arg = 'int';
            end

            %extract wind variables
            [~, ~, IntervalN, IntervalInd] = ...
                ExtractVariableTimeInterval(WindData.(AnemometerName_profile{j}),...
                StartTimesExport(i),EndTimesExport(i),'u',wind_arg,wind_arg); %horizontal wind and time

            %do export if there are data in time interval
            if ~isempty(IntervalN)

                %keep only longest component
                length_IntervalInd = cellfun('length',IntervalInd);
                ind_IntervalInd = find(length_IntervalInd==max(length_IntervalInd));
                IntervalN = IntervalN(ind_IntervalInd);
                IntervalInd = IntervalInd{ind_IntervalInd};

                %get variable values
                t = WindData.(AnemometerName_profile{j})(IntervalN).t.(wind_arg)(IntervalInd); %time
                u = WindData.(AnemometerName_profile{j})(IntervalN).u.(wind_arg)(IntervalInd); %horizontal wind
                v = WindData.(AnemometerName_profile{j})(IntervalN).v.(wind_arg)(IntervalInd); %transverse wind
                w = WindData.(AnemometerName_profile{j})(IntervalN).w.(wind_arg)(IntervalInd); %vertical wind
                T = WindData.(AnemometerName_profile{j})(IntervalN).T.(wind_arg)(IntervalInd); %temperature
                z = WindData.(AnemometerName_profile{j})(IntervalN).z.z; %anemometer height

                %create wind data table
                if strcmp(WindInstrumentType,'Sonic') %include diagnostic flag if wind data is of 'Sonic' type
                    diag = double(WindData.(AnemometerName_profile{j})(IntervalN).diag.(wind_arg)(IntervalInd)); %diagnostic variable
                    wind_data_table = table(hour(t),minute(t),second(t),u,v,w,T,diag,'VariableNames',{'Hour','Minute','Second','u','v','w','T','diag'});
                else
                    wind_data_table = table(hour(t),minute(t),second(t),u,v,w,T,'VariableNames',{'Hour','Minute','Second','u','v','w','T'});           
                end

                %write table to text file
                if use_raw_wind==1 %if using raw data
                    name_wind_data = [name_header,'_Wind_',AnemometerName_profile{j},'_z',int2str(1000*z),'mm.txt'];
                else %if using interpolated data
                    name_wind_data = [name_header,'_Wind_',AnemometerName_profile{j},'_z',int2str(1000*z),'mm_Interpolated.txt'];
                end
                writetable(wind_data_table,[folder_SaveData,name_wind_data]);
            end
        end
    
        %export wind metadata values into text file
        if use_raw_wind==1 %if using raw data
            name_wind_metadata = [name_header,'_WindMetadata.txt'];
        else
            name_wind_metadata = [name_header,'_WindMetadata_Interpolated.txt'];
        end
        path_wind_metadata = [folder_SaveData,name_wind_metadata];
        fileID = fopen(path_wind_metadata,'w');
        if use_raw_wind==1 %if using raw data
            fprintf(fileID,'u (m/s) \r');
            fprintf(fileID,'v (m/s) \r');
            fprintf(fileID,'w (m/s) \r');
            fprintf(fileID,'T (C) \r');
        else
            fprintf(fileID,'u (m/s) - interpolated \r');
            fprintf(fileID,'v (m/s) - interpolated \r');
            fprintf(fileID,'w (m/s) - interpolated \r');
            fprintf(fileID,'T (C) - interpolated \r');
        end
        if strcmp(WindInstrumentType,'Sonic') %include diagnostic flag if wind data is of 'Sonic' type
            fprintf(fileID,'error diagnostic \r');
        end
    end
    
    %% export flux variables
    %find data interval
    [~, ~, IntervalN, IntervalInd] = ...
        ExtractVariableTimeInterval(WenglorData,...
        StartTimesExport(i),EndTimesExport(i),'qz','qz','t'); %get times
    
    %keep only longest component
    length_IntervalInd = cellfun('length',IntervalInd);
    ind_IntervalInd = find(length_IntervalInd==max(length_IntervalInd));
    IntervalN = IntervalN(ind_IntervalInd);
    IntervalInd = IntervalInd{ind_IntervalInd};
    
    %get values
    t = WenglorData(IntervalN).t.t(IntervalInd); %get times
    q = WenglorData(IntervalN).qz.(flux_arg)(IntervalInd,:); %get fluxes
    z = WenglorData(IntervalN).qz.z(IntervalInd,:); %get heights

    %limit to values with z>z_min
    ind_zW_usable = find(min(z)>=zW_min); %minimum Wenglor height for inclusion in calculations
    z = z(:,ind_zW_usable);
    q = q(:,ind_zW_usable);
    
    %compute total fluxes
    qz = WenglorData(IntervalN).qz.qz(IntervalInd,ind_zW_usable); %get fluxes
    n = WenglorData(IntervalN).qz.n(IntervalInd,ind_zW_usable); %get counts
    Cqn = WenglorData(IntervalN).qz.Cqn(IntervalInd,ind_zW_usable); %calibration coefficients
    sigma_Cqn = WenglorData(IntervalN).qz.sigma_Cqn(IntervalInd,ind_zW_usable); %calibration coefficients uncertainty
    sigma_n = sqrt(n); %get counts uncertainty
    sigma_qz = sqrt(sigma_Cqn.^2.*n.^2 + sigma_n.^2.*Cqn.^2); %get flux uncertainty
    [~,~,ind_BSNE] = IntervalsStartsWithin(StartTimesExport(i),EndTimesExport(i),[BSNEData.StartTime],[BSNEData.EndTime]); %get index of BSNE for zq
    if ~isempty(ind_BSNE)
        zq_BSNE = BSNEData(ind_BSNE).z.zq; %get zq for BSNE
        sigma_zq_BSNE = BSNEData(ind_BSNE).z.sigma_zq; %get sigma_zq for BSNE
    else
        zq_BSNE = zq_Site;
        sigma_zq_BSNE = sigma_zq_Site;
    end
        
    N_qz = length(qz);
    Qsum = zeros(N_qz,1);
    sigma_Qsum = zeros(N_qz,1);
    for j = 1:N_qz
        [Qsum_j,sigma_Qsum_j] = qz_summation_Wenglor(qz(j,:), z(j,:), sigma_qz(j,:), zq_BSNE, sigma_zq_BSNE);
        Qsum(j) = Qsum_j;
        sigma_Qsum(j) = sigma_Qsum_j;
    end
    
    %get information for flux table
    N_q = size(q,2); %get number of fluxes in profile
    variable_names_partialflux = cell(1,N_q); %initialize list of variable names
    table_names_partialflux = cell(1,N_q); %initialize list of value name for creating table
    for j = 1:N_q
        variable_names_partialflux{j} = [flux_label,int2str(j)]; %get variable name
        table_names_partialflux{j} = ['q(:,',int2str(j),'),']; %get value name for creating table
    end
    VariableHeaders = [{'Hour','Minute','Second'},variable_names_partialflux]; %create list of headers for table
    
    %create and export partial flux table
    TableCreationString = ['flux_data_table = table(hour(t),minute(t),second(t),',[table_names_partialflux{1,:}],'''VariableNames'',VariableHeaders);'];
    eval(TableCreationString);
    name_flux_data = [name_header,'_PartialFlux.txt'];
    writetable(flux_data_table,[folder_SaveData,name_flux_data]);

    %create and export total flux table
    totalflux_data_table = table(hour(t),minute(t),second(t),Qsum, sigma_Qsum,'VariableNames',{'Hour', 'Minute', 'Second', 'Q', 'sigma_Q'});
    name_flux_data = [name_header,'_TotalFlux.txt'];
    writetable(totalflux_data_table,[folder_SaveData,name_flux_data]);
    
    %get information for flux height table
    for j = 1:N_q
        variable_names_partialflux{j} = ['z',int2str(j)]; %get variable name
        table_names_partialflux{j} = ['z(:,',int2str(j),'),']; %get value name for creating table
    end
    VariableHeaders = [{'Hour','Minute','Second'},variable_names_partialflux]; %create list of headers for table
    
    %create flux height table
    TableCreationString = ['flux_height_table = table(hour(t),minute(t),second(t),',[table_names_partialflux{1,:}],'''VariableNames'',VariableHeaders);'];
    eval(TableCreationString);
    
    %export flux height table
    name_flux_heights = [name_header,'_FluxHeights_m.txt'];
    writetable(flux_height_table,[folder_SaveData,name_flux_heights]);
   
    %export flux metadata file
    name_flux_metadata = [name_header,'_FluxMetadata.txt'];
    path_flux_metadata = [folder_SaveData,name_flux_metadata];
    fileID = fopen(path_flux_metadata,'w');
    fprintf(fileID,'Q (g/m/s) \r'); %total flux
    fprintf(fileID,'sigma_Q (g/m/s) \r'); %total flux uncertainty
    for j = 1:N_q
        metadata_text = [flux_label,int2str(j),'(z',int2str(j),') ',flux_units,'\r'];
        fprintf(fileID,metadata_text);
    end
end