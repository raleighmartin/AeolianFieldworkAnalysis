%% SCRIPT TO EXPORT SALTATION AND WIND DATA INTO ASCII FILES FOR ANALYSIS OUTSIDE MATLAB
% FUNCTION DEPENDENCIES: IntersectingTimeIntervals.m, ExtractVariableTimeInterval.m

%% initialize
clearvars;

%% information about where to load/save data
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Export/Daily/'; %folder for storing outputs of this extraction

%info about which data to use
use_raw_wind = 1; %use raw (non-interpolated) wind data? 0 for no, 1 for yes
use_calibrated_flux = 1; %use calibrated flux values? 0 for no (raw Wenglor counts), 1 for yes
use_time_outside_BSNE = 1; %use times outside BSNE intervals (affects calibration)? 0 for no, 1 for yes

%% information for full data extraction - Jericoacoara
Site = 'Jericoacoara'; %name of field site
Dates = [datetime(2014,11,13);...
    datetime(2014,11,14);...
    datetime(2014,11,18);...
    datetime(2014,11,20)]; %dates for data extraction
WindInstrumentType = 'Ultrasonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)

% %% information for full data extraction - Rancho Guadalupe
% Site = 'RanchoGuadalupe'; %name of field site
% Dates = [datetime(2015,3,23);...
%     datetime(2015,3,24)]; %dates for data extraction
% WindInstrumentType = 'Ultrasonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)

% %% information for full data extraction - Oceano
% Site = 'Oceano'; %name of field site
% Dates = [datetime(2015,5,15);...
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
%     datetime(2015,6,2);...
%     datetime(2015,6,3);...
%     datetime(2015,6,4)]; %dates for data extraction
% WindInstrumentType = 'Sonic'; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% load processed data and extract relevant variables
LoadData_Path = strcat(folder_LoadData,'ProcessedData_',Site);
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
ind_Export = find(ismember(DatesIntersecting,Dates)); %get indices of time intervals for export
StartTimesExport = StartTimesIntersecting(ind_Export); %get list of start times for export
EndTimesExport = EndTimesIntersecting(ind_Export); %get list of end times for export

%% specify arguments for export depending on user choices

%set argument for wind variables
if use_raw_wind==1 %if using raw data
    wind_arg = 'raw';
else %if using interpolated data
    wind_arg = 'int';
end

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

%% go through time intervals to extract and export data
for i = 1:length(ind_Export)

    %set header name for exported data files
    name_header = [Site,'_',datestr(StartTimesExport(i),'yyyymmdd_HHMM'),'_',datestr(EndTimesExport(i),'HHMM')];
    
    %% export wind variables
    
    %go through each anemometer height
    for j = 1:N_Anemometers_profile
        
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
            name_wind_data = [name_header,'_Wind_',AnemometerName_profile{j},'_z',int2str(1000*z),'mm.txt'];
            writetable(wind_data_table,[folder_SaveData,name_wind_data]);
        end
    end
    
    %export wind metadata values into text file
    name_wind_metadata = [name_header,'_WindMetadata.txt'];
    path_wind_metadata = [folder_SaveData,name_wind_metadata];
    fileID = fopen(path_wind_metadata,'w');
    fprintf(fileID,'u (m/s) \r');
    fprintf(fileID,'v (m/s) \r');
    fprintf(fileID,'w (m/s) \r');
    fprintf(fileID,'T (C) \r');
    if strcmp(WindInstrumentType,'Sonic') %include diagnostic flag if wind data is of 'Sonic' type
        fprintf(fileID,'error diagnostic \r');
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

    %get information for flux table
    N_q = size(q,2); %get number of fluxes in profile
    variable_names = cell(1,N_q); %initialize list of variable names
    table_names = cell(1,N_q); %initialize list of value name for creating table
    for j = 1:N_q
        variable_names{j} = [flux_label,int2str(j)]; %get variable name
        table_names{j} = ['q(:,',int2str(j),'),']; %get value name for creating table
    end
    VariableHeaders = [{'Hour','Minute','Second'},variable_names]; %create list of headers for table
    
    %create flux table
    TableCreationString = ['flux_data_table = table(hour(t),minute(t),second(t),',[table_names{1,:}],'''VariableNames'',VariableHeaders);'];
    eval(TableCreationString);
    
    %export flux table
    name_flux_data = [name_header,'_Flux.txt'];
    writetable(flux_data_table,[folder_SaveData,name_flux_data]);
    
    %get information for flux height table
    for j = 1:N_q
        variable_names{j} = ['z',int2str(j)]; %get variable name
        table_names{j} = ['z(:,',int2str(j),'),']; %get value name for creating table
    end
    VariableHeaders = [{'Hour','Minute','Second'},variable_names]; %create list of headers for table
    
    %create flux height table
    TableCreationString = ['flux_height_table = table(hour(t),minute(t),second(t),',[table_names{1,:}],'''VariableNames'',VariableHeaders);'];
    eval(TableCreationString);
    
    %export flux height table
    name_flux_heights = [name_header,'_FluxHeights_m.txt'];
    writetable(flux_height_table,[folder_SaveData,name_flux_heights]);
   
    %export flux metadata file
    name_flux_metadata = [name_header,'_FluxMetadata.txt'];
    path_flux_metadata = [folder_SaveData,name_flux_metadata];
    fileID = fopen(path_flux_metadata,'w');
    for j = 1:N_q
        metadata_text = [flux_label,int2str(j),'(z',int2str(j),') ',flux_units,'\r'];
        fprintf(fileID,metadata_text);
    end
end