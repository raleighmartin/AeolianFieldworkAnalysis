%% initialize
clearvars;

%% information for data extraction
Sites_Export = {'Jericoacoara'}; %names of field sites for export
Dates_Export = {[...
    datetime(2014,11,13);...
    datetime(2014,11,14);...
    datetime(2014,11,18);...
    datetime(2014,11,20)]}; %dates for export
WindTypes_Export = {'Ultrasonic'}; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)

% Sites_Export = {'RanchoGuadalupe'}; %names of field sites for export
% Dates_Export = {[...
%     [datetime(2015,3,23);...
%     datetime(2015,3,24)]}; %dates for export
% WindTypes_Export = {'Ultrasonic'}; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)

% Sites_Export = {'Oceano'}; %names of field sites for export
% Dates_Export = {[...
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
%     datetime(2015,6,2);...
%     datetime(2015,6,3);...
%     datetime(2015,6,4)]}; %dates for export
% WindTypes_Export = {'Sonic'}; %choose instrument type for wind data ('Sonic' for Campbells, 'Ultrasonic' for RM Young)

use_raw_wind = 1; %use raw (non-interpolated) wind data? 0 for no, 1 for yes
use_calibrated_flux = 1; %use calibrated flux values? 0 for no (raw Wenglor counts), 1 for yes
use_time_outside_BSNE = 1; %use times outside BSNE intervals (affects calibration)? 0 for no, 1 for yes

%% folders for loading and saving data
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_WindowData = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
folder_SaveData = '../../../../Google Drive/Data/AeolianFieldwork/Export/30 Minute/'; %folder for storing outputs of this extraction
folder_Functions = '../Functions/'; %folder with functions

%% information about where to load data and save plots - unrestricted windows
LoadWindowData_Path = strcat(folder_WindowData,'DataWindowCalcs_30min_Unrestricted'); %path for loading unbinned data
ExportSummaryTable_PartialPath = strcat(folder_WindowData,'SummaryTable_30min_Unrestricted'); %path for exporting summary table

% %% information about where to load data and save plots - restricted windows
% LoadWindowData_Path = strcat(folder_WindowData,'DataWindowCalcs_30min_Restricted'); %path for loading unbinned data
% ExportSummaryTable_PartialPath = strcat(folder_WindowData,'SummaryTable_30min_Restricted'); %path for exporting summary table

%% load functions and window data
addpath(folder_Functions); %point MATLAB to location of functions
load(LoadWindowData_Path); %load data windows

%% get site info
N_Sites = length(Sites);

%% CREATE AND EXPORT SUMMARY TABLE
for i = 1:N_Sites
    StartTimes = StartTimes_all{i};
    EndTimes = EndTimes_all{i};
    tau = tauRe_all{i};
    theta = theta_all{i};
    zL = zL_all{i};
    Q = Q_all{i};
    fQ = fQ_all{i};
    
    % get predominant wind direction for each site from mean of 30-minute windows with full transport
    theta_Site = mean(theta_all{i}(fQ_all{i}==1));

    % adjust subwindow theta values according to mean for site  
    theta = theta-theta_Site;

    %generate list of variables
    TableVars = ...
        {'Date',...
        'StartTime',...
        'EndTime',...
        'tau',...
        'theta',...
        'zL',...
        'Q',...
        'fQ'};
    N_Vars = length(TableVars);

    SummaryTable = table(...
        cellstr(datestr(StartTimes,'yyyy-mm-dd')),...
        cellstr(datestr(StartTimes,'HH:MM')),...
        cellstr(datestr(EndTimes,'HH:MM')),...
        cellstr(num2str(tau,'%.3f')),...
        cellstr(num2str(theta,'%.1f')),...
        cellstr(num2str(zL,'%.3f')),...
        cellstr(num2str(Q,'%.2f')),...
        cellstr(num2str(fQ,'%.2f')),...
        'VariableNames',TableVars);

    ExportData_Path = [ExportSummaryTable_PartialPath,'_',Sites{i},'.csv'];
    writetable(SummaryTable,ExportData_Path);
end

%% go through time intervals to extract and export data

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

[~,ind_Sites] = intersect(Sites,Sites_Export); %get indices for field sites

%% go through sites for export
for i = 1:length(Sites_Export)

    %% load processed data for site, extract relevant variables
    LoadProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{ind_Sites(i)});
    load(LoadProcessedData_Path); %load processed data 
    WenglorData = ProcessedData.FluxWenglor; %Wenglor flux data
    BSNEData = ProcessedData.FluxBSNE; %BSNE data
    WindData = ProcessedData.(WindTypes_Export{i}); %data for all anemometers
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
    
    %% get start times and end times for wind observations for profile of anemometers
    AnemometerName_profile = fieldnames(WindData); %anemometer names
    N_Anemometers_profile = length(AnemometerName_profile); %number of anemometers
    
    %% figure out which time intervals to export
    ind_Export = find(ismember(Date_all{ind_Sites(i)},Dates_Export{i})); %get indices of time intervals for export
    
    %% go through time intervals
    for j = 1:length(ind_Export)
        StartTime_Export = StartTimes_all{ind_Sites(i)}(ind_Export(j));
        EndTime_Export = EndTimes_all{ind_Sites(i)}(ind_Export(j));
        
        %set header name for exported data files
        name_header = strcat(Sites_Export{i},'_',datestr(StartTime_Export,'yyyymmdd_HHMM'),'_',datestr(EndTime_Export,'HHMM'));

        %% export wind variables      
        %go through each anemometer height
        for k = 1:N_Anemometers_profile

            %extract wind variables
            [~, ~, IntervalN, IntervalInd] = ...
                ExtractVariableTimeInterval(WindData.(AnemometerName_profile{k}),...
                StartTime_Export,EndTime_Export,'u',wind_arg,wind_arg); %horizontal wind and time
            
            %do export if there are data in time interval
            if ~isempty(IntervalN)
            
                %keep only longest component
                length_IntervalInd = cellfun('length',IntervalInd);
                ind_IntervalInd = find(length_IntervalInd==max(length_IntervalInd));
                IntervalN = IntervalN(ind_IntervalInd);
                IntervalInd = IntervalInd{ind_IntervalInd};

                %get variable values
                t = WindData.(AnemometerName_profile{k})(IntervalN).t.(wind_arg)(IntervalInd); %time
                u = WindData.(AnemometerName_profile{k})(IntervalN).u.(wind_arg)(IntervalInd); %horizontal wind
                v = WindData.(AnemometerName_profile{k})(IntervalN).v.(wind_arg)(IntervalInd); %transverse wind
                w = WindData.(AnemometerName_profile{k})(IntervalN).w.(wind_arg)(IntervalInd); %vertical wind
                T = WindData.(AnemometerName_profile{k})(IntervalN).T.(wind_arg)(IntervalInd); %temperature
                z = WindData.(AnemometerName_profile{k})(IntervalN).z.z; %anemometer height

                %get diagnostic values
                if strcmp(WindTypes_Export{i},'Sonic') %directly from diagnostic column if Sonic
                    diag = double(WindData.(AnemometerName_profile{k})(IntervalN).diag.(wind_arg)(IntervalInd)); %diagnostic variable
                else %otherwise, look for error times, assign these as diag = -1
                    diag = zeros(size(u));
                    [~, ind_wind_err, ~] = intersect(t,WindData.(AnemometerName_profile{k})(IntervalN).t.err);
                    diag(ind_wind_err)=-1;                   
                end

                %create wind data table
%                 if strcmp(WindTypes_Export{i},'Sonic')
                    wind_data_table = table(hour(t),minute(t),second(t),u,v,w,T,diag,'VariableNames',{'Hour','Minute','Second','u','v','w','T','diag'});
%                 else
%                     wind_data_table = table(hour(t),minute(t),second(t),u,v,w,T,'VariableNames',{'Hour','Minute','Second','u','v','w','T'});
%                 end
                    
                %write table to text file
                name_wind_data = strcat(name_header,'_Wind_',AnemometerName_profile{k},'_z',int2str(1000*z),'mm.txt');
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
        fprintf(fileID,'error diagnostic \r');
        
        %% export flux variables
        if Q_all{ind_Sites(i)}(ind_Export(j))>0
            
            %find data interval
            [~, ~, IntervalN, IntervalInd] = ...
                ExtractVariableTimeInterval(WenglorData,...
                StartTime_Export,EndTime_Export,'qz','qz','t'); %get times
           
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
            for k = 1:N_q
                variable_names{k} = [flux_label,int2str(k)]; %get variable name
                table_names{k} = ['q(:,',int2str(k),'),']; %get value name for creating table
            end
            VariableHeaders = [{'Hour','Minute','Second'},variable_names]; %create list of headers for table

            %create flux table
            TableCreationString = ['flux_data_table = table(hour(t),minute(t),second(t),',[table_names{1,:}],'''VariableNames'',VariableHeaders);'];
            eval(TableCreationString);

            %export flux table
            name_flux_data = [name_header,'_Flux.txt'];
            writetable(flux_data_table,[folder_SaveData,name_flux_data]);

            %get information for flux height table
            for k = 1:N_q
                variable_names{k} = ['z',int2str(k)]; %get variable name
                table_names{k} = ['z(:,',int2str(k),'),']; %get value name for creating table
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
            for k = 1:N_q
                metadata_text = [flux_label,int2str(k),'(z',int2str(k),') ',flux_units,'\r'];
                fprintf(fileID,metadata_text);
            end
        end
    end
end