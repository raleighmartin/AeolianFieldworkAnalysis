%% initialize
clearvars;

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,1); %sampling interval for analysis

%% folders for loading and saving data
folder_WindowData = '../../AnalysisData/Windowing/'; %folder for retrieving processed data
folder_AnalysisData = '../../AnalysisData/Thresholds/'; %folder for retrieving analysis data
folder_Functions = '../Functions/'; %folder with functions

%% information about where to load data and save plots
LoadWindowData_Path = strcat(folder_WindowData,'DataWindowCalcs_30min_Unrestricted'); %path for loading unbinned data
LoadSubwindowData_Path = strcat(folder_WindowData,'DataSubwindowCalcs_30min_Unrestricted'); %path for loading unbinned data
LoadThresholdData_Path = strcat(folder_AnalysisData,'ThresholdAnalysisData'); %path for loading binned threshold data
LoadPredictionData_Path = strcat(folder_AnalysisData,'ActivityPredictionData'); %path for loading predicted activities
ExportData_PartialPath = strcat(folder_AnalysisData,'ThresholdData'); %path for saving unbinned data

%% load data and point to functions
load(LoadWindowData_Path); %load data windows
load(LoadSubwindowData_Path); %load data subwindows
load(LoadThresholdData_Path); %load threshold data
load(LoadPredictionData_Path); %load activity prediction data
addpath(folder_Functions); %point MATLAB to location of functions

%% get indices for timescale
ind_deltat = find(deltat_all==deltat_analysis);
ind_Deltat = find(Deltat_all==Deltat_analysis);

%% go through sites
N_Sites = length(Sites);

%% CREATE TABLE OF VALUES AT EACH SITE
for i = 1:N_Sites
    
    % get data - unbinned
    StartTime = StartTime_subwindow_all{i}{ind_Deltat,ind_deltat};
    EndTime = EndTime_subwindow_all{i}{ind_Deltat,ind_deltat};
    fQ = fQ_subwindow_all{i}{ind_Deltat,ind_deltat};
    uth = uth_subwindow_all{i}{ind_Deltat,ind_deltat};
    fQpred_ft = fQpred_ft_all{i};
    fQpred_it = fQpred_it_all{i};
    fQpred_hyst = fQpred_hyst_all{i};
    
    % get predominant wind direction for each site from mean of 30-minute windows with full transport
    theta_Site = mean(theta_all{i}(fQ_all{i}==1));
    
    % adjust subwindow theta values according to mean for site  
    theta = theta_subwindow_all{i}{ind_Deltat,ind_deltat}-theta_Site;
    
    % get z/L from associated 30-minute window
    zL = zL_all{i}(ind_window_subwindow_all{i}{ind_Deltat,ind_deltat});
   
    % for Jeri, convert fQ NaN values to zeros, then sort values by date
    if strcmp(Sites{i},'Jericoacoara')
        fQ(isnan(fQ))=0;
        [~, ind_sort] = sort(StartTime);
        StartTime = StartTime(ind_sort);
        EndTime = EndTime(ind_sort);
        fQ = fQ(ind_sort);
        uth = uth(ind_sort);
        zL = zL(ind_sort);
        theta = theta(ind_sort);
        fQpred_ft = fQpred_ft(ind_sort);
        fQpred_it = fQpred_it(ind_sort);
        fQpred_hyst = fQpred_hyst(ind_sort);
    else % otherwise, remove fQ NaN values
        ind_defined = find(~isnan(fQ));
        StartTime = StartTime(ind_defined);
        EndTime = EndTime(ind_defined);
        fQ = fQ(ind_defined);
        uth = uth(ind_defined);
        zL = zL(ind_defined);
        theta = theta(ind_defined);
        fQpred_ft = fQpred_ft(ind_defined);
        fQpred_it = fQpred_it(ind_defined);
        fQpred_hyst = fQpred_hyst(ind_defined);
    end
       
    %generate list of variables
    TableVars = cell(1,10);
    TableVars{1,1} = 'Date';
    TableVars{1,2} = 'StartTime';
    TableVars{1,3} = 'EndTime';
    TableVars{1,4} = 'f_Q';
    TableVars{1,5} = 'uth';
    TableVars{1,6} = 'f_Q_ft';
    TableVars{1,7} = 'f_Q_it';
    TableVars{1,8} = 'f_Q_dual';
    TableVars{1,9} = 'zL';
    TableVars{1,10} = 'theta';

    ThresholdTable = table(...
        cellstr(datestr(StartTime,'yyyy-mm-dd')),...
        cellstr(datestr(StartTime,'HH:MM')),...
        cellstr(datestr(EndTime,'HH:MM')),...
        cellstr(num2str(fQ,'%.2f')),...
        cellstr(num2str(uth,'%.2f')),...
        cellstr(num2str(fQpred_ft,'%.2f')),...
        cellstr(num2str(fQpred_it,'%.2f')),...
        cellstr(num2str(fQpred_hyst,'%.2f')),...
        cellstr(num2str(zL,'%.3f')),...
        cellstr(num2str(theta,'%.1f')),...
        'VariableNames',TableVars);

    ExportData_Path = [ExportData_PartialPath,'_',Sites{i},'.csv'];
    writetable(ThresholdTable,ExportData_Path);
end