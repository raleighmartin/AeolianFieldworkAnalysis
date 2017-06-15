%% Master script for reading in raw logger data and producing useful interpolated data files, but without interpretation

%% 0. Initialize
clearvars;

%% 1. Information about input and output folder and files with data and metadata
folder_GeneralMetadata = '../../Metadata/'; %folder with general metadata
folder_SiteMetadata = '../../Metadata/RanchoGuadalupe/'; %folder with site metadata
folder_LoggerRawData = '../../../../Google Drive/Data/AeolianFieldwork/Raw/Datalogger/RanchoGuadalupe/'; %folder with logger raw data tables
folder_GrainSize = '../../../../Google Drive/Data/AeolianFieldwork/Raw/GrainSize/'; %folder with grain size data
folder_DataOutput = '../../../../Google Drive/Data/AeolianFieldwork/Processed/RanchoGuadalupe/'; %folder for storing data output
folder_Functions = '../Functions/'; %folder with functions

addpath(folder_Functions); %point MATLAB to location of functions

file_InstrumentCalibration = 'InstrumentCalibration.xlsx'; %file with instrument calibration values
file_LoggerTables = 'LoggerTables_RanchoGuadalupe.xlsx'; %file with logger tables
file_LoggerTimes = 'LoggerTimes_RanchoGuadalupe.xlsx'; %file with logger times
file_InstrumentVariables = 'InstrumentVariables_RanchoGuadalupe.xlsx'; %file with types of instruments and variables
file_InstrumentMetadata = 'InstrumentMetadata_RanchoGuadalupe.xlsx'; %file with logger tables
file_WeightBSNE = 'WeightBSNE_RanchoGuadalupe.xlsx'; %file with BSNE weights and times
file_GrainSizeMetadata = 'GrainSizeMetadata_RanchoGuadalupe.xlsx'; %file with grain size metadata information

Metadata_Path = strcat(folder_DataOutput,'Metadata_RanchoGuadalupe'); %get path to saving metadata
RawData_Path = strcat(folder_DataOutput,'RawData_RanchoGuadalupe'); %get path to saving raw data
InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_RanchoGuadalupe'); %get path to saving interpolated data
GrainSize_Path = strcat(folder_DataOutput,'GrainSize_RanchoGuadalupe'); % get path to saving grain size data

%% 2. Parse spreadsheets for info about loggers and site
LoggerTables = ParseLoggerTables(folder_SiteMetadata,file_LoggerTables); %extract info about logger tables from .xlsx file
LoggerTimes = ParseLoggerTimes(folder_SiteMetadata,file_LoggerTimes); %extract info about logger start/end times from .xlsx file
InstrumentVariables = ParseInstrumentVariables(folder_SiteMetadata,file_InstrumentVariables); %!!extract info about types of instruments and associated variables from .xlsx file
InstrumentMetadata = ParseInstrumentMetadata(folder_SiteMetadata,file_InstrumentMetadata); %extract info about instruments from .xlsx file
InstrumentCalibration = ParseInstrumentCalibration(folder_GeneralMetadata,file_InstrumentCalibration); %!!extract info about instrument calibrations from .xlsx file
WeightBSNE = ParseBSNE(folder_SiteMetadata,file_WeightBSNE); %extract info about BSNE weights from .xlsx spreadsheet
[GrainSizeMetadata_Surface, GrainSizeMetadata_BSNE] = ParseGrainSizeMetadata(folder_SiteMetadata,file_GrainSizeMetadata); %extract info about instruments from .xlsx file

%% 3. Import and aggregate data from grain size files
GrainSize_Surface = GetGrainSize(GrainSizeMetadata_Surface,folder_GrainSize);
GrainSize_BSNE = GetGrainSize(GrainSizeMetadata_BSNE,folder_GrainSize);

%% 4. Save all metadata together
save(Metadata_Path,'LoggerTables','LoggerTimes','InstrumentVariables',...
   'InstrumentMetadata','InstrumentCalibration','WeightBSNE',...
   'GrainSize_Surface','GrainSize_BSNE',...
   'GrainSizeMetadata_Surface','GrainSizeMetadata_BSNE','-v7.3');

%% 5. Import and aggregate data from logger files
RawData = GetRawData(LoggerTables,LoggerTimes,InstrumentVariables,InstrumentMetadata,folder_LoggerRawData);
save(RawData_Path,'RawData','-v7.3'); %save raw data

%% 6. Perform interpolation on raw data
InterpolatedData = InterpolateData(RawData, InstrumentVariables); %interpolate raw data, flag errors

%% 7. Apply calibration factor to data based on unique identifier
InterpolatedData = ParseCalibrationFactor(InterpolatedData,InstrumentVariables,InstrumentCalibration);
save(InterpolatedData_Path,'InterpolatedData','-v7.3'); %save data

%% restore function path to default value
restoredefaultpath;