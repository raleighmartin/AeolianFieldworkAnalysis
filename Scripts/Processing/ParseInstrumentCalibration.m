%% Function to read in metadata from 'file_InstrumentVariables' (a .xlsx spreadsheet) that
% describes information about the variables associated with each instrument
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: NONE
% Used by: Processing_Master

function InstrumentCalibration = ParseInstrumentCalibration(folder_InstrumentMetadata,file_InstrumentCalibration)
    
%% Get file location
filePath = strcat(folder_InstrumentMetadata,file_InstrumentCalibration);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
InstrumentType = {raw(:,1)};
InstrumentID = {raw(:,2)};
VarNameGeneral = {raw(:,3)};
CalibrationFactor = cell2mat(raw(:,4));
CalibrationIntercept = cell2mat(raw(:,5));

%% Create structured array with spreadsheet information
InstrumentCalibration = struct(...
    'InstrumentType',InstrumentType,...
    'InstrumentID',InstrumentID,...
    'VarNameGeneral',VarNameGeneral,...
    'CalibrationFactor',CalibrationFactor,...
    'CalibrationIntercept',CalibrationIntercept);