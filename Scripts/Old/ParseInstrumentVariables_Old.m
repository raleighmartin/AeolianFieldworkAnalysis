%% Function to read in metadata from 'file_InstrumentVariables' (a .xlsx spreadsheet) that
% describes information about the variables associated with each instrument
% Function outputs a structured array with elements from spreadsheet columns
% Dependencies: NONE
% Used by: Processing_Master

function InstrumentVariables = ParseInstrumentVariables(folder_InstrumentMetadata,file_InstrumentVariables)
    
%% Get file location
filePath = strcat(folder_InstrumentMetadata,file_InstrumentVariables);

%% Import the data
[~, ~, raw] = xlsread(filePath,'Sheet1');
raw = raw(2:end,:);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Allocate imported array to column variable names
InstrumentType = {raw(:,1)};
Instrument = {raw(:,2)};
VarNameSpecific = {raw(:,3)};
VarNameGeneral = {raw(:,4)};
DataType = {raw(:,5)};
Units = {raw(:,6)};
Calibration = {raw(:,7)};
ErrorValues = {raw(:,8)};

%convert error values to vectors
ErrorValues = ErrorValues{1};
for i = 1:length(ErrorValues)
    if isstr(ErrorValues{i})
        ErrorValues_split = strsplit(ErrorValues{i},',');
        N_ErrorValues = length(ErrorValues_split);
        ErrorValues{i} = zeros(N_ErrorValues,1);
        for j = 1:N_ErrorValues
            ErrorValues{i}(j) = str2num(ErrorValues_split{j});
        end
    end
end
ErrorValues = {ErrorValues};

%% Create structured array with spreadsheet information
InstrumentVariables = struct(...
    'InstrumentType',InstrumentType,...
    'Instrument',Instrument,...
    'VarNameSpecific',VarNameSpecific,...
    'VarNameGeneral',VarNameGeneral,...
    'DataType',DataType,...
    'Units',Units,...
    'Calibration',Calibration,...
    'ErrorValues',ErrorValues);