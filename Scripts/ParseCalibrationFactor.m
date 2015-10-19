%% Function to read in 'Data' and add calibration values based on
% 'InstrumentVariables' and 'InstrumentCalibration'
% Dependencies: NONE
% Used by: DataExtraction_Jericoacoara2014

function DataWithCalibration = ParseCalibrationFactor(Data,InstrumentVariables,InstrumentCalibration)

InstrumentTypes = fieldnames(Data);
N_InstrumentTypes = length(InstrumentTypes);
for i = 1:N_InstrumentTypes
    Instruments = fieldnames(Data.(InstrumentTypes{i}));
    N_Instruments = length(Instruments);
    for j = 1:N_Instruments
        N_Intervals = length(Data.(InstrumentTypes{i}).(Instruments{j}));
        Variables = InstrumentVariables.VarNameGeneral(intersect(...
            find(strcmp(InstrumentVariables.InstrumentType,InstrumentTypes{i})),...
            find(strcmp(InstrumentVariables.Instrument,Instruments{j}))));
        N_Variables = length(Variables);
        for k = 1:N_Intervals
            InstrumentID = Data.(InstrumentTypes{i}).(Instruments{j})(k).InstrumentID;
            for l = 1:N_Variables
                if InstrumentVariables.Calibration(intersect(...
                        find(strcmp(InstrumentVariables.Instrument,Instruments{j})),...
                        find(strcmp(InstrumentVariables.VarNameGeneral,Variables{l}))))==1
                    Data.(InstrumentTypes{i}).(Instruments{j})(k).(Variables{l}).CalibrationFactor = ...
                        InstrumentCalibration.CalibrationFactor(intersect(...
                        find(strcmp(InstrumentCalibration.InstrumentID,InstrumentID)),...
                        find(strcmp(InstrumentCalibration.VarNameGeneral,Variables{l}))));
                    Data.(InstrumentTypes{i}).(Instruments{j})(k).(Variables{l}).CalibrationIntercept = ...
                        InstrumentCalibration.CalibrationIntercept(intersect(...
                        find(strcmp(InstrumentCalibration.InstrumentID,InstrumentID)),...
                        find(strcmp(InstrumentCalibration.VarNameGeneral,Variables{l}))));
                end
            end
        end
    end
end

DataWithCalibration = Data;