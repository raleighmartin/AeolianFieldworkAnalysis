%% 4. Move extraneous fields into "AdditionalMetadata"
function [Data_new] = MoveExtraneousMetadataFields(Data)

ExtraneousMetadataFields = {'InstrumentType','InstrumentID','StartHeight_m',...
    'EndHeight_m','HeightErr_m','HeightRef','Longitudinal_m','Spanwise_m',...
    'AngleErr_deg','ErrorCode'};
N_ExtraneousMetadataFields = length(ExtraneousMetadataFields);
Data_new = Data;

%get instrument types, remove from consideration fields that are not actually instrument types (e.g. flux calcs)
InstrumentTypes = fieldnames(Data);
ExcludedInstrumentTypes = {'FluxBSNE','FluxWenglor'};
InstrumentTypes = setdiff(InstrumentTypes,ExcludedInstrumentTypes);
N_InstrumentTypes = length(InstrumentTypes);

%go through each instrument interval
for i = 1:N_InstrumentTypes
    Instruments = fieldnames(Data.(InstrumentTypes{i}));
    N_Instruments = length(Instruments);
    for j = 1:N_Instruments
        N_Intervals = length(Data.(InstrumentTypes{i}).(Instruments{j}));
        for k = 1:N_Intervals
            Data_new.(InstrumentTypes{i}).(Instruments{j})(k).AdditionalMetadata = struct(); %initialize "additional metadata" structured array within instrument interval
            
            %go through each extraneous metadata field
            for l = 1:N_ExtraneousMetadataFields
                Data_new.(InstrumentTypes{i}).(Instruments{j})(k).AdditionalMetadata.(ExtraneousMetadataFields{l})=...
                    Data.(InstrumentTypes{i}).(Instruments{j})(k).(ExtraneousMetadataFields{l});
            end
        end
        %remove extraneous fields
        Data_new.(InstrumentTypes{i}).(Instruments{j})=rmfield(Data_new.(InstrumentTypes{i}).(Instruments{j}),ExtraneousMetadataFields);
    end
end