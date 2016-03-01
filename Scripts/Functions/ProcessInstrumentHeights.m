%% function to determine absolute heights for instruments
 
function DataWithHeights = ProcessInstrumentHeights(Data)

%ignore baseline uncertainty in reference height, becuase this is also captured in other uncertainties... but could change later
RefHeightErr_baseline_m = 0;

%% GO THROUGH INSTRUMENT TYPES, DETERMINE ABSOLUTE HEIGHTS BASED ON DISTANCE SENSOR READINGS
InstrumentTypes = fieldnames(Data); %get list of instrument types
N_InstrumentTypes = length(InstrumentTypes); %how many instrument types

for i = 1:N_InstrumentTypes
    InstrumentType = InstrumentTypes{i}; %get current instrument type 
    Instruments = fieldnames(Data.(InstrumentType)); %get list of instruments
    N_Instruments = length(Instruments); %how many instruments of this type
    
    for j = 1:N_Instruments
        Instrument = Instruments{j} %get current instrument
        N_Intervals = length(Data.(InstrumentType).(Instrument)); %how many time intervals for this instrument
        RefHeight_List = zeros(N_Intervals,1); %create list of reference heights (m)
        RefHeightErr_List = zeros(N_Intervals,1); %create list of reference heights errors (m)
        
        for k = 1:N_Intervals
            StartTime = Data.(InstrumentType).(Instrument)(k).StartTime; %get start time
            EndTime = Data.(InstrumentType).(Instrument)(k).EndTime; %get end time
            StartHeight_m = Data.(InstrumentType).(Instrument)(k).StartHeight_m; %get start height
            EndHeight_m = Data.(InstrumentType).(Instrument)(k).EndHeight_m; %get end height
            HeightRef = Data.(InstrumentType).(Instrument)(k).HeightRef; %get reference variable for height
            HeightErr_m = Data.(InstrumentType).(Instrument)(k).HeightErr_m; %get height error for instrument
            
            %compute reference height
            %if it is 0, RefHeight_m and RefHeightErr_m are 0
            if strcmp(HeightRef,'0')
                RefHeight_m = 0;
                RefHeightErr_m = 0;
            
            %otherwise extract data from distance sensor specified by this string
            else
               
                %extract distance sensor heights (mm) for this interval
                z_interval_mm = ExtractVariableTimeInterval(Data.Distance.(HeightRef),StartTime,EndTime,'z','int','int');
                
                %if values are contained in this interval, compute mean of them 
                if ~isempty(z_interval_mm)
                    RefHeight_m = mean(z_interval_mm)/1000; %compute mean, convert to meters
                    RefHeightErr_sensor_m = std(z_interval_mm)/1000; %sensor error is bed elevation standard deviation, convert to meters
                    RefHeightErr_m = sqrt(RefHeightErr_baseline_m.^2+RefHeightErr_sensor_m.^2); %total error from sensor fluctuations and baseline error
                    
                %if no values contained in interval, use values from previous interval
                else
                    RefHeight_m = RefHeight_List(k-1);
                    RefHeightErr_m = RefHeightErr_List(k-1);
                end
            end

            %assign values to lists
            RefHeight_List(k) = RefHeight_m;
            RefHeightErr_List(k) = RefHeightErr_m;
            
            %compute instrument height (gives mean height for the entire interval)
            z_Instrument = RefHeight_m+mean([StartHeight_m; EndHeight_m]);
            
            %compute error as combination of measurement error, reference
            %height error, and difference of start / end heights
            sigma_z_Instrument = sqrt(HeightErr_m^2 + RefHeightErr_m^2 + (StartHeight_m-EndHeight_m).^2);
            
            %assign height to structured array
            Data.(InstrumentType).(Instrument)(k).z.z = z_Instrument;
            Data.(InstrumentType).(Instrument)(k).z.sigma_z = sigma_z_Instrument;
            Data.(InstrumentType).(Instrument)(k).z.Units = 'm';
        end
    end
end

%% now, rename "Data" as "DataWithHeights"
DataWithHeights = Data;