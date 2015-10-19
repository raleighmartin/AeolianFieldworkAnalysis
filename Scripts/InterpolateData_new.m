%% function to interpolate and flag error values in raw data
% 'RawData' contains raw observations grouped by instrument
% Output 'InterpolatedData' is built on 'RawData', but adds fields to account for interpolations 
% Dependencies: NONE
% Used by: Processing_Master
 
function [InterpolatedData] = InterpolateData(RawData, InstrumentVariables)

%go through instrument types, first perform basic processing common to all instruments
InstrumentTypes = fieldnames(RawData); %get list of instrument types
N_InstrumentTypes = length(InstrumentTypes); %how many instrument types

for i = 1:N_InstrumentTypes
    InstrumentType = InstrumentTypes{i}; %get current instrument type 
    Instruments = fieldnames(RawData.(InstrumentType)); %get list of instruments
    N_Instruments = length(Instruments); %how many instruments of this type
    
    for j = 1:N_Instruments
        Instrument = Instruments{j} %get current instrument
        N_Intervals = length(RawData.(InstrumentType).(Instrument)); %how many time intervals for this instrument
        
        Variables = unique(InstrumentVariables.VarNameGeneral(strcmp(InstrumentVariables.Instrument,Instrument))); %get variables for instrument
        N_Variables = length(Variables);
        
        for k = 1:N_Intervals
            t_raw = RawData.(InstrumentType).(Instrument)(k).t.raw; %get raw times
            
            %deal with repeated times
            [t_unique, ind_t_unique, ind_t_all] = unique(t_raw); %find unique times and indices
            t_repeat = t_raw(ind_t_all(diff(ind_t_all))); %get list of repeat times
            
            %deal with time gaps
            dt = mode(diff(t_unique)); %get expected timestep
            StartTime = min(t_unique); %initial time
            EndTime = max(t_unique); %final time
            t_int = (StartTime:dt:EndTime)'; %create list of interpolated times
            N_t = length(t_int); %number of interpolated timesteps
            [~, ind_nongap] = intersect(t_int, t_unique); %record indices of nongap times
            [t_gap, ind_gap] = setdiff(t_int,t_raw); %get times and indices of gap times
            
            %now, interpolate individual variables
            for l = 1:N_Variables
                values_raw = RawData.(InstrumentType).(Instrument)(k).(Variables{l}).raw; %get raw values
                values_t_unique = values_raw(ind_t_unique); %remove repeats
                values_int = zeros(N_t)*NaN; %initialize vector of interpolated values as list of NaNs
                values_int(ind_nongap) = values_t_unique; %use values from associated with unique times for nongap locations
                good_ind = intersect(find(values_int~=-999),find(isnan(values_int)==0)); %points without -999 or NaN are good
                error_ind = union(find(values_int==-999),find(isnan(values_int)==1)); %find error points based on -999 or NaN
                N_error = length(error_ind); %how many error points?
                t_err = [t_err; t_int(error_ind)]; %add to list of error times
                t_err = unique(t_err); %keep only unique times in list
                
                %get previous indices for interpolation
                prev_ind = error_ind-1; %take previous indices for interpolation to be error indices minus one
                [prev_error_intersect, prev_error_intersect_i] = intersect(prev_ind,error_ind);
                while(length(prev_error_intersect)>0); %run loop if there is intersection
                    prev_ind(prev_error_intersect_i) = prev_ind(prev_error_intersect_i)-1;
                    [prev_error_intersect, prev_error_intersect_i] = intersect(prev_ind,error_ind);
                end
                prev_ind(prev_ind<1)=1; %prev_ind cannot be less than one
                
                %get next indices for interpolation
                next_ind = error_ind+1; %take next indices for interpolation to be error indices plus one
                [next_error_intersect, next_error_intersect_i] = intersect(next_ind,error_ind);
                while(~isempty(next_error_intersect)); %run loop if there is intersection
                    next_ind(next_error_intersect_i) = next_ind(next_error_intersect_i)+1;
                    [next_error_intersect, next_error_intersect_i] = intersect(next_ind,error_ind);
                end
                next_ind(next_ind>N_t)=N_t; %next_ind cannot be greater than length of timeseries
                
                %set first and last interpolation values
                values_int(1) = values_int(min(good_ind)); %set first point 
                values_int(end) = values_int(max(good_ind)); %set last point

                %make calculations for interpolation values and weighting
                prev_value = values_int(prev_ind);
                next_value = values_int(next_ind);
                prev_weight = (next_ind-error_ind)./(next_ind-prev_ind);
                next_weight = (error_ind-prev_ind)./(next_ind-prev_ind);
                
                %set remaining interpolation values based on prev/next values and their weights
                values_int(error_ind)=prev_value.*prev_weight+next_value.*next_weight;
                
%                 %go through each error point and interpolate
%                 values_int(1) = values_int(min(good_ind)); %set first point 
%                 values_int(end) = values_int(max(good_ind)); %set last point
%                 for m = 1:N_error
%                     ind = error_ind(m);
%                     if ind~=1&&ind~=length(values_int) %ignore first and last point (set above)
%                         ind_prev = max([good_ind(good_ind<ind);1]); %get index of previous good point (or 1)
%                         ind_next = min([good_ind(good_ind>ind);length(values_int)]); %get index of next good point (or last point)
%                         value_prev = values_int(ind_prev); %get value of last good point
%                         value_next = values_int(ind_next); %get value of next good point
%                         values_int(ind) = ... %perform interpolation
%                             value_prev*(ind_next-ind)/(ind_next-ind+1)+...
%                             value_next/(ind_next-ind+1);
%                     end
%                 end
                RawData.(InstrumentType).(Instrument)(k).(Variables{l}).int = values_int; %add interpolation field to variable values
            end
            RawData.(InstrumentType).(Instrument)(k).t = struct('raw',t_raw,'int',t_int,'err',t_err); %record interpolated times (with raw and error times) in structured array
        end 
    end
end
    
%now that everything is done, rename "RawData" as "InterpolatedData"
InterpolatedData = RawData;