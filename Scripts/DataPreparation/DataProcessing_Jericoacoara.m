% Script for taking interpolated data and making useful basic calculations
% for subsequent analysis, including fitting profiles to BSNEs

%% 0. Initialize
clearvars;

%% 1. Information about input and output folder and files with data and metadata, load these
folder_DataInput = '../../../../Google Drive/Data/AeolianFieldwork/Processed/Jericoacoara/'; %folder for loading data input
folder_Functions = '../Functions/'; %folder with functions
% folder_ProcessingFunctions = '../../AeolianFieldworkAnalysis/Scripts/Processing/'; %folder for data processing functions
% folder_GeneralFunctions = '../../AeolianFieldworkAnalysis/Scripts/Functions/'; %folder with functions
folder_DataOutput = '../../../../Google Drive/Data/AeolianFieldwork/Processed/Jericoacoara/'; %folder for storing data output
% folder_DataOutput = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output

addpath(folder_Functions); %point MATLAB to location of processing functions
% addpath(folder_ProcessingFunctions); %point MATLAB to location of processing functions
% addpath(folder_GeneralFunctions); %point MATLAB to location of general functions

file_Metadata = 'Metadata_Jericoacoara.mat'; %file with aggregated metadata for site
file_InterpolatedData = 'InterpolatedData_Jericoacoara.mat'; %file with aggregated interpolated data for site
Metadata_Path = strcat(folder_DataInput,file_Metadata); %get path to metadata
InterpolatedData_Path = strcat(folder_DataInput,file_InterpolatedData); %get path to interpolated data
% Metadata_Path = strcat(folder_DataOutput,'Metadata_Jericoacoara'); %get path to metadata
% InterpolatedData_Path = strcat(folder_DataOutput,'InterpolatedData_Jericoacoara'); %get path to interpolated data

load(Metadata_Path); %load metadata
load(InterpolatedData_Path); %load interpolated data

file_BSNEData = 'FluxBSNE_Jericoacoara.mat'; %filename for saving BSNE data
file_ProcessedData = 'ProcessedData_Jericoacoara.mat'; %filename for saving processed data
BSNEData_Path = strcat(folder_DataOutput,file_BSNEData); %path for saving processed data
ProcessedData_Path = strcat(folder_DataOutput,file_ProcessedData); %create path for saving BSNE data
% BSNEData_Path = strcat(folder_DataOutput,'FluxBSNE_Jericoacoara'); %path for saving output data
% ProcessedData_Path = strcat(folder_DataOutput,'ProcessedData_Jericoacoara'); %create path to processed data

zq_estimated_Site_m = 0.097; %estimated saltation layer height (m) for site (for BSNE flux uncertainty estimation - from Martin & Kok 2017, Science Advances)

%% 2. Process instrument heights
ProcessedData = ProcessInstrumentHeights(InterpolatedData); %run function to process instrument heights, output generates "ProcessedData" structured array

%% 3. Process BSNE profiles
FluxBSNE = ProcessBSNEs(WeightBSNE,GrainSize_BSNE,zq_estimated_Site_m); %run function to process BSNE information, generate BSNE fluxes
save(BSNEData_Path,'FluxBSNE'); %save BSNE data as its own file
ProcessedData.FluxBSNE = FluxBSNE; %also add BSNE flux data to ProcessedData structured array

%% 4. Process Wenglors
[ProcessedWenglors, FluxWenglor] = ProcessWenglors(ProcessedData, FluxBSNE, InstrumentMetadata);
ProcessedData.Wenglor = ProcessedWenglors;
ProcessedData.FluxWenglor = FluxWenglor;

% 5. Move extraneous metadata fields out of file
ProcessedData = MoveExtraneousMetadataFields(ProcessedData);

% 6. Modify wind velocities so that u is oriented with setup (Jeri only)
Anemometers = fieldnames(ProcessedData.Ultrasonic);
N_Anemometers = length(Anemometers);
for i = 1:N_Anemometers
    U_old = ProcessedData.Ultrasonic.(Anemometers{i}); %get old anemometer data
    U_new = U_old; %initialize new as old
    N_windows = length(U_old); %go through each time interval
    for j = 1:N_windows
        U_new(j).u = U_old(j).v; %replace u with v
        U_new(j).v = U_old(j).u; %replace v with u
        U_new(j).v.raw = -U_new(j).v.raw; %take negative for v
        U_new(j).v.int = -U_new(j).v.int; %take negative for v
    end
    ProcessedData.Ultrasonic.(Anemometers{i}) = U_new; %replace new U in Processed Data
end
    
%% 7. Save processed data
save(ProcessedData_Path,'ProcessedData','-v7.3'); %save data

%% Restore function path to default value
restoredefaultpath;