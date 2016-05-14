%% SCRIPT TO TO ESTIMATE RIPPLE MIGRATION FLUX
% 
%% initialize
clearvars;
close all;

%% parameter values
phi = 0.65; %volume fraction of soil

%% information about field sites and distance sensors
Sites = {'RanchoGuadalupe';'Oceano'};
UpwindSensor = {'L1';'L2'};
DownwindSensor = {'L2';'L1'};
delta_x_Sensor = [0.08 0.06]; %distance separation of upwind/downwind distance sensors (m)
N_Sites = length(Sites);

%% information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../AnalysisData/'; %folder for storing outputs of this analysis
SaveFullData_Path = strcat(folder_ProcessedData,'StressFluxWindows_all');
SaveData_Path = strcat(folder_AnalysisData,'StressFluxWindows_all');

%% load processed and metadata for each site, add to structured arrays of all data and metadata
Data = cell(N_Sites,1);
Metadata = cell(N_Sites,1);
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    Metadata_Path = strcat(folder_ProcessedData,'Metadata_',Sites{i});
    Data{i} = load(ProcessedData_Path); %load processed data
    Metadata{i} = load(Metadata_Path); %load metadata
end

%% extract distance sensor data
LU = cell(1,N_Sites); %upwind sensor data
LD = cell(1,N_Sites); %downwind sensor data
for i = 1:N_Sites
    LU{i} = Data{i}.ProcessedData.Distance.(UpwindSensor{i});
    LD{i} = Data{i}.ProcessedData.Distance.(DownwindSensor{i});
end

%% plot timeseries
for i = 1:N_Sites
    N_Windows = length(LU{i})
    for j = 1:N_Windows
        zU = mean(LU{i}(j).z.int)-LU{i}(j).z.int;
        zD = mean(LD{i}(j).z.int)-LD{i}(j).z.int;
        tU = LU{i}(j).t.int;
        tD = LD{i}(j).t.int;
        figure;
        plot(tU,zU,tD,zD);
        legend('Upwind','Downwind','Location','SouthEast');
        ylabel('z (mm)');
        title([Sites{i},', ',datestr(LU{i}(j).Date)]);
        set(gca,'FontSize',16);
    end
end