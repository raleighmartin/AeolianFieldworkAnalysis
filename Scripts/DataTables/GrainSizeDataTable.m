%% initialize
clearvars;

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% plotting information
PlotFont = 12;
LineWidth_Surface = 1;

%% load grain size / processed data
FluxBSNE_all = cell(N_Sites,1);
for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path);
    FluxBSNE_all{i} = ProcessedData.FluxBSNE;
    clear ProcessedData;
end

%% get BSNE heights
BSNE_dates = cell(N_Sites,1);
BSNE_z_min = cell(N_Sites,1);
BSNE_z_max = cell(N_Sites,1);

for i = 1:3
    dates_all = [FluxBSNE_all{i}.Date]';
    dates_unique = unique(dates_all);
    N_dates = length(dates_unique);
    BSNE_dates{i} = dates_unique;
    BSNE_z_min{i} = zeros(N_dates,1);
    BSNE_z_max{i} = zeros(N_dates,1);

    for j = 1:N_dates
        BSNE_ind_date = find(dates_all == dates_unique(j));
        BSNE_z_min_date = zeros(size(BSNE_ind_date));
        BSNE_z_max_date = zeros(size(BSNE_ind_date));
        for k = 1:length(BSNE_ind_date)
            BSNE_z_date = FluxBSNE_all{i}(BSNE_ind_date(k)).z.z;
            BSNE_z_min_date(k) = min(BSNE_z_date(BSNE_z_date>0));
            BSNE_z_max_date(k) = max(BSNE_z_date(BSNE_z_date>0));
        end
        BSNE_z_min{i}(j) = min(BSNE_z_min_date);
        BSNE_z_max{i}(j) = max(BSNE_z_max_date);
    end
end