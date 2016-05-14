%% initialize
clearvars;

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder containing BSNE data
folder_BSNE = '../../AnalysisData/BSNE/'; %folder with BSNE data

%% information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% load BSNE data for each site
FluxBSNE = cell(N_Sites,1);
for i = 1:N_Sites
    FluxBSNE_Site = load(strcat(folder_ProcessedData,'FluxBSNE_',Sites{i})); 
    FluxBSNE{i} = FluxBSNE_Site.FluxBSNE;
end

%% save data
path_SaveData = strcat(folder_BSNE,'FluxBSNE');
save(path_SaveData,'FluxBSNE');