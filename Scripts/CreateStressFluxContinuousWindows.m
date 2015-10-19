%% ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output

%load flux stress window data
load(strcat(folder_ProcessedData,'StressFluxWindows'));

%get info about number of sites
N_Sites = length(Sites);

%% create subset lists for fQ >= threshold
fQ_threshold = 1; %transport frequency to include data in analysis
Q_continuous = cell(N_Sites,1);
ze_continuous = cell(N_Sites,1);
fQ_continuous = cell(N_Sites,1);
q_continuous = cell(N_Sites,1);
zq_continuous = cell(N_Sites,1);
fq_continuous = cell(N_Sites,1);
qcal_continuous = cell(N_Sites,1);
ustRe_raw_continuous = cell(N_Sites,1);
tauRe_raw_continuous = cell(N_Sites,1);
ustRe_cal_continuous = cell(N_Sites,1);
tauRe_cal_continuous = cell(N_Sites,1);
d10_continuous = cell(N_Sites,1);
d50_continuous = cell(N_Sites,1);
d90_continuous = cell(N_Sites,1);
date_continuous = cell(N_Sites,1);
StartTime_continuous = cell(N_Sites,1);
EndTime_continuous = cell(N_Sites,1);

for i=1:N_Sites
    ind_continuous = intersect(find(fQ_list{i}>=fQ_threshold),find(Q_list{i}>0));
    Q_continuous{i} = Q_list{i}(ind_continuous);
    ze_continuous{i} = ze_list{i}(ind_continuous);
    fQ_continuous{i} = fQ_list{i}(ind_continuous);
    q_continuous{i} = q_list{i}{ind_continuous};
    zq_continuous{i} = zq_list{i}{ind_continuous};
    fq_continuous{i} = fq_list{i}{ind_continuous};
    qcal_continuous{i} = qcal_list{i}{ind_continuous};
    ustRe_raw_continuous{i} = ustRe_raw_list{i}(ind_continuous);
    tauRe_raw_continuous{i} = tauRe_raw_list{i}(ind_continuous);
    ustRe_cal_continuous{i} = ustRe_cal_list{i}(ind_continuous);
    tauRe_cal_continuous{i} = tauRe_cal_list{i}(ind_continuous);
    d10_continuous{i} = d10_list{i}(ind_continuous);
    d50_continuous{i} = d50_list{i}(ind_continuous);
    d90_continuous{i} = d90_list{i}(ind_continuous);
    date_continuous{i} = date_list{i}(ind_continuous);
    StartTime_continuous{i} = StartTime_list{i}(ind_continuous);
    EndTime_continuous{i} = EndTime_list{i}(ind_continuous);
end

savepath = strcat(folder_ProcessedData,'StressFluxContinuousWindows');
save(savepath,'Sites','*continuous');