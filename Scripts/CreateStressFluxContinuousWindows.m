%% KEEP ONLY TIME INTERVALS OF CONTINUOUS TRANSPORT

%initialize
clearvars;

%information about where to load data and save data
folder_AnalysisData = '../AnalysisData/'; %folder for storing data output

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows_all'));

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
u_bar_raw_continuous = cell(N_Sites,1); %raw mean u values
u_std_raw_continuous = cell(N_Sites,1); %standard deviation u values
u2_raw_continuous = cell(N_Sites,1); %mean u^2
u2_ex_raw_continuous = cell(N_Sites,1); %mean u^2 - u^2th
eta_raw_continuous = cell(N_Sites,1); %mean u2ex/mean u2
u_bar_cal_continuous = cell(N_Sites,1); %raw mean u values
u_std_cal_continuous = cell(N_Sites,1); %standard deviation u values
u2_cal_continuous = cell(N_Sites,1); %mean u^2
u2_ex_cal_continuous = cell(N_Sites,1); %mean u^2 - u^2th
eta_cal_continuous = cell(N_Sites,1); %mean u2ex/mean u2
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
u_raw_profile_continuous = cell(N_Sites,1); %raw velocity profile
u_cal_profile_continuous = cell(N_Sites,1); %calibrated velocity profile
z_profile_continuous = cell(N_Sites,1); %heights for profile
Anemometer_profile_continuous = cell(N_Sites,1); %names of anemometers for profile
u_raw_lowest_continuous = cell(N_Sites,1); %raw velocity profile - lowest anemometers
u_cal_lowest_continuous = cell(N_Sites,1); %calibrated velocity profile - lowest anemometers
z_lowest_continuous = cell(N_Sites,1); %heights of lowest anemometers
Anemometer_lowest_continuous = cell(N_Sites,1); %names of lowest anemometers for profile
ustLog_raw_continuous = cell(N_Sites,1); %u* calculated from mean of raw velocities for lowest anemometers in profile
z0_raw_continuous = cell(N_Sites,1); %z0 calculated from mean of raw velocities for lowest anemometers in profile
ustLog_cal_continuous = cell(N_Sites,1); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
z0_cal_continuous = cell(N_Sites,1); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile

for i=1:N_Sites
    %ind_continuous = intersect(find(fQ_all{i}>=fQ_threshold),find(Q_all{i}>0));
    ind_continuous = find(fQ_all{i}>=fQ_threshold);
    Q_continuous{i} = Q_all{i}(ind_continuous);
    ze_continuous{i} = ze_all{i}(ind_continuous);
    fQ_continuous{i} = fQ_all{i}(ind_continuous);
    q_continuous{i} = q_all{i}(ind_continuous);
    zq_continuous{i} = zq_all{i}(ind_continuous);
    fq_continuous{i} = fq_all{i}(ind_continuous);
    qcal_continuous{i} = qcal_all{i}(ind_continuous);
    u_bar_raw_continuous{i} = u_bar_raw_all{i}(ind_continuous);
    u_std_raw_continuous{i} = u_std_raw_all{i}(ind_continuous);
    u2_raw_continuous{i} = u2_raw_all{i}(ind_continuous);
    u2_ex_raw_continuous{i} = u2_ex_raw_all{i}(ind_continuous);
    eta_raw_continuous{i} = eta_raw_all{i}(ind_continuous);
    u_bar_cal_continuous{i} = u_bar_cal_all{i}(ind_continuous);
    u_std_cal_continuous{i} = u_std_cal_all{i}(ind_continuous);
    u2_cal_continuous{i} = u2_cal_all{i}(ind_continuous);
    u2_ex_cal_continuous{i} = u2_ex_cal_all{i}(ind_continuous);
    eta_cal_continuous{i} = eta_cal_all{i}(ind_continuous);
    ustRe_raw_continuous{i} = ustRe_raw_all{i}(ind_continuous);
    tauRe_raw_continuous{i} = tauRe_raw_all{i}(ind_continuous);
    ustRe_cal_continuous{i} = ustRe_cal_all{i}(ind_continuous);
    tauRe_cal_continuous{i} = tauRe_cal_all{i}(ind_continuous);
    d10_continuous{i} = d10_all{i}(ind_continuous);
    d50_continuous{i} = d50_all{i}(ind_continuous);
    d90_continuous{i} = d90_all{i}(ind_continuous);
    date_continuous{i} = date_all{i}(ind_continuous);
    StartTime_continuous{i} = StartTime_all{i}(ind_continuous);
    EndTime_continuous{i} = EndTime_all{i}(ind_continuous);
    u_raw_profile_continuous{i} = u_raw_profile_all{i}(ind_continuous); %raw velocity profile
    u_cal_profile_continuous{i} = u_cal_profile_all{i}(ind_continuous); %calibrated velocity profile
    z_profile_continuous{i} = z_profile_all{i}(ind_continuous); %heights for profile
    Anemometer_profile_continuous{i} = Anemometer_profile_all{i}(ind_continuous); %names of anemometers for profile
    u_raw_lowest_continuous{i} = u_raw_lowest_all{i}(ind_continuous); %raw velocity profile - lowest anemometers
    u_cal_lowest_continuous{i} = u_cal_lowest_all{i}(ind_continuous); %calibrated velocity profile - lowest anemometers
    z_lowest_continuous{i} = z_lowest_all{i}(ind_continuous); %heights of lowest anemometers
    Anemometer_lowest_continuous{i} = Anemometer_lowest_all{i}(ind_continuous); %names of lowest anemometers for profile
    ustLog_raw_continuous{i} = ustLog_raw_all{i}(ind_continuous); %u* calculated from mean of raw velocities for lowest anemometers in profile
    z0_raw_continuous{i} = z0_raw_all{i}(ind_continuous); %z0 calculated from mean of raw velocities for lowest anemometers in profile
    ustLog_cal_continuous{i} = ustLog_cal_all{i}(ind_continuous); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
    z0_cal_continuous{i} = z0_cal_all{i}(ind_continuous); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile
end

savepath = strcat(folder_AnalysisData,'StressFluxWindows_continuous');
save(savepath,'Sites','*continuous');