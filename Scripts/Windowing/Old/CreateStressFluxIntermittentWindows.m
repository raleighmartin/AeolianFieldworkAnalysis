%% KEEP ONLY TIME INTERVALS OF INTERMITTENT TRANSPORT

%initialize
clearvars;

%information about where to load data and save data
folder_AnalysisData = '../AnalysisData/'; %folder for storing data output

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows_all'));

%get info about number of sites
N_Sites = length(Sites);

%% create subset lists for fQ < threshold
fQ_threshold = 1; %max transport frequency to include data in analysis
Q_intermittent = cell(N_Sites,1);
ze_intermittent = cell(N_Sites,1);
fQ_intermittent = cell(N_Sites,1);
q_intermittent = cell(N_Sites,1);
zq_intermittent = cell(N_Sites,1);
fq_intermittent = cell(N_Sites,1);
qcal_intermittent = cell(N_Sites,1);
u_bar_raw_intermittent = cell(N_Sites,1); %raw mean u values
u_std_raw_intermittent = cell(N_Sites,1); %standard deviation u values
u2_raw_intermittent = cell(N_Sites,1); %mean u^2
u2_ex_raw_intermittent = cell(N_Sites,1); %mean u^2 - u^2th
eta_raw_intermittent = cell(N_Sites,1); %mean u2ex/mean u2
u_bar_cal_intermittent = cell(N_Sites,1); %raw mean u values
u_std_cal_intermittent = cell(N_Sites,1); %standard deviation u values
u2_cal_intermittent = cell(N_Sites,1); %mean u^2
u2_ex_cal_intermittent = cell(N_Sites,1); %mean u^2 - u^2th
eta_cal_intermittent = cell(N_Sites,1); %mean u2ex/mean u2
ustRe_raw_intermittent = cell(N_Sites,1);
tauRe_raw_intermittent = cell(N_Sites,1);
ustRe_cal_intermittent = cell(N_Sites,1);
tauRe_cal_intermittent = cell(N_Sites,1);
d10_intermittent = cell(N_Sites,1);
d50_intermittent = cell(N_Sites,1);
d90_intermittent = cell(N_Sites,1);
date_intermittent = cell(N_Sites,1);
StartTime_intermittent = cell(N_Sites,1);
EndTime_intermittent = cell(N_Sites,1);
u_raw_profile_intermittent = cell(N_Sites,1); %raw velocity profile
u_cal_profile_intermittent = cell(N_Sites,1); %calibrated velocity profile
z_profile_intermittent = cell(N_Sites,1); %heights for profile
Anemometer_profile_intermittent = cell(N_Sites,1); %names of anemometers for profile
u_raw_lowest_intermittent = cell(N_Sites,1); %raw velocity profile - lowest anemometers
u_cal_lowest_intermittent = cell(N_Sites,1); %calibrated velocity profile - lowest anemometers
z_lowest_intermittent = cell(N_Sites,1); %heights of lowest anemometers
Anemometer_lowest_intermittent = cell(N_Sites,1); %names of lowest anemometers for profile
ustLog_raw_intermittent = cell(N_Sites,1); %u* calculated from mean of raw velocities for lowest anemometers in profile
z0_raw_intermittent = cell(N_Sites,1); %z0 calculated from mean of raw velocities for lowest anemometers in profile
ustLog_cal_intermittent = cell(N_Sites,1); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
z0_cal_intermittent = cell(N_Sites,1); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile

for i=1:N_Sites
    %ind_intermittent = intersect(find(fQ_all{i}<fQ_threshold),find(Q_all{i}>0));
    ind_intermittent = find(fQ_all{i}<fQ_threshold);
    Q_intermittent{i} = Q_all{i}(ind_intermittent);
    ze_intermittent{i} = ze_all{i}(ind_intermittent);
    fQ_intermittent{i} = fQ_all{i}(ind_intermittent);
    q_intermittent{i} = q_all{i}(ind_intermittent);
    zq_intermittent{i} = zq_all{i}(ind_intermittent);
    fq_intermittent{i} = fq_all{i}(ind_intermittent);
    qcal_intermittent{i} = qcal_all{i}(ind_intermittent);
    u_bar_raw_intermittent{i} = u_bar_raw_all{i}(ind_intermittent);
    u_std_raw_intermittent{i} = u_std_raw_all{i}(ind_intermittent);
    u2_raw_intermittent{i} = u2_raw_all{i}(ind_intermittent);
    u2_ex_raw_intermittent{i} = u2_ex_raw_all{i}(ind_intermittent);
    eta_raw_intermittent{i} = eta_raw_all{i}(ind_intermittent);
    u_bar_cal_intermittent{i} = u_bar_cal_all{i}(ind_intermittent);
    u_std_cal_intermittent{i} = u_std_cal_all{i}(ind_intermittent);
    u2_cal_intermittent{i} = u2_cal_all{i}(ind_intermittent);
    u2_ex_cal_intermittent{i} = u2_ex_cal_all{i}(ind_intermittent);
    eta_cal_intermittent{i} = eta_cal_all{i}(ind_intermittent);
    ustRe_raw_intermittent{i} = ustRe_raw_all{i}(ind_intermittent);
    tauRe_raw_intermittent{i} = tauRe_raw_all{i}(ind_intermittent);
    ustRe_cal_intermittent{i} = ustRe_cal_all{i}(ind_intermittent);
    tauRe_cal_intermittent{i} = tauRe_cal_all{i}(ind_intermittent);
    d10_intermittent{i} = d10_all{i}(ind_intermittent);
    d50_intermittent{i} = d50_all{i}(ind_intermittent);
    d90_intermittent{i} = d90_all{i}(ind_intermittent);
    date_intermittent{i} = date_all{i}(ind_intermittent);
    StartTime_intermittent{i} = StartTime_all{i}(ind_intermittent);
    EndTime_intermittent{i} = EndTime_all{i}(ind_intermittent);
    u_raw_profile_intermittent{i} = u_raw_profile_all{i}(ind_intermittent); %raw velocity profile
    u_cal_profile_intermittent{i} = u_cal_profile_all{i}(ind_intermittent); %calibrated velocity profile
    z_profile_intermittent{i} = z_profile_all{i}(ind_intermittent); %heights for profile
    Anemometer_profile_intermittent{i} = Anemometer_profile_all{i}(ind_intermittent); %names of anemometers for profile
    u_raw_lowest_intermittent{i} = u_raw_lowest_all{i}(ind_intermittent); %raw velocity profile - lowest anemometers
    u_cal_lowest_intermittent{i} = u_cal_lowest_all{i}(ind_intermittent); %calibrated velocity profile - lowest anemometers
    z_lowest_intermittent{i} = z_lowest_all{i}(ind_intermittent); %heights of lowest anemometers
    Anemometer_lowest_intermittent{i} = Anemometer_lowest_all{i}(ind_intermittent); %names of lowest anemometers for profile
    ustLog_raw_intermittent{i} = ustLog_raw_all{i}(ind_intermittent); %u* calculated from mean of raw velocities for lowest anemometers in profile
    z0_raw_intermittent{i} = z0_raw_all{i}(ind_intermittent); %z0 calculated from mean of raw velocities for lowest anemometers in profile
    ustLog_cal_intermittent{i} = ustLog_cal_all{i}(ind_intermittent); %u* calculated from mean of calibrated velocities for lowest anemometers in profile
    z0_cal_intermittent{i} = z0_cal_all{i}(ind_intermittent); %z0 calculated from mean of calibrated velocities for lowest anemometers in profile
end

savepath = strcat(folder_AnalysisData,'StressFluxWindows_intermittent');
save(savepath,'Sites','*intermittent');