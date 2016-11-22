%% initialize
clearvars;

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/Thresholds/'; %folder for outputs
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for outputs
PrimaryTable_Path = strcat(folder_AnalysisData,'PrimaryTable.csv'); %path for saving primary table
AveragingIntervalTable_Path = strcat(folder_AnalysisData,'AveragingIntervalTable.csv'); %path for averaging interval sensitivity analysis
AnalysisIntervalTable_Path = strcat(folder_AnalysisData,'AnalysisIntervalTable.csv'); %path for analysis interval sensitivity analysis
DiurnalTable_Path = strcat(folder_AnalysisData,'DiurnalTable.csv'); %path for saving diurnal interval sensitivity analysis
DateTable_Path = strcat(folder_AnalysisData,'DateTable.csv'); %path for saving date range sensitivity analysis

%% load data
load(strcat(folder_AnalysisData,'ThresholdAnalysisData')); %load analysis windows
load(strcat(folder_GrainSizeData,'MeanGrainSize')); %load grain size data


%% Primary table
N_Sites = length(Sites);

%generate list of variables
PrimaryTableVars = cell(1,6);
PrimaryTableVars{1,1} = 'Site';
PrimaryTableVars{1,2} = 'd50_surface';
PrimaryTableVars{1,3} = 'tau_ft';
PrimaryTableVars{1,4} = 'tau_it';
PrimaryTableVars{1,5} = 'ustitft_ratio';
PrimaryTableVars{1,6} = 'tauth_flux';

%insert table entries
PrimaryTable = cell(N_Sites,6);
for i = 1:N_Sites
    PrimaryTable{i,1} = Sites{i};
    PrimaryTable{i,2} = [num2str(d50_surface_site(i),'%.3f'), '+/-',num2str(sigma_d50_surface_site(i),'%.3f')];
    PrimaryTable{i,3} = [num2str(tauft_all(i),'%.3f'), '+/-',num2str(sigma_tauft_all(i),'%.3f')];
    PrimaryTable{i,4} = [num2str(tauit_all(i),'%.3f'), '+/-',num2str(sigma_tauit_all(i),'%.3f')];
    PrimaryTable{i,5} = [num2str(ustitftratio_all(i),'%.3f'), '+/-',num2str(sigma_ustitftratio_all(i),'%.3f')];
    PrimaryTable{i,6} = [num2str(tauit_intercept(i),'%.3f'), '+/-',num2str(sigma_tauit_intercept(i),'%.3f')];
end

%convert to table and write
PrimaryTable = cell2table(PrimaryTable,'VariableNames',PrimaryTableVars);
writetable(PrimaryTable,PrimaryTable_Path);


%% AveragingInterval table
N_AveragingInterval = length(deltat_all);

%generate list of variables
AveragingIntervalTableVars = cell(1,4);
AveragingIntervalTableVars{1,1} = 'delta_t';
AveragingIntervalTableVars{1,2} = 'tau_ft';
AveragingIntervalTableVars{1,3} = 'tau_it';
AveragingIntervalTableVars{1,4} = 'ustitft_ratio';

%insert table entries
AveragingIntervalTable = cell(N_AveragingInterval,4);
for i = 1:N_AveragingInterval
    AveragingIntervalTable{i,1} = num2str(seconds(deltat_all(i)),'%.2f');
    AveragingIntervalTable{i,2} = [num2str(tauft_deltat(i),'%.3f'), '+/-',num2str(sigma_tauft_deltat(i),'%.3f')];
    AveragingIntervalTable{i,3} = [num2str(tauit_deltat(i),'%.3f'), '+/-',num2str(sigma_tauit_deltat(i),'%.3f')];
    AveragingIntervalTable{i,4} = [num2str(ustitftratio_deltat(i),'%.3f'), '+/-',num2str(sigma_ustitftratio_deltat(i),'%.3f')];
end

%convert to table and write
AveragingIntervalTable = cell2table(AveragingIntervalTable,'VariableNames',AveragingIntervalTableVars);
writetable(AveragingIntervalTable,AveragingIntervalTable_Path);


%% AnalysisInterval table
N_AnalysisInterval = length(Deltat_all);

%generate list of variables
AnalysisIntervalTableVars = cell(1,4);
AnalysisIntervalTableVars{1,1} = 'Delta_t';
AnalysisIntervalTableVars{1,2} = 'tau_ft';
AnalysisIntervalTableVars{1,3} = 'tau_it';
AnalysisIntervalTableVars{1,4} = 'ustitft_ratio';

%insert table entries
AnalysisIntervalTable = cell(N_AnalysisInterval,4);
for i = 1:N_AnalysisInterval
    AnalysisIntervalTable{i,1} = num2str(seconds(Deltat_all(i)),'%.2f');
    AnalysisIntervalTable{i,2} = [num2str(tauft_Deltat(i),'%.3f'), '+/-',num2str(sigma_tauft_Deltat(i),'%.3f')];
    AnalysisIntervalTable{i,3} = [num2str(tauit_Deltat(i),'%.3f'), '+/-',num2str(sigma_tauit_Deltat(i),'%.3f')];
    AnalysisIntervalTable{i,4} = [num2str(ustitftratio_Deltat(i),'%.3f'), '+/-',num2str(sigma_ustitftratio_Deltat(i),'%.3f')];
end

%convert to table and write
AnalysisIntervalTable = cell2table(AnalysisIntervalTable,'VariableNames',AnalysisIntervalTableVars);
writetable(AnalysisIntervalTable,AnalysisIntervalTable_Path);

%% Diurnal table
N_Diurnal = length(diurnalrange_starthour);

%generate list of variables
DiurnalTableVars = cell(1,4);
DiurnalTableVars{1,1} = 'times';
DiurnalTableVars{1,2} = 'tau_ft';
DiurnalTableVars{1,3} = 'tau_it';
DiurnalTableVars{1,4} = 'ustitft_ratio';

%insert table entries
DiurnalTable = cell(N_Diurnal,4);
for i = 1:N_Diurnal
    DiurnalTable{i,1} = [int2str(diurnalrange_starthour(i)),'-',int2str(diurnalrange_endhour(i))];
    DiurnalTable{i,2} = [num2str(tauft_diurnalrange(i),'%.3f'), '+/-',num2str(sigma_tauft_diurnalrange(i),'%.3f')];
    DiurnalTable{i,3} = [num2str(tauit_diurnalrange(i),'%.3f'), '+/-',num2str(sigma_tauit_diurnalrange(i),'%.3f')];
    DiurnalTable{i,4} = [num2str(ustitftratio_diurnalrange(i),'%.3f'), '+/-',num2str(sigma_ustitftratio_diurnalrange(i),'%.3f')];
end

%convert to table and write
DiurnalTable = cell2table(DiurnalTable,'VariableNames',DiurnalTableVars);
writetable(DiurnalTable,DiurnalTable_Path);


%% Date table
N_Dates = length(daterange_startdate);

%generate list of variables
DateTableVars = cell(1,5);
DateTableVars{1,1} = 'dates';
DateTableVars{1,2} = 'd_50';
DateTableVars{1,3} = 'tau_ft';
DateTableVars{1,4} = 'tau_it';
DateTableVars{1,5} = 'ustitft_ratio';

%insert table entries
DateTable = cell(N_Dates,5);
for i = 1:N_Dates
    DateTable{i,1} = [datestr(daterange_startdate(i),'mmm dd'),'-',datestr(daterange_enddate(i),'DD')];
    DateTable{i,2} = [num2str(d50_surface_cluster(i),'%.3f'), '+/-',num2str(sigma_d50_surface_cluster(i),'%.3f')];
    DateTable{i,3} = [num2str(tauft_daterange(i),'%.3f'), '+/-',num2str(sigma_tauft_daterange(i),'%.3f')];
    DateTable{i,4} = [num2str(tauit_daterange(i),'%.3f'), '+/-',num2str(sigma_tauit_daterange(i),'%.3f')];
    DateTable{i,5} = [num2str(ustitftratio_daterange(i),'%.3f'), '+/-',num2str(sigma_ustitftratio_daterange(i),'%.3f')];
end

%convert to table and write
DateTable = cell2table(DateTable,'VariableNames',DateTableVars);
writetable(DateTable,DateTable_Path);