%% initialize
clearvars;

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for outputs
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for outputs
folder_DataBSNE = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder containing BSNE data
SiteTable_Path = strcat(folder_AnalysisData,'SiteTable.csv'); %path for saving site table
DailyTable_Path = strcat(folder_AnalysisData,'DailyTable.csv'); %path for saving daily table
FitTable_Path = strcat(folder_AnalysisData,'FitTable.csv'); %path for saving fit table
StressFluxTable_Path = strcat(folder_AnalysisData,'StressFluxTable.csv'); %path for saving stress flux table

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);
N_sigma_tauit = 2; %std deviations from impact threshold for fitting / plotting
fQ_Q_tau = 0.1; %minimum mean fQ for Q versus tau comparison (and determination of threshold)
LitNames = {'Greeley'; 'Namikas';'Farrell'};
d50_Lit = {'0.23','0.25','NaN'};
N_Lit = length(LitNames);

%% plotting information
PlotFont = 12;

%% load analysis windows
load(strcat(folder_AnalysisData,'StressFluxWindows_Analysis'));

%% load grain size data
load(strcat(folder_GrainSizeData,'MeanGrainSize'));

%% load BSNE data
FluxBSNE = cell(N_Sites,1);
for i = 1:N_Sites
    FluxBSNE_Site = load(strcat(folder_DataBSNE,'FluxBSNE_',Sites{i})); 
    FluxBSNE{i} = FluxBSNE_Site.FluxBSNE;
end
    
%% SITE TABLE
SiteTable = cell(6,N_Sites+1);
SiteTable{1,1} = 'd10_surface';
SiteTable{2,1} = 'd50_surface';
SiteTable{3,1} = 'd90_surface';
SiteTable{4,1} = 'd10_airborne';
SiteTable{5,1} = 'd50_airborne';
SiteTable{6,1} = 'd90_airborne';
for i = 1:N_Sites
    SiteTable{1,i+1} = num2str(d10_surface_site(i),'%.2f');
    SiteTable{2,i+1} = num2str(d50_surface_site(i),'%.2f');
    SiteTable{3,i+1} = num2str(d90_surface_site(i),'%.2f');
    SiteTable{4,i+1} = num2str(d10_airborne_site(i),'%.2f');
    SiteTable{5,i+1} = num2str(d50_airborne_site(i),'%.2f');
    SiteTable{6,i+1} = num2str(d90_airborne_site(i),'%.2f');    
end
SiteTableVars = cell(1,4);
SiteTableVars{1} = 'variable';
for i = 1:N_Sites
    SiteTableVars{i+1} = Sites{i};
end

%% convert to table
SiteTable = cell2table(SiteTable,'VariableNames',SiteTableVars);


%% FIT VALUE TABLE
site_fit = cell(N_Sites+N_Lit,1);
d50_fit = cell(N_Sites+N_Lit,1);
zq_a_fit = cell(N_Sites+N_Lit,1);
zq_b_fit = cell(N_Sites+N_Lit,1);
zqnorm_bar_fit = cell(N_Sites+N_Lit,1);
tauit_linearfit = cell(N_Sites+N_Lit,1);
ustit_linearfit = cell(N_Sites+N_Lit,1);
Qhatalt_bar_fit = cell(N_Sites+N_Lit,1);
Qhat_bar_fit = cell(N_Sites+N_Lit,1);
for i = 1:N_Sites
    site_fit{i} = Sites{i};
    d50_fit{i} = num2str(d50_surface_site(i),'%.2f');
    zq_a_fit{i} = [num2str(slope_zqustfit_all(i),'%.3f'),' +/- ',num2str(sigma_slope_zqustfit_all(i),'%.3f')];
    zq_b_fit{i} = [num2str(intercept_zqustfit_all(i),'%.3f'),' +/- ',num2str(sigma_intercept_zqustfit_all(i),'%.3f')];
    zqnorm_bar_fit{i} = [int2str(zqnorm_bar_all(i)),' +/- ',int2str(sigma_zqnorm_bar_all(i))];
    tauit_linearfit{i} = [num2str(tauit_linearfit_all(i),'%.3f'),' +/- ',num2str(tauit_sigma_linearfit_all(i),'%.3f')];
    ustit_linearfit{i} = [num2str(ustit_linearfit_all(i),'%.3f'),' +/- ',num2str(ustit_sigma_linearfit_all(i),'%.3f')];
    Qhatalt_bar_fit{i} = [num2str(Qhatalt_bar_all(i),'%.2f'),' +/- ',num2str(sigma_Qhatalt_bar_all(i),'%.2f')];
    Qhat_bar_fit{i} = [num2str(Qhat_bar_all(i),'%.2f'),' +/- ',num2str(sigma_Qhat_bar_all(i),'%.2f')];
end
for i = 1:N_Lit
    site_fit{i+N_Sites} = LitNames{i};
    d50_fit{i+N_Sites} = d50_Lit{i};
    zq_a_fit{i+N_Sites} = [num2str(slope_zqust_lit_all(i),'%.3f'),' +/- ',num2str(sigma_slope_zqust_lit_all(i),'%.3f')];
    zq_b_fit{i+N_Sites} = [num2str(intercept_zqust_lit_all(i),'%.3f'),' +/- ',num2str(sigma_intercept_zqust_lit_all(i),'%.3f')];
    if i ~= 3
        zqnorm_bar_fit{i+N_Sites} = [int2str(zqnorm_bar_lit_all(i)),' +/- ',int2str(sigma_zqnorm_bar_lit_all(i))];
    elseif i == 3
        zqnorm_bar_fit{i+N_Sites} = 'NaN'; %leave last entry (for Namikas) empty
    end
    tauit_linearfit{i+N_Sites} = 'NaN';
    ustit_linearfit{i+N_Sites} = 'NaN';
    Qhat_bar_fit{i+N_Sites} = 'NaN';
end

%create structured array
FitTable = struct(...
    'site',site_fit,...
    'd50',d50_fit,...
    'a',zq_a_fit,...
    'b',zq_b_fit,...
    'zqnorm',zqnorm_bar_fit,...
    'tauit',tauit_linearfit,...
    'ustit',ustit_linearfit,...
    'Qhatalt',Qhatalt_bar_fit,...
    'Qhat',Qhat_bar_fit);

%% convert to table
FitTable = struct2table(FitTable);


%% STRESS FLUX TABLE
site_fit = cell(N_Sites,1);
d50_fit = cell(N_Sites,1);
C_linearfit = cell(N_Sites,1);
tauit_linearfit = cell(N_Sites,1);
Chi2nu_linearfit = cell(N_Sites,1);
C_nonlinearfit = cell(N_Sites,1);
n_nonlinearfit = cell(N_Sites,1);
tauit_nonlinearfit = cell(N_Sites,1);
Chi2nu_nonlinearfit = cell(N_Sites,1);
for i = 1:N_Sites
    site_fit{i} = Sites{i};
    d50_fit{i} = num2str(d50_surface_site(i),'%.2f');
    C_linearfit{i} = [num2str(C_linearfit_all(i),'%.0f'),' +/- ',num2str(C_sigma_linearfit_all(i),'%.0f')];
    tauit_linearfit{i} = [num2str(tauit_linearfit_all(i),'%.3f'),' +/- ',num2str(tauit_sigma_linearfit_all(i),'%.3f')];
    Chi2nu_linearfit{i} = num2str(Chi2_linearfit_all(i)./df_linearfit_all(i),'%.2f');
    C_nonlinearfit{i} = [num2str(C_nonlinearfit_all(i),'%.0f'),' [',num2str(C_range_nonlinearfit_all(i,1),'%.0f'),' ',num2str(C_range_nonlinearfit_all(i,2),'%.0f'),']'];
    n_nonlinearfit{i} = [num2str(n_nonlinearfit_all(i),'%.2f'),' [',num2str(n_range_nonlinearfit_all(i,1),'%.2f'),' ',num2str(n_range_nonlinearfit_all(i,2),'%.2f'),']'];
    tauit_nonlinearfit{i} = [num2str(tauit_nonlinearfit_all(i),'%.3f'),' [',num2str(tauit_range_nonlinearfit_all(i,1),'%.3f'),' ',num2str(tauit_range_nonlinearfit_all(i,2),'%.3f'),']'];
    Chi2nu_nonlinearfit{i} = num2str(Chi2_nonlinearfit_all(i)./df_nonlinearfit_all(i),'%.2f');
end
    
%create structured array
StressFluxTable = struct(...
    'site',site_fit,...
    'd50',d50_fit,...
    'C_linearfit',C_linearfit,...
    'tauit_linearfit',tauit_linearfit,...
    'Chi2nu_linearfit',Chi2nu_linearfit,...
    'C_nonlinearfit',C_nonlinearfit,...
    'n_nonlinearfit',n_nonlinearfit,... 
    'tauit_nonlinearfit',tauit_nonlinearfit,...
    'Chi2nu_nonlinearfit',Chi2nu_nonlinearfit);

%% convert to table
StressFluxTable = struct2table(StressFluxTable);


%% DAILY TABLE
%% gather together all dates and analysis windows
dates_site = cell(N_Sites,1); %get dates for site
ind_analysis = cell(N_Sites,1); %get indices of windows for analysis at site
%(must either have above threshold stress or nonzero flux)

%go through each site to populate cell arrays
for i = 1:N_Sites
    ind_abovethreshold = find(tauRe_all{i}>tauit_linearfit_all(i)+tauit_sigma_linearfit_all(i)*N_sigma_tauit);
    ind_transport = find(Q_all{i}>0);
    ind_analysis{i} = union(ind_abovethreshold,ind_transport);
    dates_site{i} = unique(Date_all{i}(ind_analysis{i}));
end
N_dates_site = cellfun(@length,dates_site); %get number of dates for each site


%% get daily information for site
N_windows_daily = cell(N_Sites,1); %number of time windows per day
zU_daily = cell(N_Sites,1); %get minimum anemometer height by day
zW_min_daily = cell(N_Sites,1); %get minimum Wenglor height by day
zW_max_daily = cell(N_Sites,1); %get maximum Wenglor height by day
Q_min_daily = cell(N_Sites,1); %get minimum flux by day
Q_max_daily = cell(N_Sites,1); %get maximum flux by day
fQ_min_daily = cell(N_Sites,1); %get minimum flux frequency by day
fQ_max_daily = cell(N_Sites,1); %get maximum flux frequency by day
ust_min_daily = cell(N_Sites,1); %get minimum u* by day
ust_max_daily = cell(N_Sites,1); %get maximum u* by day

%go through each site to populate cell arrays
for i = 1:N_Sites
    N_windows_daily{i} = zeros(N_dates_site(i),1);
    zU_daily{i} = zeros(N_dates_site(i),1);
    zW_min_daily{i} = zeros(N_dates_site(i),1);
    zW_max_daily{i} = zeros(N_dates_site(i),1);
    Q_min_daily{i} = zeros(N_dates_site(i),1);
    Q_max_daily{i} = zeros(N_dates_site(i),1);
    fQ_min_daily{i} = zeros(N_dates_site(i),1);
    fQ_max_daily{i} = zeros(N_dates_site(i),1);
    ust_min_daily{i} = zeros(N_dates_site(i),1);
    ust_max_daily{i} = zeros(N_dates_site(i),1);
    for j = 1:N_dates_site(i)
        ind_date = intersect(find(Date_all{i}==dates_site{i}(j)),ind_analysis{i});
       
        N_windows_daily{i}(j) = length(ind_date);

        zU_daily{i}(j) = mode(zU_all{i}(ind_date));

        zW_min_daily{i}(j) = min(zW_min_all{i}(ind_date));
        zW_max_daily{i}(j) = max(zW_max_all{i}(ind_date));
        
        Q_min_daily{i}(j) = min(Q_all{i}(ind_date));
        Q_max_daily{i}(j) = max(Q_all{i}(ind_date));
        
        fQ_min_daily{i}(j) = min(fQ_all{i}(ind_date));
        fQ_max_daily{i}(j) = max(fQ_all{i}(ind_date));
        
        ust_min_daily{i}(j) = min(ustRe_all{i}(ind_date));
        ust_max_daily{i}(j) = max(ustRe_all{i}(ind_date));
    end
end


%% get BSNE heights
zBSNE_min_daily = cell(N_Sites,1); %get minimum BSNE height by day
zBSNE_max_daily = cell(N_Sites,1); %get maximum BSNE height by day

%go through each site to populate cell arrays
for i = 1:N_Sites
    dates_BSNE = [FluxBSNE{i}.Date]';
    zBSNE_min_daily{i} = zeros(N_dates_site(i),1);
    zBSNE_max_daily{i} = zeros(N_dates_site(i),1);

    for j = 1:N_dates_site(i)
        ind_date = find(dates_BSNE == dates_site{i}(j));
        zBSNE_min_profile = zeros(size(ind_date));
        zBSNE_max_profile = zeros(size(ind_date));
        for k = 1:length(ind_date)
            zBSNE_profile = FluxBSNE{i}(ind_date(k)).z.z;
            zBSNE_min_profile(k) = min(zBSNE_profile);
            zBSNE_max_profile(k) = max(zBSNE_profile);
        end
        zBSNE_min_daily{i}(j) = min(zBSNE_min_profile);
        zBSNE_max_daily{i}(j) = max(zBSNE_max_profile);
    end
end


%% get grain sizes
d10_surface_daily = cell(N_Sites,1); %get surface d10s for site
d50_surface_daily = cell(N_Sites,1); %get surface d50s for site
d90_surface_daily = cell(N_Sites,1); %get surface d90s for site
d10_airborne_daily = cell(N_Sites,1); %get airborne d10s for site
d50_airborne_daily = cell(N_Sites,1); %get airborne d50s for site
d90_airborne_daily = cell(N_Sites,1); %get airborne d90s for site

%go through each site to populate cell arrays
for i = 1:N_Sites
    d10_surface_daily{i} = zeros(N_dates_site(i),1)*NaN; %get surface d10s for site
    d50_surface_daily{i} = zeros(N_dates_site(i),1)*NaN; %get surface d50s for site
    d90_surface_daily{i} = zeros(N_dates_site(i),1)*NaN; %get surface d90s for site
    d10_airborne_daily{i} = zeros(N_dates_site(i),1)*NaN; %get airborne d10s for site
    d50_airborne_daily{i} = zeros(N_dates_site(i),1)*NaN; %get airborne d50s for site
    d90_airborne_daily{i} = zeros(N_dates_site(i),1)*NaN; %get airborne d90s for site

    for j = 1:N_dates_site(i)
        ind_date_surface = find(dates_surface{i} == dates_site{i}(j));
        ind_date_airborne = find(dates_airborne{i} == dates_site{i}(j));
        if ~isempty(ind_date_surface)
            d10_surface_daily{i}(j) = d10_surface_date{i}(ind_date_surface);
            d50_surface_daily{i}(j) = d50_surface_date{i}(ind_date_surface);
            d90_surface_daily{i}(j) = d90_surface_date{i}(ind_date_surface);
        end
        if ~isempty(ind_date_airborne)
            d10_airborne_daily{i}(j) = d10_airborne_date{i}(ind_date_airborne);
            d50_airborne_daily{i}(j) = d50_airborne_date{i}(ind_date_airborne);
            d90_airborne_daily{i}(j) = d90_airborne_date{i}(ind_date_airborne);
        end
    end
end


%% CREATE STRUCTURED ARRAY

% create list of sites
N_rows_table = sum(N_dates_site); %get number of rows in table
sites_table = cell(N_rows_table,1); %get list of sites for table
dates_table = cell(N_rows_table,1); %get list of dates for table
N_windows_table = cell(N_rows_table,1); %get list of number of windows for table
Q_table = cell(N_rows_table,1); %get list of fluxes
fQ_table = cell(N_rows_table,1); %get list of flux frequencies
ust_table = cell(N_rows_table,1); %get list of u*
zU_table = cell(N_rows_table,1); %get list of anemometer lower heights for table
zBSNE_table = cell(N_rows_table,1); %get list of BSNE heights
zW_table = cell(N_rows_table,1); %get list of Wenglor heights
d10_surface_table = cell(N_rows_table,1); %get list of surface d10s
d50_surface_table = cell(N_rows_table,1); %get list of surface d10s
d90_surface_table = cell(N_rows_table,1); %get list of surface d10s
d10_airborne_table = cell(N_rows_table,1); %get list of airborne d10s
d50_airborne_table = cell(N_rows_table,1); %get list of airborne d10s
d90_airborne_table = cell(N_rows_table,1); %get list of airborne d10s

ind = 0;
for i = 1:N_Sites
    for j = 1:N_dates_site(i)
        ind = ind+1;
        sites_table{ind} = Sites{i};
        dates_table{ind} = datestr(dates_site{i}(j));
        N_windows_table{ind} = int2str(N_windows_daily{i}(j));
        Q_table{ind} = [num2str(Q_min_daily{i}(j),'%.1f'),'-',num2str(Q_max_daily{i}(j),'%.1f')];
        fQ_table{ind} = [num2str(fQ_min_daily{i}(j),'%.2f'),'-',num2str(fQ_max_daily{i}(j),'%.2f')];
        ust_table{ind} = [num2str(ust_min_daily{i}(j),'%.2f'),'-',num2str(ust_max_daily{i}(j),'%.2f')];
        zU_table{ind} = num2str(round(zU_daily{i}(j),2));
        zBSNE_table{ind} = [num2str(zBSNE_min_daily{i}(j),'%.2f'),'-',num2str(zBSNE_max_daily{i}(j),'%.2f')];
        zW_table{ind} = [num2str(zW_min_daily{i}(j),'%.2f'),'-',num2str(zW_max_daily{i}(j),'%.2f')];
        d10_surface_table{ind} = num2str(d10_surface_daily{i}(j),'%.2f');
        d50_surface_table{ind} = num2str(d50_surface_daily{i}(j),'%.2f');
        d90_surface_table{ind} = num2str(d90_surface_daily{i}(j),'%.2f');
        d10_airborne_table{ind} = num2str(d10_airborne_daily{i}(j),'%.2f');
        d50_airborne_table{ind} = num2str(d50_airborne_daily{i}(j),'%.2f');
        d90_airborne_table{ind} = num2str(d90_airborne_daily{i}(j),'%.2f');
    end
end

%create structured array
DailyTable = struct(...
    'Site',sites_table,...
    'Date',dates_table,...
    'N',N_windows_table,...
    'Q',Q_table,...
    'fQ',fQ_table,...
    'ust',ust_table,...
    'zU',zU_table,...
    'zBSNE',zBSNE_table,...
    'zW',zW_table,...
    'd10_surface',d10_surface_table,...
    'd50_surface',d50_surface_table,...
    'd90_surface',d90_surface_table,...
    'd10_airborne',d10_airborne_table,...
    'd50_airborne',d50_airborne_table,...
    'd90_airborne',d90_airborne_table);

%% convert to table
DailyTable = struct2table(DailyTable);


%% SAVE TABLES
writetable(SiteTable,SiteTable_Path);
writetable(FitTable,FitTable_Path);
writetable(StressFluxTable,StressFluxTable_Path);
writetable(DailyTable,DailyTable_Path);