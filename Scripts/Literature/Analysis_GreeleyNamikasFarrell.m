%% SCRIPT TO ANALYZE GREELEY (1996), NAMIKAS (2003), AND FARRELL (2012) FLUX VALUES FROM TEXT FILE
clearvars; %initialize

%% assumptions
ust_relerr = 0.1; %10% relative error for u*
zq_Greeley96_assumed = 0.05; %assumed saltation e-folding height (m), from flux law paper
zq_Namikas03_assumed = 0.049; %assumed saltation e-folding height (m), from flux law paper
zq_Namikas06_assumed = 0.049; %assumed saltation e-folding height (m), from Namikas 03
zq_Farrell12_assumed = 0.081; %assumed saltation e-folding height (m), from flux law paper

%% ust ranges for binning
ust_lower_Farrell12 = [0.41; 0.47; 0.50];
ust_upper_Farrell12 = [0.45; 0.49; 0.53];
N_ustbins_Farrell12 = 3;

%% filepath information
folder_LoadData = '../../AnalysisData/Literature/'; %folder for loading inputs of this analysis
folder_SaveData = '../../AnalysisData/Literature/'; %folder for storing outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions
LoadData_Path = strcat(folder_LoadData,'LitData'); %path for loading data
SaveData_Path = strcat(folder_SaveData,'LitAnalysis'); %path for saving data

%% load data
load(LoadData_Path);

%% Calculate saltation flux and height - Greeley

% Make calculations
zq_Greeley96 = zeros(N_Greeley96,1); %(m)
sigma_zq_Greeley96 = zeros(N_Greeley96,1); %(m)
Q_fit_Greeley96 = zeros(N_Greeley96,1); %(g/m/s)
sigma_Q_fit_Greeley96 = zeros(N_Greeley96,1); %(g/m/s)
for i=1:N_Greeley96
    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q] = qz_profilefit_exponential(q_Greeley96{i}, z_Greeley96{i}, sigma_q_Greeley96{i}, sigma_z_Greeley96{i}, zq_Greeley96_assumed); %fit will assume equal uncertainties
    zq_Greeley96(i) = zq;
    sigma_zq_Greeley96(i) = sigma_zq; %flux height uncertainty
    Q_fit_Greeley96(i) = Q;
    sigma_Q_fit_Greeley96(i) = sigma_Q;
end

% Remove outlier value(s) from data
N_Greeley96 = length(GoodInd_Greeley96);
q_Greeley96 = q_Greeley96(GoodInd_Greeley96);
Q_pub_Greeley96 = Q_pub_Greeley96(GoodInd_Greeley96);
Q_fit_Greeley96 = Q_fit_Greeley96(GoodInd_Greeley96);
sigma_Q_fit_Greeley96 = sigma_Q_fit_Greeley96(GoodInd_Greeley96);
RunName_Greeley96 = RunName_Greeley96(GoodInd_Greeley96);
ust_Greeley96 = ust_Greeley96(GoodInd_Greeley96);
z_Greeley96 = z_Greeley96(GoodInd_Greeley96);
zq_Greeley96 = zq_Greeley96(GoodInd_Greeley96);
sigma_zq_Greeley96 = sigma_zq_Greeley96(GoodInd_Greeley96);


%% Calculate saltation flux and height - Namikas
zq_Namikas03 = zeros(N_Namikas03,1); %e-folding height (m)
sigma_zq_Namikas03 = zeros(N_Namikas03,1); %e-folding height (m)
Q_fit_Namikas03 = zeros(N_Namikas03,1); %total flux (g/m/s)
sigma_Q_fit_Namikas03 = zeros(N_Namikas03,1); %total flux (g/m/s)

% Make calculations
for i=1:N_Namikas03
    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q] = qz_profilefit_exponential(q_Namikas03{i}, z_Namikas03{i}, sigma_q_Namikas03{i}, sigma_z_Namikas03{i}, zq_Namikas03_assumed); %fit will assume equal uncertainties
    zq_Namikas03(i) = zq;
    sigma_zq_Namikas03(i) = sigma_zq;
    Q_fit_Namikas03(i) = Q;
    sigma_Q_fit_Namikas03(i) = sigma_Q;
end

%remove outlier values from data
N_Namikas03 = length(GoodInd_Namikas03);
q_Namikas03 = q_Namikas03(GoodInd_Namikas03);
Q_pub_Namikas03 = Q_pub_Namikas03(GoodInd_Namikas03);
Q_fit_Namikas03 = Q_fit_Namikas03(GoodInd_Namikas03);
sigma_Q_fit_Namikas03 = sigma_Q_fit_Namikas03(GoodInd_Namikas03);
RunName_Namikas03 = RunName_Namikas03(GoodInd_Namikas03);
ust_Namikas03 = ust_Namikas03(GoodInd_Namikas03);
z_Namikas03 = z_Namikas03(GoodInd_Namikas03);
zq_Namikas03 = zq_Namikas03(GoodInd_Namikas03);
sigma_zq_Namikas03 = sigma_zq_Namikas03(GoodInd_Namikas03);


%% Calculate size-selective saltation heights - Namikas
zqi_Namikas06 = zeros(N_Namikas06,N_d_Namikas06); %e-folding height (m)
sigma_zqi_Namikas06 = zeros(N_Namikas06,N_d_Namikas06); %e-folding height (m)
Qi_fit_Namikas06 = zeros(N_Namikas03,N_d_Namikas06); %total flux (g/m/s)
sigma_Qi_fit_Namikas06 = zeros(N_Namikas03,N_d_Namikas06); %total flux (g/m/s)
% z_Namikas06 = cell(N_Namikas06,1); %initialize array of calculated trap heights

% Make calculations
for i=1:N_Namikas06
    
%     %initialize matrix of calculated trap heights
%     z_Namikas06{i} = zeros(N_z_Namikas06, N_d_Namikas06);
    
    for j=1:N_d_Namikas06
        %iterative fit to profile to optimize trap heights for fitting
        [z_profile,q0,zq,Q,~,sigma_zq,sigma_Q,~,~,~,~,z_profile_geomean,q0_geomean,zq_geomean,Q_geomean] = ...
            BSNE_profilefit_exponential(qi_Namikas06{i}(:,j), z_bottom_Namikas06, H_Namikas06,...
            sigma_qi_Namikas06{i}(:,j), sigma_z_bottom_Namikas06, zq_Namikas06_assumed);
    
%         %add final list of profile heights to array
%         z_Namikas06{i}(:,j) = z_profile;
        
        %add fluxes to array
        zqi_Namikas06(i,j) = zq; %e-folding height (m)
        sigma_zqi_Namikas06(i,j) = sigma_zq; %e-folding height uncertainty (m)
        Qi_fit_Namikas06(i,j) = Q; %total flux (g/m/s)
        sigma_Qi_fit_Namikas06(i,j) = sigma_Q; %total flux (g/m/s)
    end
end

%% Calculate height-specific PSDs - Namikas
z_Namikas06 = z_profile_calc_exponential(z_bottom_Namikas06,H_Namikas06,zq_Namikas06_assumed); %calculate trap midpoint heights based on assumed zq
znorm_Namikas06 = z_Namikas06/zq_Namikas06_assumed; %calculate z/zq
znorm_bottom_Namikas06 = z_bottom_Namikas06/zq_Namikas06_assumed; %calculate z/zq for bottom of trap
znorm_top_Namikas06 = z_top_Namikas06/zq_Namikas06_assumed; %calculate z/zq for top of trap
dlogd_Namikas06 = log(d_upper_Namikas06)-log(d_lower_Namikas06); %calculate dlogD
dVdlogd_Namikas06 = cell(N_Namikas06,1); %initialize array for size distributions
dbar_airborne_Namikas06 = cell(N_Namikas06,1); %compute mean grain size at each height
d50_airborne_Namikas06 = cell(N_Namikas06,1); %compute median grain size at each height
d90_airborne_Namikas06 = cell(N_Namikas06,1); %compute d90 at each height

% Make calculations
for i=1:N_Namikas06
    dVdlogd_Namikas06{i} = zeros(N_z_Namikas06, N_d_Namikas06); %initialize matrix for size distributions
    dbar_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for mean sizes
    d50_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for median sizes
    d90_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for d90
    for j=1:N_z_Namikas06
        dV = (qi_Namikas06{i}(j,:)/sum(qi_Namikas06{i}(j,:)));
        dVdlogd_Namikas06{i}(j,:) = dV./dlogd_Namikas06;
        dbar_airborne_Namikas06{i}(j) = sum(dV.*d_Namikas06);
        [d10, d50, d90] = ReferenceGrainSizes(dV, d_lower_Namikas06, d_upper_Namikas06);
        d50_airborne_Namikas06{i}(j) = d50;
        d90_airborne_Namikas06{i}(j) = d90;
    end
end

%% Calculate saltation flux and height - Farrell
z_Farrell12 = cell(N_Farrell12,1); %trap height (m)
znorm_Farrell12 = cell(N_Farrell12,1); %normalized trap height (m) - use mean zq for site
zq_Farrell12 = zeros(N_Farrell12,1); %e-folding height (m)
sigma_zq_Farrell12 = zeros(N_Farrell12,1); %e-folding height (m)
Q_fit_Farrell12 = zeros(N_Farrell12,1); %total flux (g/m/s)
sigma_Q_fit_Farrell12 = zeros(N_Farrell12,1); %total flux (g/m/s)

% Allocate imported array to column variable names and cell arrays
for i = 1:N_Farrell12
    
    %iterative fit to profile to optimize trap heights for fitting
    [z_profile,q0,zq,Q,~,sigma_zq,sigma_Q,~,~,~,~,z_profile_geomean,q0_geomean,zq_geomean,Q_geomean] = ...
        BSNE_profilefit_exponential(q_Farrell12{i}, z_bottom_Farrell12{i}, H_Farrell12{i},...
        sigma_q_Farrell12{i}, sigma_z_bottom_Farrell12{i}, zq_Farrell12_assumed);
    
    %add final list of profile heights to array
    z_Farrell12{i} = z_profile;
    znorm_Farrell12{i} = z_profile/zq_Farrell12_assumed;
    
    %add fluxes to array
    zq_Farrell12(i) = zq; %e-folding height (m)
    sigma_zq_Farrell12(i) = sigma_zq; %e-folding height uncertainty (m)
    Q_fit_Farrell12(i) = Q; %total flux (g/m/s)
    sigma_Q_fit_Farrell12(i) = sigma_Q; %total flux (g/m/s)
end

% estimate d90 values from mu and sigma
d90_airborne_Farrell12 = cell(N_Farrell12,1);
for i = 1:N_Farrell12
    d90_airborne_Farrell12{i} = zeros(size(dbar_airborne_Farrell12{i}));
    for j=1:length(dbar_airborne_Farrell12{i})
        d90_log = norminv(0.9,log10(dbar_airborne_Farrell12{i}(j)),log10(dsigma_airborne_Farrell12{i}(j))); %estimate d90 for lognormal distribution
        d90 = 10^d90_log; %convert d90 out of log units
        d90_airborne_Farrell12{i}(j) = d90;
    end
end

%remove outlier values
N_Farrell12 = length(GoodInd_Farrell12);
q_Farrell12 = q_Farrell12(GoodInd_Farrell12);
Q_fit_Farrell12 =Q_fit_Farrell12(GoodInd_Farrell12);
sigma_Q_fit_Farrell12 =sigma_Q_fit_Farrell12(GoodInd_Farrell12);
RunName_Farrell12 = RunName_Farrell12(GoodInd_Farrell12);
ust_Farrell12 = ust_Farrell12(GoodInd_Farrell12);
sigma_ust_Farrell12 = sigma_ust_Farrell12(GoodInd_Farrell12);
z_Farrell12 = z_Farrell12(GoodInd_Farrell12);
znorm_Farrell12 = znorm_Farrell12(GoodInd_Farrell12);
z_bottom_Farrell12 = z_bottom_Farrell12(GoodInd_Farrell12);
H_Farrell12 = H_Farrell12(GoodInd_Farrell12);
zq_Farrell12 = zq_Farrell12(GoodInd_Farrell12);
sigma_zq_Farrell12 = sigma_zq_Farrell12(GoodInd_Farrell12);
dbar_airborne_Farrell12 = dbar_airborne_Farrell12(GoodInd_Farrell12);
d90_airborne_Farrell12 = d90_airborne_Farrell12(GoodInd_Farrell12);

%group dbar and d90 and z by u*
z_airborne_ustbin_Farrell12 = cell(N_ustbins_Farrell12,1);
znorm_airborne_ustbin_Farrell12 = cell(N_ustbins_Farrell12,1);
dbar_airborne_ustbin_Farrell12 = cell(N_ustbins_Farrell12,1);
d90_airborne_ustbin_Farrell12 = cell(N_ustbins_Farrell12,1);
for i = 1:N_ustbins_Farrell12
    ind_ustbin = intersect(find(ust_Farrell12>=ust_lower_Farrell12(i)),...
        find(ust_Farrell12<=ust_upper_Farrell12(i))); %indices of profiles in ust bin
    z_airborne_ustbin_Farrell12{i} = zeros(N_z_Farrell12,1);
    dbar_airborne_ustbin_Farrell12{i} = zeros(N_z_Farrell12,1);
    d90_airborne_ustbin_Farrell12{i} = zeros(N_z_Farrell12,1);
    
    %go through each height for this bin
    for j = 1:N_z_Farrell12
        z_ustbin = []; %intialize list of zs for this bin
        dbar_ustbin = []; %initialize list of dbars for this bin
        d90_ustbin = []; %initialize list of d90s for this bin
        
        %go through each profile for this height
        for k = 1:length(ind_ustbin)
            z_profile = z_Farrell12{ind_ustbin(k)}; %get z profiles for u* bin
            ind_z = intersect(find(z_profile>z_bottom_Farrell12{1}(j)),...
                find(z_profile<z_bottom_Farrell12{1}(j)+H_Farrell12{1}(j))); %get ind of z for height
            z_ustbin = [z_ustbin; z_profile(ind_z)]; %add z to list (if it exists)
            dbar_ustbin = [dbar_ustbin; dbar_airborne_Farrell12{ind_ustbin(k)}(ind_z)]; %add dbar to list (if it exists)
            d90_ustbin = [d90_ustbin; d90_airborne_Farrell12{ind_ustbin(k)}(ind_z)]; %add d90 to list (if it exists)
        end
        z_airborne_ustbin_Farrell12{i}(j) = mean(z_ustbin(~isnan(dbar_ustbin))); %get mean z for ustbin
        dbar_airborne_ustbin_Farrell12{i}(j) = mean(dbar_ustbin(~isnan(dbar_ustbin))); %get mean dbar for ustbin
        d90_airborne_ustbin_Farrell12{i}(j) = mean(d90_ustbin(~isnan(d90_ustbin))); %get mean d90 for ustbin
    end
    znorm_airborne_ustbin_Farrell12{i} = z_airborne_ustbin_Farrell12{i}/zq_Farrell12_assumed;
end

%% save useful data
save(SaveData_Path,'*Greeley96','*Namikas03','*Namikas06','*Farrell12');