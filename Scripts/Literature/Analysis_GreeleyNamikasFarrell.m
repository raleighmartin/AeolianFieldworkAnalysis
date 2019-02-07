%% SCRIPT TO ANALYZE GREELEY (1996), NAMIKAS (2003), AND FARRELL (2012) FLUX VALUES FROM TEXT FILE
clearvars; %initialize

%% assumptions
ust_relerr = 0.1; %10% relative error for u*
zq_Greeley96_assumed = 0.05; %assumed saltation e-folding height (m), from flux law paper
zq_Namikas03_assumed = 0.049; %assumed saltation e-folding height (m), from flux law paper
zq_Namikas06_assumed = 0.049; %assumed saltation e-folding height (m), from Namikas 03
zq_Farrell12_assumed = 0.081; %assumed saltation e-folding height (m), from flux law paper
rho_Namikas06_assumed = 1.22; %kg/m^3 - use Oceano value
rho_Farrell12_assumed = 1.16; %kg/m^3 - use Jeri value
d50_surface_Namikas06_assumed = 0.40; %mm - use Oceano value
d50_surface_Farrell12_assumed = 0.53; %mm - use Jeri Value
sigma_d50_surface_Namikas06_assumed = 0.04; %mm - use Oceano value
sigma_d50_surface_Farrell12_assumed = 0.07; %mm - use Oceano value
tauit_Namikas06_assumed = 0.094; %Pa - use Oceano value
tauit_Farrell12_assumed = 0.135; %Pa - use Jeri value
sigma_tauit_Namikas06_assumed = 0.006; %Pa - use Oceano value
sigma_tauit_Farrell12_assumed = 0.015; %Pa - use Jeri value

%% ust ranges for binning (Farrell only)
ust_lower_Farrell12 = [0.41; 0.47; 0.50]; %m/s
ust_upper_Farrell12 = [0.45; 0.49; 0.53]; %m/s
N_ustbins_Farrell12 = 3;

%% tau ranges for binning (Farrell only)
tau_lower_Farrell12 = [0.19; 0.25; 0.28]; %Pa
tau_upper_Farrell12 = [0.25; 0.28; 0.33]; %Pa
N_taubins_Farrell12 = 3;

%% taunorm ranges for binning (from size-selective analysis)
taunorm_min_bins = [0.0; 1.0; 1.5; 2.2; 3.0]; %Pa
taunorm_max_bins = [1.0; 1.5; 2.2; 3.0; 4.0]; %Pa
N_taunorm_bins = 5;

%% filepath information
folder_LoadData = '../../AnalysisData/Literature/'; %folder for loading inputs of this analysis
folder_SaveData = '../../AnalysisData/Literature/'; %folder for storing outputs of this analysis
folder_Plots = '../../PlotOutput/SizeSelective/'; %folder for printing plots
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
sigma_ust_Namikas03 = sigma_ust_Namikas03(GoodInd_Namikas03);
z_Namikas03 = z_Namikas03(GoodInd_Namikas03);
zq_Namikas03 = zq_Namikas03(GoodInd_Namikas03);
sigma_zq_Namikas03 = sigma_zq_Namikas03(GoodInd_Namikas03);


%% Calculate size-selective saltation heights - Namikas
zqi_Namikas06 = zeros(N_Namikas06,N_d_Namikas06); %e-folding height (m)
sigma_zqi_Namikas06 = zeros(N_Namikas06,N_d_Namikas06); %e-folding height (m)
Qi_fit_Namikas06 = zeros(N_Namikas06,N_d_Namikas06); %total flux (g/m/s)
sigma_Qi_fit_Namikas06 = zeros(N_Namikas06,N_d_Namikas06); %total flux (g/m/s)

%calculate equivalent tau and taunorm values
tau_Namikas06 = rho_Namikas06_assumed.*ust_Namikas06.^2; %Pa
sigma_tau_Namikas06 = 2*rho_Namikas06_assumed.*ust_Namikas06.*sigma_ust_Namikas06; %Pa
taunorm_Namikas06 = tau_Namikas06/tauit_Namikas06_assumed;
sigma_taunorm_Namikas06 = (1/tauit_Namikas06_assumed)*...
    sqrt(sigma_tau_Namikas06.^2 + sigma_tauit_Namikas06_assumed^2*taunorm_Namikas06.^2); %need to double check this last calc

% Make calculations
for i=1:N_Namikas06
    
    %matrix of calculated trap heights for diagnostics
    z_profile_all = zeros(size(qi_Namikas06{i}));
    
    %list of calculated q0 for diagnostics
    q0_all = zeros(1,N_d_Namikas06);
    
    for j=1:N_d_Namikas06
        ind_fit = find(qi_Namikas06{i}(:,j)>0); %fit only values with qi>0
        qi_fit = qi_Namikas06{i}(ind_fit,j);
        z_bottom_fit = z_bottom_Namikas06(ind_fit);
        H_fit = H_Namikas06(ind_fit);
        sigma_qi_fit = sigma_qi_Namikas06{i}(ind_fit,j);
        sigma_z_bottom_fit = sigma_z_bottom_Namikas06(ind_fit);
        
        %iterative fit to profile to optimize trap heights for fitting
        [z_profile,q0,zq,Q,~,sigma_zq,sigma_Q,~,~,~,~,z_profile_geomean,q0_geomean,zq_geomean,Q_geomean] = ...
            BSNE_profilefit_exponential(qi_fit, z_bottom_fit, H_fit,...
            sigma_qi_fit, sigma_z_bottom_fit, zq_Namikas06_assumed);
        
        %add fluxes to array
        zqi_Namikas06(i,j) = zq; %e-folding height (m)
        sigma_zqi_Namikas06(i,j) = sigma_zq; %e-folding height uncertainty (m)
        Qi_fit_Namikas06(i,j) = Q; %total flux (g/m/s)
        sigma_Qi_fit_Namikas06(i,j) = sigma_Q; %total flux (g/m/s)
        
        %add to matrix of calculated trap heights and for diagnostics
        z_profile_all(ind_fit,j) = z_profile;
        q0_all(j) = q0;
    end
    
    %plot profiles and fits for diagnostic purposes
    figure(i); clf; hold on;
    legend_items = cell(N_d_Namikas06,1);
    for j = 1:N_d_Namikas06 %plot raw profiles
        plot(qi_Namikas06{i}(:,j),z_profile_all(:,j));
        legend_items{j} = ['d = ',num2str(d_Namikas06(j),2),' mm'];
    end
    for j = 1:N_d_Namikas06 %plot fit profiles
        qi_fit = q0_all(j)*exp(-z_profile_all(:,j)/zqi_Namikas06(i,j));
        plot(qi_fit,z_profile_all(:,j),'k--');
    end
    set(gca,'xscale','log','box','on')
    xlabel('$$q_{i}$$ (kg/m$$^2$$/s)','Interpreter','Latex');
    ylabel('$$z$$ (m)','Interpreter','Latex');
    if i==1
        legend(legend_items,'Location','NorthEast');
    else
        legend(legend_items,'Location','SouthWest');
    end
    title(['u_* = ',num2str(ust_Namikas06(i),2),' m/s']);
    print([folder_Plots,'qi_Namikas06_u*_',num2str(ust_Namikas06(i),2),'.png'],'-dpng')
end

%plot zq versus d for diagnostic purposes
figure(4); clf; hold on;
legend_items = cell(N_Namikas06,1);
for i = 1:N_Namikas06
    plot(d_Namikas06,zqi_Namikas06(i,:));
    legend_items{i} = ['u_* = ',num2str(ust_Namikas06(i),2),' m/s'];
end
xlabel('d (mm)');
ylabel('$$z_{q,i}$$ (m)','Interpreter','Latex');
legend(legend_items);
print([folder_Plots,'zqi_di_Namikas06.png'],'-dpng')


%% Calculate height-specific PSDs - Namikas
z_Namikas06 = z_profile_calc_exponential(z_bottom_Namikas06,H_Namikas06,zq_Namikas06_assumed); %calculate trap midpoint heights based on assumed zq
znorm_Namikas06 = z_Namikas06/zq_Namikas06_assumed; %calculate z/zq
znorm_bottom_Namikas06 = z_bottom_Namikas06/zq_Namikas06_assumed; %calculate z/zq for bottom of trap
znorm_top_Namikas06 = z_top_Namikas06/zq_Namikas06_assumed; %calculate z/zq for top of trap
dlogd_profile_Namikas06 = log(d_upper_Namikas06)-log(d_lower_Namikas06); %calculate dlogD
dVdlogd_profile_airborne_Namikas06 = cell(N_Namikas06,1); %initialize array for size distributions
dbar_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute mean grain size at each height
d10_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute d10 grain size at each height
d50_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute median grain size at each height
d90_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute d90 at each height
d10norm_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute normalized d10 grain size at each height
d50norm_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute normalized median grain size at each height
d90norm_profile_airborne_Namikas06 = cell(N_Namikas06,1); %compute normalized d90 at each height

% Make calculations
for i=1:N_Namikas06
    dVdlogd_profile_airborne_Namikas06{i} = zeros(N_z_Namikas06, N_d_Namikas06); %initialize matrix for size distributions
    dbar_profile_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for mean sizes
    d10_profile_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for d10
    d50_profile_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for median sizes
    d90_profile_airborne_Namikas06{i} = zeros(N_z_Namikas06, 1); %initialize vector for d90
    for j=1:N_z_Namikas06
        dV = (qi_Namikas06{i}(j,:)/sum(qi_Namikas06{i}(j,:)));
        dVdlogd_profile_airborne_Namikas06{i}(j,:) = dV./dlogd_profile_Namikas06;
        dbar_profile_airborne_Namikas06{i}(j) = sum(dV.*d_Namikas06);
        [d10, d50, d90] = ReferenceGrainSizes(dV, d_lower_Namikas06, d_upper_Namikas06);
        d10_profile_airborne_Namikas06{i}(j) = d10;
        d50_profile_airborne_Namikas06{i}(j) = d50;
        d90_profile_airborne_Namikas06{i}(j) = d90;
    end
    d10norm_profile_airborne_Namikas06{i} = d10_profile_airborne_Namikas06{i}/d50_surface_Namikas06_assumed; %compute normalized d10 grain size at each height
    d50norm_profile_airborne_Namikas06{i} = d50_profile_airborne_Namikas06{i}/d50_surface_Namikas06_assumed; %compute normalized median grain size at each height
    d90norm_profile_airborne_Namikas06{i} = d90_profile_airborne_Namikas06{i}/d50_surface_Namikas06_assumed; %compute normalized d90 at each height
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

%calculate equivalent tau and taunorm values
tau_Farrell12 = rho_Farrell12_assumed.*ust_Farrell12.^2; %Pa
sigma_tau_Farrell12 = 2*rho_Farrell12_assumed.*ust_Farrell12.*sigma_ust_Farrell12; %Pa
taunorm_Farrell12 = tau_Farrell12/tauit_Farrell12_assumed;
sigma_taunorm_Farrell12 = (1/tauit_Farrell12_assumed)*...
    sqrt(sigma_tau_Farrell12.^2 + sigma_tauit_Farrell12_assumed^2*taunorm_Farrell12.^2);

% estimate d10, d50, d90 values from mu and sigma
d10_profile_airborne_Farrell12 = cell(N_Farrell12,1);
d50_profile_airborne_Farrell12 = cell(N_Farrell12,1);
d90_profile_airborne_Farrell12 = cell(N_Farrell12,1);
for i = 1:N_Farrell12
    d10_profile_airborne_Farrell12{i} = zeros(size(dbar_airborne_Farrell12{i}));
    d50_profile_airborne_Farrell12{i} = zeros(size(dbar_airborne_Farrell12{i}));
    d90_profile_airborne_Farrell12{i} = zeros(size(dbar_airborne_Farrell12{i}));   
    for j=1:length(dbar_airborne_Farrell12{i})
        d10_log = norminv(0.1,log10(dbar_airborne_Farrell12{i}(j)),log10(dsigma_airborne_Farrell12{i}(j))); %estimate d10 for lognormal distribution
        d50_log = norminv(0.5,log10(dbar_airborne_Farrell12{i}(j)),log10(dsigma_airborne_Farrell12{i}(j))); %estimate d50 for lognormal distribution
        d90_log = norminv(0.9,log10(dbar_airborne_Farrell12{i}(j)),log10(dsigma_airborne_Farrell12{i}(j))); %estimate d90 for lognormal distribution
        d10 = 10^d10_log; %convert d10 out of log units
        d50 = 10^d50_log; %convert d50 out of log units
        d90 = 10^d90_log; %convert d90 out of log units
        d10_profile_airborne_Farrell12{i}(j) = d10;
        d50_profile_airborne_Farrell12{i}(j) = d50;
        d90_profile_airborne_Farrell12{i}(j) = d90;
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
tau_Farrell12 = tau_Farrell12(GoodInd_Farrell12);
sigma_tau_Farrell12 = sigma_tau_Farrell12(GoodInd_Farrell12);
taunorm_Farrell12 = taunorm_Farrell12(GoodInd_Farrell12);
sigma_taunorm_Farrell12 = sigma_taunorm_Farrell12(GoodInd_Farrell12);
z_Farrell12 = z_Farrell12(GoodInd_Farrell12);
znorm_Farrell12 = znorm_Farrell12(GoodInd_Farrell12);
z_bottom_Farrell12 = z_bottom_Farrell12(GoodInd_Farrell12);
H_Farrell12 = H_Farrell12(GoodInd_Farrell12);
zq_Farrell12 = zq_Farrell12(GoodInd_Farrell12);
sigma_zq_Farrell12 = sigma_zq_Farrell12(GoodInd_Farrell12);
dbar_airborne_Farrell12 = dbar_airborne_Farrell12(GoodInd_Farrell12);
d10_profile_airborne_Farrell12 = d10_profile_airborne_Farrell12(GoodInd_Farrell12);
d50_profile_airborne_Farrell12 = d50_profile_airborne_Farrell12(GoodInd_Farrell12);
d90_profile_airborne_Farrell12 = d90_profile_airborne_Farrell12(GoodInd_Farrell12);

%group values by taunorm
z_taunorm_Farrell12 = cell(N_taunorm_bins,1);
znorm_taunorm_Farrell12 = cell(N_taunorm_bins,1);
dbar_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
d10_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
d50_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
d90_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
d10norm_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
d50norm_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
d90norm_profile_airborne_taunorm_Farrell12 = cell(N_taunorm_bins,1);
for i = 1:N_taunorm_bins
    ind_taunorm = intersect(find(taunorm_Farrell12>=taunorm_min_bins(i)),...
        find(taunorm_Farrell12<=taunorm_max_bins(i))); %indices of profiles in taunorm bin
    z_taunorm_Farrell12{i} = zeros(N_z_Farrell12,1);
    dbar_profile_airborne_taunorm_Farrell12{i} = zeros(N_z_Farrell12,1);
    d10_profile_airborne_taunorm_Farrell12{i} = zeros(N_z_Farrell12,1);
    d50_profile_airborne_taunorm_Farrell12{i} = zeros(N_z_Farrell12,1);
    d90_profile_airborne_taunorm_Farrell12{i} = zeros(N_z_Farrell12,1);
    
    %go through each height for this bin
    for j = 1:N_z_Farrell12
        z_taunorm = []; %intialize list of zs for this bin
        dbar_taunorm = []; %initialize list of dbars for this bin
        d10_taunorm = []; %initialize list of d10s for this bin
        d50_taunorm = []; %initialize list of d50s for this bin
        d90_taunorm = []; %initialize list of d90s for this bin
        
        %go through each profile for this height
        for k = 1:length(ind_taunorm)
            z_profile = z_Farrell12{ind_taunorm(k)}; %get z profiles for u* bin
            ind_z = intersect(find(z_profile>z_bottom_Farrell12{1}(j)),...
                find(z_profile<z_bottom_Farrell12{1}(j)+H_Farrell12{1}(j))); %get ind of z for height
            z_taunorm = [z_taunorm; z_profile(ind_z)]; %add z to list (if it exists)
            dbar_taunorm = [dbar_taunorm; dbar_airborne_Farrell12{ind_taunorm(k)}(ind_z)]; %add dbar to list (if it exists)
            d10_taunorm = [d10_taunorm; d10_profile_airborne_Farrell12{ind_taunorm(k)}(ind_z)]; %add d10 to list (if it exists)
            d50_taunorm = [d50_taunorm; d50_profile_airborne_Farrell12{ind_taunorm(k)}(ind_z)]; %add d50 to list (if it exists)
            d90_taunorm = [d90_taunorm; d90_profile_airborne_Farrell12{ind_taunorm(k)}(ind_z)]; %add d90 to list (if it exists)
        end
        z_taunorm_Farrell12{i}(j) = mean(z_taunorm(~isnan(dbar_taunorm))); %get mean z for taunormbin
        dbar_profile_airborne_taunorm_Farrell12{i}(j) = mean(dbar_taunorm(~isnan(dbar_taunorm))); %get mean dbar for taunormbin
        d10_profile_airborne_taunorm_Farrell12{i}(j) = mean(d10_taunorm(~isnan(d10_taunorm))); %get mean d10 for taunormbin
        d50_profile_airborne_taunorm_Farrell12{i}(j) = mean(d50_taunorm(~isnan(d50_taunorm))); %get mean d50 for taunormbin
        d90_profile_airborne_taunorm_Farrell12{i}(j) = mean(d90_taunorm(~isnan(d90_taunorm))); %get mean d90 for taunormbin
    end
    znorm_taunorm_Farrell12{i} = z_taunorm_Farrell12{i}/zq_Farrell12_assumed;
    d10norm_profile_airborne_taunorm_Farrell12{i} = d10_profile_airborne_taunorm_Farrell12{i}/d50_surface_Farrell12_assumed;
    d50norm_profile_airborne_taunorm_Farrell12{i} = d50_profile_airborne_taunorm_Farrell12{i}/d50_surface_Farrell12_assumed;
    d90norm_profile_airborne_taunorm_Farrell12{i} = d90_profile_airborne_taunorm_Farrell12{i}/d50_surface_Farrell12_assumed;
end

%% save useful data
save(SaveData_Path,'*Greeley96','*Namikas03','*Namikas06','*Farrell12');