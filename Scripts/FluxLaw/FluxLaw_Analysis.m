%% ANALYZE SALTATION FLUX VERSUS SHEAR STRESS RELATIONSHIP

%%
%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%

%initialize
clearvars;
close all;

%parameters
g = 9.8; %gravity (m/s^2)
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
fQ_Qtau_fit_min = 0.1; %minimum mean fQ for Q versus tau comparison (and determination of threshold)
N_sigma_tauit = 2; %std deviations from impact threshold for minimum tauex

%binning information
bin_N_min = 3; %minimum number of entries for bin
tau_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
fQ_bin_minrange = 0.1; %mininum range of fQ for binning
fQ_bin_maxrange = 0.2; %maximum range of fQ for binning

%% folders for loading and saving data
folder_LoadData = '../../AnalysisData/Windowing/'; %folder for retrieving data for this analysis
folder_SaveData = '../../AnalysisData/FluxLaw/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
folder_GrainSizeData = '../../AnalysisData/BSNE/'; %folder for grain size data
folder_LitData = '../../AnalysisData/Literature/'; %folder for loading lit data

%% paths for loading and saving data - restricted
LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'FluxLawCalcs_30min_Restricted'); %path for 30 minute data

% %% paths for loading and saving data - restricted - alt values
% LoadData_Path = strcat(folder_LoadData,'DataWindowCalcs_30min_Restricted_alt'); %path for 30 minute data
% SaveData_Path = strcat(folder_SaveData,'FluxLawCalcs_30min_Restricted_alt'); %path for 30 minute data

%load data
load(LoadData_Path);
addpath(folder_Functions); %point MATLAB to location of functions

%load grain size data
load(strcat(folder_GrainSizeData,'MeanGrainSize'));
d50_Site = d50_surface_site;
sigma_d50_Site = sigma_d50_surface_site;

%load external data and aggregate it together
load(strcat(folder_LitData,'LitData')); %Literature data
N_Lit = 3; %number of literature sites
LitNames = {'Greeley et al. (1996)'; 'Namikas (2003)';'Farrell et al. (2012)'};
ust_Lit = {ust_Greeley96, ust_Namikas03, ust_Farrell12};
sigma_ust_Lit = {sigma_ust_Greeley96, sigma_ust_Namikas03, sigma_ust_Farrell12};
zq_Lit = {zq_Greeley96, zq_Namikas03, zq_Farrell12};
sigma_zq_Lit = {sigma_zq_Greeley96, sigma_zq_Namikas03, sigma_zq_Farrell12};
zqnorm_Lit = {1e3*zq_Greeley96/d50_Greeley96, 1e3*zq_Namikas03/d50_Namikas03, 1e3*zq_Farrell12/d50_Farrell12};
sigma_zqnorm_Lit = {1e3*sigma_zq_Greeley96/d50_Greeley96, 1e3*sigma_zq_Namikas03/d50_Namikas03, 1e3*sigma_zq_Farrell12/d50_Farrell12};
d50_Lit = [d50_Greeley96 d50_Namikas03 d50_Farrell12];
zqnorm_Lit{3} = zqnorm_Lit{3}*NaN; %change zqnorm_Lit so that Farrell doesn't show up

%% Notes on wind tunnel data for comparison:
% Creyssels et al (2009): d50 = 0.242 mm, zq/d50 = 40 (from exponential fit)
% Ho et al. (2011): d50 = 0.230 mm, zq/d50 = 50 (from exponential fit)
% Ho et al. (2014): d50 = 0.23, 0.63 mm, zq/d50 = 22, 12.5
% (however, these are quoted as being from mean of lognormal distribution, but no equations are given - they refer to appendix but don't actually talk about this in the appendix)

%specify wind tunnel data for comparison
WindTunnelNames = {'Creyssels et al. (2009)';'Ho et al. (2011)';'Ho et al. (2014)'};
d50_WindTunnel = {[0.242];[0.230];[0.230,0.630]};
zqnorm_WindTunnel = {[40];[50];[22,12.5]};
sigma_zqnorm_WindTunnel = {[2];[];[]};


%%
%%%%%%%%%%%%%%%%%%%%%
% REMOVE BAD POINTS %
%%%%%%%%%%%%%%%%%%%%%
variable_list = who('variables','*all');
N_variables = length(variable_list);

% remove points with undefined tau
for i = 1:N_Sites
    ind_tauRe = find(~isnan(tauRe_all{i}));
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_tauRe);']);
    end
end

% remove points with undefined Q
for i = 1:N_Sites
    ind_Q = find(~isnan(Q_all{i}));
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_Q);']);
    end
end

%%
%%%%%%%%%%%%%%%%%%%
% PRIMARY BINNING %
%%%%%%%%%%%%%%%%%%%

%% FOR EACH SITE, PERFORM BINNING BY TAU

%initialize number of elements in bin
bin_N_all = cell(N_Sites,1);

%initialize lists of tau bins by Site
tau_bin_values_all = cell(N_Sites,1); %tau's in each bin
sigma_tau_bin_values_all = cell(N_Sites,1); %sigma_tau's in each bin
tau_bin_avg_all = cell(N_Sites,1); %avg tau in bin
tau_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average tau for bin
tau_bin_std_all = cell(N_Sites,1); %std dev tau in bin
tau_bin_SE_all = cell(N_Sites,1); %SE tau in bin
tau_bin_sigma_all = cell(N_Sites,1); %total estimated tau uncertainty for bin
tau_bin_min_all = cell(N_Sites,1); %get minimum tau for bin

%initialize lists of equivalent ust bins by Site
ust_bin_values_all = cell(N_Sites,1); %ust's in each bin
sigma_ust_bin_values_all = cell(N_Sites,1); %sigma_ust's in each bin
ust_bin_avg_all = cell(N_Sites,1); %avg ust in bin
ust_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average ust for bin
ust_bin_std_all = cell(N_Sites,1); %std dev ust in bin
ust_bin_SE_all = cell(N_Sites,1); %SE ust in bin
ust_bin_sigma_all = cell(N_Sites,1); %total estimated ust uncertainty for bin

%separate fQs into tau bins
fQ_bin_values_all = cell(N_Sites,1); %fQs into tau bins
fQ_bin_avg_all = cell(N_Sites,1); %get average fQ for tau bins
fQ_bin_std_all = cell(N_Sites,1); %get std fQ for tau bins
fQ_bin_stdmin_all = zeros(N_Sites,1); %get mininum std dev for fQ in tau bins
fQ_bin_stdmed_all = zeros(N_Sites,1); %get median std dev for fQ in tau bins
fQ_bin_SE_all = cell(N_Sites,1); %get SE fQ for tau bins
fQ_bin_SEalt_all = cell(N_Sites,1); %get alternate SE fQ for sparse tau bins
fQ_bin_sigma_all = cell(N_Sites,1); %get total uncertainty for fQ in tau bins
fQ_bin_max_all = cell(N_Sites,1); %get maximum fQ in bin

%separate Qs into tau bins
Q_bin_values_all = cell(N_Sites,1); %Q into tau bins
sigma_Q_bin_values_all = cell(N_Sites,1); %sigma_Q into tau bins
Q_bin_avg_all = cell(N_Sites,1); %get average Q for tau bins
Q_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average Q for tau bins
Q_bin_std_all = cell(N_Sites,1); %get std Q for tau bins
Q_bin_stdmin_all = zeros(N_Sites,1); %get mininum std dev for Q in tau bins
Q_bin_stdmed_all = zeros(N_Sites,1); %get median std dev for Q in tau bins
Q_bin_SE_all = cell(N_Sites,1); %get SE Q for tau bins
Q_bin_SEalt_all = cell(N_Sites,1); %get alternate SE Q for sparse tau bins
Q_bin_sigma_all = cell(N_Sites,1); %total estimated Q uncertainty for bin
Q_bin_min_all = cell(N_Sites,1); %get minimum Q in bin
Q_bin_avg_stdmed_all = cell(N_Sites,1); %subset of Q_bin_avg included in stdmed analysis
Q_bin_std_stdmed_all = cell(N_Sites,1); %subset of Q_bin_std included in stdmed analysis

%separate zq's into tau bins
zq_bin_values_all = cell(N_Sites,1); %zq into tau bins
sigma_zq_bin_values_all = cell(N_Sites,1); %sigma_zq into tau bins
zq_bin_avg_all = cell(N_Sites,1); %get average zq for tau bins
zq_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average zq for tau bins
zq_bin_std_all = cell(N_Sites,1); %get std zq for tau bins
zq_bin_stdmin_all = zeros(N_Sites,1); %compute minimum std dev for zq in tau bins
zq_bin_stdmed_all = zeros(N_Sites,1); %compute median std dev for zq in tau bins
zq_bin_SE_all = cell(N_Sites,1); %get SE zq for tau bins
zq_bin_SEalt_all = cell(N_Sites,1); %get alternate SE zq for sparse tau bins
zq_bin_sigma_all = cell(N_Sites,1); %total estimated zq uncertainty for bin
zq_bin_sigmaust_all = cell(N_Sites,1); %additional ust contribution to zq uncertainty
zq_bin_sigmatotal_all = cell(N_Sites,1);  %total uncertainty in zq from zq and ust contributions

%get binned values for all Sites
for i = 1:N_Sites
    
    %get values for binning
    tau = tauRe_all{i};
    sigma_tau = sigma_tauRe_all{i};
    ust = ustRe_all{i};
    sigma_ust = sigma_ustRe_all{i};
    fQ = fQ_all{i};
    Q = Q_all{i};
    sigma_Q = sigma_Q_all{i};
    zq = zq_all{i};
    sigma_zq = sigma_zq_all{i};
    
    %get tau bins
    [tau_bin_values, bin_N, tau_bin_min, tau_bin_max, tau_bin_avg, tau_bin_SE] = ...
        Binning(tau, tau_bin_minrange, tau_bin_maxrange, bin_N_min);
    tau_bin_values_all{i} = tau_bin_values; %taus in each bin
    tau_bin_avg_all{i} = tau_bin_avg; %avg tau in bin
    tau_bin_std_all{i} = tau_bin_SE.*sqrt(bin_N); %std dev of tau in bin
    tau_bin_SE_all{i} = tau_bin_SE; %SE tau in bin
    tau_bin_min_all{i} = tau_bin_min; %min tau in bin
    
    %number of elements in each bin
    bin_N_all{i} = bin_N;

    %get number of bins
    N_bins = length(tau_bin_values);
    
    %initialize number of elements in zq bins (may be less than others because excluding Q=0 values)
    bin_N_zq = zeros(N_bins,1);
    
    %initialize bin values
    sigma_tau_bin_values_all{i} = cell(N_bins,1); %initialize sigma_tau's
    ust_bin_values_all{i} = cell(N_bins,1); %initialize usts
    sigma_ust_bin_values_all{i} = cell(N_bins,1); %initialize sigma_usts
    fQ_bin_values_all{i} = cell(N_bins,1); %initialize fQ's
    Q_bin_values_all{i} = cell(N_bins,1); %initialize Q's
    sigma_Q_bin_values_all{i} = cell(N_bins,1); %initialize sigma_Q's
    Q_bin_min_all{i} = zeros(N_bins,1); %initialize minimum Q's in bins
    fQ_bin_max_all{i} = zeros(N_Sites,1); %initialize maximum Q's in bins
    zq_bin_values_all{i} = cell(N_bins,1); %initialize zq's
    sigma_zq_bin_values_all{i} = cell(N_bins,1); %initialize sigma_zq's
    
    %go through bins to get values
    for j = 1:N_bins
        bin_ind = find(tau>=tau_bin_min(j)&tau<=tau_bin_max(j)); %get indices of values in bin
        
        sigma_tau_bin_values_all{i}{j} = sigma_tau(bin_ind); %sigma_tau's
        
        ust_bin_values_all{i}{j} = ust(bin_ind); %ust's
        sigma_ust_bin_values_all{i}{j} = sigma_ust(bin_ind); %sigma_ust's
        
        fQ_bin_values_all{i}{j} = fQ(bin_ind); %fQ's
        fQ_bin_max_all{i}(j) = max(fQ(bin_ind)); %minimum Q in bin
        
        Q_bin_values_all{i}{j} = Q(bin_ind); %Q's
        sigma_Q_bin_values_all{i}{j} = sigma_Q(bin_ind); %sigma_Q's
        Q_bin_min_all{i}(j) = min(Q(bin_ind)); %minimum Q in bin
        
        bin_ind_zq = intersect(bin_ind,find(Q>0)); %indices for zq also limited by elements with nonzero flux
        bin_N_zq(j) = length(bin_ind_zq); %number of elements in zq bin
        zq_bin_values_all{i}{j} = zq(bin_ind_zq); %zq's
        sigma_zq_bin_values_all{i}{j} = sigma_zq(bin_ind_zq); %sigma_zq's
    end
    
    %get ust avg, std, and SE
    ust_bin_avg_all{i} = cellfun(@mean,ust_bin_values_all{i});
    ust_bin_std_all{i} = cellfun(@std,ust_bin_values_all{i});
    ust_bin_SE_all{i} = ust_bin_std_all{i}./sqrt(bin_N);    

    %get indices for computing mininum standard deviation
    ind_full = find(bin_N>=3); %indices for full bins
    ind_transport = find(fQ_bin_max_all{i}>=fQ_Qtau_fit_min); %indices for bins with transport, as defined by at least one element in bin with fQ>fQ_Qtau_fit_min
    ind_stdmed_Q = intersect(ind_full,ind_transport); %indices for computing median standard deviation
    ind_stdmed_zq = find(bin_N_zq>=3); %indices for computing median standard deviation for zq

    %get fQ avg, std, SE, stdmin, stdmed, and SEalt
    fQ_bin_avg_all{i} = cellfun(@mean,fQ_bin_values_all{i});
    fQ_bin_std_all{i} = cellfun(@std,fQ_bin_values_all{i});
    fQ_bin_SE_all{i} = fQ_bin_std_all{i}./sqrt(bin_N);
    fQ_bin_stdmin_all(i) = min(fQ_bin_std_all{i}(ind_stdmed_Q));
    fQ_bin_stdmed_all(i) = median(fQ_bin_std_all{i}(ind_stdmed_Q));
    fQ_bin_SEalt_all{i} = zeros(N_bins,1);
    fQ_bin_SEalt_all{i}(ind_transport) = fQ_bin_stdmed_all(i)./sqrt(bin_N(ind_transport));
    
    %get Q avg, std, and SE, stdmin, stdmed and SEalt
    Q_bin_avg_all{i} = cellfun(@mean,Q_bin_values_all{i});
    Q_bin_std_all{i} = cellfun(@std,Q_bin_values_all{i});
    Q_bin_SE_all{i} = Q_bin_std_all{i}./sqrt(bin_N);
    Q_bin_stdmin_all(i) = min(Q_bin_std_all{i}(ind_stdmed_Q));
    Q_bin_stdmed_all(i) = median(Q_bin_std_all{i}(ind_stdmed_Q));
    Q_bin_SEalt_all{i} = zeros(N_bins,1);
    Q_bin_SEalt_all{i}(ind_transport) = Q_bin_stdmed_all(i)./sqrt(bin_N(ind_transport));
        
    %revise SEalt to make sure it does not exceed value of Q_bin_avg
    for j = 1:length(ind_transport)
        Q_bin_SEalt_all{i}(ind_transport(j)) = ...
            min([Q_bin_avg_all{i}(ind_transport(j)),Q_bin_SEalt_all{i}(ind_transport(j))]);
    end

    %get subset of Q_bin_avg and Q_bin_std included in stdmed analysis (for later plotting)
    Q_bin_avg_stdmed_all{i} = Q_bin_avg_all{i}(ind_stdmed_Q);
    Q_bin_std_stdmed_all{i} = Q_bin_std_all{i}(ind_stdmed_Q);
    
    %get zQ avg, std, and SE, stdmin, stdmed and SEalt
    zq_bin_avg_all{i} = cellfun(@mean,zq_bin_values_all{i});
    zq_bin_std_all{i} = cellfun(@std,zq_bin_values_all{i});
    zq_bin_SE_all{i} = zq_bin_std_all{i}./sqrt(bin_N);
    zq_bin_stdmin_all(i) = min(zq_bin_std_all{i}(ind_stdmed_zq));
    zq_bin_stdmed_all(i) = median(zq_bin_std_all{i}(ind_stdmed_zq));
    zq_bin_SEalt_all{i} = zeros(N_bins,1);
    zq_bin_SEalt_all{i}(ind_transport) = zq_bin_stdmed_all(i)./sqrt(bin_N_zq(ind_transport));

    %compute sigmaavg (uncertainy on mean of random errors)
    tau_bin_sigmaavg_all{i} = zeros(N_bins,1); %initialize list of uncertainties for mean tau in bin
    ust_bin_sigmaavg_all{i} = zeros(N_bins,1); %initialize list of uncertainties for mean ust in bin
    Q_bin_sigmaavg_all{i} = zeros(N_bins,1); %initialize list of uncertainties for mean Q in bin
    zq_bin_sigmaavg_all{i} = zeros(N_bins,1); %initialize list of uncertainties for mean zq in bin
    for j = 1:N_bins
        [~, tau_bin_sigmaavg] = MeanUncertainty(tau_bin_values_all{i}{j},sigma_tau_bin_values_all{i}{j}); %get uncertainty in avg tau for bin
        tau_bin_sigmaavg_all{i}(j) = tau_bin_sigmaavg;        
        
        [~, ust_bin_sigmaavg] = MeanUncertainty(ust_bin_values_all{i}{j},sigma_ust_bin_values_all{i}{j}); %get uncertainty in avg ust for bin
        ust_bin_sigmaavg_all{i}(j) = ust_bin_sigmaavg;
        
        [~, Q_tau_bin_sigmaavg] = MeanUncertainty(Q_bin_values_all{i}{j},sigma_Q_bin_values_all{i}{j}); %get uncertainty in avg Q for bin
        Q_bin_sigmaavg_all{i}(j) = Q_tau_bin_sigmaavg;
            
        [~, zq_tau_bin_sigmaavg] = MeanUncertainty(zq_bin_values_all{i}{j},sigma_zq_bin_values_all{i}{j}); %get uncertainty in avg zq for bin
        zq_bin_sigmaavg_all{i}(j) = zq_tau_bin_sigmaavg;        
    end
    
    %estimate total sigma for bins
    tau_bin_sigma_all{i} = max([tau_bin_SE_all{i}, tau_bin_sigmaavg_all{i}]')'; %total tau uncertainty (max of SE or uncertainty in mean)
    ust_bin_sigma_all{i} = max([ust_bin_SE_all{i}, ust_bin_sigmaavg_all{i}]')'; %total ust uncertainty (max of SE or uncertainty in mean)
    Q_bin_sigma_all{i} = max([max([Q_bin_SE_all{i},Q_bin_SEalt_all{i}]');...
        Q_bin_sigmaavg_all{i}'])'; %total Q uncertainty (max of SE/SEalt or uncertainty in mean)
    zq_bin_sigma_all{i} = max([max([zq_bin_SE_all{i},zq_bin_SEalt_all{i}]');...
        zq_bin_sigmaavg_all{i}'])'; %total zq uncertainty (max of SE/SEalt or uncertainty in mean)
    fQ_bin_sigma_all{i} = max([fQ_bin_SE_all{i},fQ_bin_SEalt_all{i}]')'; %for fQ, take total uncertainty simply as SE
end

%% COMPUTE ZQNORM
zqnorm_bin_avg_all = cell(N_Sites,1); %initialize average zqnorm for tau bins
zqnorm_bin_sigma_all = cell(N_Sites,1); %intialize uncertainty in zqnorm for tau bins
for i=1:N_Sites
    zqnorm_bin_avg_all{i} = 1000*zq_bin_avg_all{i}./d50_Site(i); % determine zqnorm based on grain size by Site
    zqnorm_bin_sigma_all{i} = (1000./d50_Site(i))*sqrt(zq_bin_sigma_all{i}.^2+(sigma_d50_Site(i)/1000).^2.*zqnorm_bin_avg_all{i}.^2); %calculate uncertainty based on d50 and zq
end

%%
%%%%%%%%%%%
% FITTING %
%%%%%%%%%%%

%% FITTING ZQ VS U*
%initialize values used in zq-ust linear fit
ust_zqustfit_all = cell(N_Sites,1);
sigma_ust_zqustfit_all = cell(N_Sites,1);
zq_zqustfit_all = cell(N_Sites,1);
sigma_zq_zqustfit_all = cell(N_Sites,1);
zqnorm_zqustfit_all = cell(N_Sites,1);
sigma_zqnorm_zqustfit_all = cell(N_Sites,1);
Chi2_zqustfit_all = zeros(N_Sites,1);
df_zqustfit_all = zeros(N_Sites,1);

%initialize values resulting from zq-ust linear fit
intercept_zqustfit_all = zeros(N_Sites,1); %intercept of fit
sigma_intercept_zqustfit_all = zeros(N_Sites,1); %uncertainty in intercept of fit
slope_zqustfit_all = zeros(N_Sites,1); %slope of fit for zq versus ustar
sigma_slope_zqustfit_all = zeros(N_Sites,1); %uncertainty in slope for zq versus ustar
zq_pred_zqustfit_all = cell(N_Sites,1); %predicted zq

%initialize mean values for zq
zq_bar_all = zeros(N_Sites,1); %mean of zq
sigma_zq_bar_all = zeros(N_Sites,1); %std dev of zq

%initialize mean values for zqnorm
zqnorm_bar_all = zeros(N_Sites,1); %mean of zqnorm
sigma_zqnorm_bar_all = zeros(N_Sites,1); %std dev of zqnorm

%go through sites for calculations
for i = 1:N_Sites
    
    %get indices of values for fit (only bins with meaningful zq values)
    ind_fit = find(~isnan(zq_bin_avg_all{i}));
    
    %get the values to use in fitting
    ust_zqustfit_all{i} = ust_bin_avg_all{i}(ind_fit);
    sigma_ust_zqustfit_all{i} = ust_bin_sigma_all{i}(ind_fit);
    zq_zqustfit_all{i} = zq_bin_avg_all{i}(ind_fit);
    sigma_zq_zqustfit_all{i} = zq_bin_sigma_all{i}(ind_fit);
    zqnorm_zqustfit_all{i} = zqnorm_bin_avg_all{i}(ind_fit);
    sigma_zqnorm_zqustfit_all{i} = zqnorm_bin_sigma_all{i}(ind_fit);
    
    %perform linear fit for each Site
    [intercept, slope, sigma_intercept, sigma_slope, zq_pred, ~] = linearfit(ust_zqustfit_all{i}, zq_zqustfit_all{i}, sigma_zq_zqustfit_all{i});
    intercept_zqustfit_all(i) = intercept; %intercept of fit
    sigma_intercept_zqustfit_all(i) = sigma_intercept; %uncertainty in intercept of fit
    slope_zqustfit_all(i) = slope; %slope of fit
    sigma_slope_zqustfit_all(i) = sigma_slope; %uncertainty in slope of fit
    zq_pred_zqustfit_all{i} = zq_pred; %predicted zq
        
    %compute Chi2
    zq_residuals = zq_pred - zq_zqustfit_all{i};
    Chi2_zqustfit_all(i) = sum((zq_residuals./sigma_zq_zqustfit_all{i}).^2); %compute total Chi2
    df_zqustfit_all(i) = length(zq_residuals)-2;
    
    %get mean, uncertainty of mean, std dev, and std error for zq
    zq_bar_all(i) = mean(zq_zqustfit_all{i});
    [~, sigma_zq_bar] = MeanUncertainty(zq_zqustfit_all{i}, sigma_zq_zqustfit_all{i});
    sigma_zq_bar_all(i) = std(zq_zqustfit_all{i});
    
    %get mean, uncertainty of mean, std dev, and std error for zqnorm
    zqnorm_bar_all(i) = mean(zqnorm_zqustfit_all{i});
    [~, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_zqustfit_all{i}, sigma_zqnorm_zqustfit_all{i});
    sigma_zqnorm_bar_all(i) = std(zqnorm_zqustfit_all{i});
end

%fit values for zq versus ust - Literature
intercept_zqustfit_Lit_all = zeros(N_Lit,1); %intercept of fit
sigma_intercept_zqustfit_Lit_all = zeros(N_Lit,1); %uncertainty in intercept of fit
slope_zqustfit_Lit_all = zeros(N_Lit,1); %slope of fit for zq versus ustar
sigma_slope_zqustfit_Lit_all = zeros(N_Lit,1); %uncertainty in slope for zq versus ustar

zq_bar_Lit_all = zeros(N_Lit,1); %mean of zq (unweighted)
sigma_zq_bar_Lit_all = zeros(N_Lit,1); %std dev of zq

zqnorm_bar_Lit_all = zeros(N_Lit,1); %mean of zqnorm (unweighted)
sigma_zqnorm_bar_Lit_all = zeros(N_Lit,1); %std dev of zqnorm

for i = 1:N_Lit
    [intercept, slope, sigma_intercept, sigma_slope, zq_pred, ~] = linearfit(ust_Lit{i},zq_Lit{i},sigma_zq_Lit{i});
    intercept_zqustfit_Lit_all(i) = intercept; %intercept of fit
    sigma_intercept_zqustfit_Lit_all(i) = sigma_intercept; %uncertainty in intercept of fit
    slope_zqustfit_Lit_all(i) = slope; %slope of fit
    sigma_slope_zqustfit_Lit_all(i) = sigma_slope; %uncertainty in slope of fit
   
    %get mean zq and zqnorm
    zq_bar_Lit_all(i) = mean(zq_Lit{i});
    [~, sigma_zq_bar] = MeanUncertainty(zq_Lit{i}, sigma_zq_Lit{i});
    sigma_zq_bar_Lit_all(i) = std(zq_Lit{i});
    
    %get mean zqnorm
    zqnorm_bar_Lit_all(i) = mean(zqnorm_Lit{i});
    [~, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_Lit{i}, sigma_zqnorm_Lit{i});
    sigma_zqnorm_bar_Lit_all(i) = std(zqnorm_Lit{i});
end

%% FITTING Q VS TAU
%values used in Q-tau fit
ust_fit_all = cell(N_Sites,1);
sigma_ust_fit_all = cell(N_Sites,1);
tau_fit_all = cell(N_Sites,1);
sigma_tau_fit_all = cell(N_Sites,1);
fQ_fit_all = cell(N_Sites,1);
sigma_fQ_fit_all = cell(N_Sites,1);
Q_fit_all = cell(N_Sites,1);
sigma_Q_fit_all = cell(N_Sites,1);

%values resulting from Q-tau linear fit
C_linearfit_all = zeros(N_Sites,1);
sigma_C_linearfit_all = zeros(N_Sites,1);
tauit_linearfit_all = zeros(N_Sites,1);
sigma_tauit_linearfit_all = zeros(N_Sites,1);
ustit_linearfit_all = zeros(N_Sites,1);
sigma_ustit_linearfit_all = zeros(N_Sites,1);
Q_pred_linearfit_all = cell(N_Sites,1); %predicted Q for linear fit
Q_sigmatau_linearfit_all = cell(N_Sites,1); %additional tau contribution to Q uncertainty for linear fit
Q_sigmatotal_linearfit_all = cell(N_Sites,1); %total uncertainty in Q from Q and tau contributions for linear fit
Q_residuals_linearfit_all = cell(N_Sites,1); %get residuals in Q for linear fit
Chi2_linearfit_all = zeros(N_Sites,1);
Chi2_contributions_linearfit_all = cell(N_Sites,1);
df_linearfit_all = zeros(N_Sites,1);

%values resulting from Q-tau three halves fit
tauit_threehalvesfit_all = zeros(N_Sites,1);
sigma_tauit_threehalvesfit_all = zeros(N_Sites,1);
C_threehalvesfit_all = zeros(N_Sites,1);
sigma_C_threehalvesfit_all = zeros(N_Sites,1);
Q_pred_threehalvesfit_all = cell(N_Sites,1);
Q_sigmatau_threehalvesfit_all = cell(N_Sites,1); %additional tau contribution to Q uncertainty for threehalves fit
Q_sigmatotal_threehalvesfit_all = cell(N_Sites,1);  %total uncertainty in Q from Q and tau contributions for threehalves fit
Q_residuals_threehalvesfit_all = cell(N_Sites,1); %get residuals in Q for threehalves fit
Chi2_threehalvesfit_all = zeros(N_Sites,1);
Chi2_contributions_threehalvesfit_all = cell(N_Sites,1);
df_threehalvesfit_all = zeros(N_Sites,1);

for i = 1:N_Sites
    %get indices of values for fit
    ind_fit = find(fQ_bin_avg_all{i}>=fQ_Qtau_fit_min);
    
    %determine which values to use in fitting
    ust_fit_all{i} = ust_bin_avg_all{i}(ind_fit);
    sigma_ust_fit_all{i} = ust_bin_sigma_all{i}(ind_fit);
    tau_fit_all{i} = tau_bin_avg_all{i}(ind_fit);
    sigma_tau_fit_all{i} = tau_bin_sigma_all{i}(ind_fit);
    fQ_fit_all{i} = fQ_bin_avg_all{i}(ind_fit);
    sigma_fQ_fit_all{i} = fQ_bin_sigma_all{i}(ind_fit);
    Q_fit_all{i} = Q_bin_avg_all{i}(ind_fit);
    sigma_Q_fit_all{i} = Q_bin_sigma_all{i}(ind_fit);

    %perform linear fit for each Site
    [intercept, slope, sigma_intercept, sigma_slope, Q_pred, ~] = ...
        linearfit(tau_fit_all{i}, Q_fit_all{i}, sigma_Q_fit_all{i});
    C_linearfit_all(i) = slope;
    sigma_C_linearfit_all(i) = sigma_slope;
    tauit_linearfit_all(i) = -intercept/slope;
    sigma_tauit_linearfit_all(i) = sigma_intercept/slope;
    ustit_linearfit_all(i) = sqrt(tauit_linearfit_all(i)/rho_a(i));
    sigma_ustit_linearfit_all(i) = sigma_tauit_linearfit_all(i)*(1/(2*rho_a(i)))*(1/ustit_linearfit_all(i));
    Q_pred_linearfit_all{i} = Q_pred;
    Q_sigmatau_linearfit_all{i} = slope*sigma_tau_fit_all{i}; %uncertainty in Q due to tau
    Q_sigmatotal_linearfit_all{i} = sqrt(sigma_Q_fit_all{i}.^2+Q_sigmatau_linearfit_all{i}.^2); %total uncertainty in Q for linear fit
    Q_residuals_linearfit_all{i} = Q_fit_all{i}-Q_pred; %get residuals in Q for linear fit
    Chi2_contributions_linearfit_all{i} = (Q_residuals_linearfit_all{i}./Q_sigmatotal_linearfit_all{i}).^2; %compute individual contributions to Chi2 (Bevington and Robinson, Eq. 8.4)
    Chi2_linearfit_all(i) = sum(Chi2_contributions_linearfit_all{i}); %compute total Chi2
    df_linearfit_all(i) = length(ind_fit)-2;

    %perform three-halves fit for each Site
    [tauit, tauit_sigma, C, C_sigma, Chi2, Chi2_contributions, Q_pred] = ...
        ThreeHalvesFluxFit(ust_fit_all{i}, sigma_ust_fit_all{i},...
        tau_fit_all{i}, sigma_tau_fit_all{i},...
        Q_fit_all{i}, sigma_Q_fit_all{i});
    tauit_threehalvesfit_all(i) = tauit;
    sigma_tauit_threehalvesfit_all(i) = tauit_sigma;
    C_threehalvesfit_all(i) = C;
    sigma_C_threehalvesfit_all(i) = C_sigma;
    Q_pred_threehalvesfit_all{i} = Q_pred;
    Q_sigmatau_threehalvesfit_all{i} = sigma_ust_fit_all{i}.*...
        C.*abs(3*tau_fit_all{i}-tauit); %uncertainty in Q due to tau
    Q_sigmatotal_threehalvesfit_all{i} = sqrt(sigma_Q_fit_all{i}.^2+Q_sigmatau_threehalvesfit_all{i}.^2); %total uncertainty in Q for linear fit
    Q_residuals_threehalvesfit_all{i} = Q_fit_all{i}-Q_pred; %get residuals in Q for threehalves fit
    Chi2_threehalvesfit_all(i) = Chi2;
    Chi2_contributions_threehalvesfit_all{i} = Chi2_contributions;
    df_threehalvesfit_all(i) = length(ind_fit)-2;

    %print out fit values
    linear_output = ['linear, \chi^2_{\nu} = ',...
        num2str(Chi2_linearfit_all(i)/df_linearfit_all(i),'%.2f')]
    threehalves_output = ['3/2, \chi^2_{\nu} = ',...
        num2str(Chi2_threehalvesfit_all(i)/df_threehalvesfit_all(i),'%.2f')]
end

%%
%%%%%%%%%%%%%%%%%%
% DERIVED VALUES %
%%%%%%%%%%%%%%%%%%

%% COMPUTE DERIVED VARIABLES
%initialize lists of tauex values for tau bins
tauex_bin_avg_all = cell(N_Sites,1); %avg tauex in bin
tauex_bin_sigma_all = cell(N_Sites,1); %total estimated tauex uncertainty for bin (combination of tau and tauit uncertainty)
tauex_bin_min_all = cell(N_Sites,1); %minimum tauex in bin

%ititialize lists of tau ratio values for tau bins
tauratio_bin_avg_all = cell(N_Sites,1); %avg tau/tauit in bin
tauratio_bin_sigma_all = cell(N_Sites,1); %uncertainty in tau/tauit for bin
tauratio_bin_min_all = cell(N_Sites,1); %determine minimum tau/tauit in bin based on threshold by Site

%initialize CQs 
CQ_bin_avg_all = cell(N_Sites,1); %get average CQ for tau bins
CQ_bin_sigma_all = cell(N_Sites,1); %total estimated CQ uncertainty for bin

%initialize Qnorm for tau bins
Qnorm_bin_avg_all = cell(N_Sites,1); %avg Qnorm in bin
Qnorm_bin_sigma_all = cell(N_Sites,1); %total estimated Qnorm uncertainty for bin (combination of Q and ustit uncertainty)

for i=1:N_Sites
    tauex_bin_avg_all{i} = tau_bin_avg_all{i} - tauit_linearfit_all(i); %determine tauex based on threshold by Site
    tauex_bin_sigma_all{i} = sqrt(tau_bin_sigma_all{i}.^2 + sigma_tauit_linearfit_all(i).^2); %uncertainty in tauex combines tau uncertainty and threshold uncertainty
    tauex_bin_min_all{i} = tau_bin_min_all{i} - tauit_linearfit_all(i); %determine minimum tauex in bin based on threshold by Site
    
    tauratio_bin_avg_all{i} = tau_bin_avg_all{i}/tauit_linearfit_all(i); %determine tauratio based on threshold by Site
    tauratio_bin_sigma_all{i} = sqrt(tau_bin_sigma_all{i}.^2 + (sigma_tauit_linearfit_all(i).*tauratio_bin_avg_all{i}).^2)/tauit_linearfit_all(i); %uncertainty in tauratio combines tau uncertainty and threshold uncertainty
    tauratio_bin_min_all{i} = tau_bin_min_all{i}/tauit_linearfit_all(i); %determine minimum tauratio in bin based on threshold by Site
    
    CQ_bin_avg_all{i} = (Q_bin_avg_all{i}*1e-3)./((1/g)*ustit_linearfit_all(i).*tauex_bin_avg_all{i}); %compute CQ
    CQ_Q_sigma = Q_bin_sigma_all{i}.*(CQ_bin_avg_all{i}./Q_bin_avg_all{i}); %contribution of Q to CQ uncertainty
    CQ_ustit_sigma = sigma_ustit_linearfit_all(i).*(CQ_bin_avg_all{i}./ustit_linearfit_all(i)); %contribution of ustit to CQ uncertainty
    CQ_tauex_sigma = tauex_bin_sigma_all{i}.*(CQ_bin_avg_all{i}./tauex_bin_avg_all{i}); %contribution of tauex to CQ uncertainty
    CQ_bin_sigma_all{i} = sqrt(CQ_Q_sigma.^2+CQ_ustit_sigma.^2+CQ_tauex_sigma.^2); %total CQ uncertainty
    
    Qnorm_bin_avg_all{i} = (Q_bin_avg_all{i}*1e-3).*g./ustit_linearfit_all(i); %avg Qnorm in bin
    Qnorm_bin_sigma_all{i} = Qnorm_bin_avg_all{i}.*sqrt((Q_bin_sigma_all{i}./Q_bin_avg_all{i}).^2+(sigma_ustit_linearfit_all(i)./ustit_linearfit_all(i)).^2); %total estimated Qnorm uncertainty for bin (combination of Q and ustit uncertainty)
    
    %set uncertainties to zero if associated Qnorm=0
    Qnorm_bin_sigma_all{i}(Qnorm_bin_avg_all{i}==0)=0;
end


%% COMPUTE PARAMETERS ASSOCIATED WITH DERIVED VALUES

CQ_all = zeros(N_Sites,1); %CQ for site
sigma_CQ_all = zeros(N_Sites,1); %uncertainty in CQ for site
Ct_all = zeros(N_Sites,1); %Ct/(1-e) for site
sigma_Ct_all = zeros(N_Sites,1); %uncertainty in Ct/(1-e) for site
for i = 1:N_Sites
    ind_CQ = find(tauex_bin_min_all{i}>=N_sigma_tauit*sigma_tauit_linearfit_all(i)); %min tau must exceed threshold by at least 2 sigma
    Qnorm_over_tauex = Qnorm_bin_avg_all{i}(ind_CQ)./tauex_bin_avg_all{i}(ind_CQ); %get ratios of Qnorm / tauex
    CQ_all(i) = mean(Qnorm_over_tauex); %get CQ for site
    sigma_CQ_all(i) = sqrt((1/length(ind_CQ))*sum((Qnorm_over_tauex-mean(Qnorm_over_tauex)).^2)); %get uncertainty in CQ for site
    Q_over_tauex = (Q_bin_avg_all{i}(ind_CQ)*1e-3)./tauex_bin_avg_all{i}(ind_CQ); %get ratios of Q / tauex
    Ct_all(i) = sqrt(g/zq_bar_all(i))*mean(Q_over_tauex); %get Ct/(1-e) for site
    sigma_Ct_all(i) = sqrt((1/length(ind_CQ))*(g/zq_bar_all(i))*sum((Q_over_tauex-mean(Q_over_tauex)).^2)); %get uncertainty in Ct/(1-e) for site
end

%compute fit values for tauex, tauratio, Qnorm, and CQ
tauex_fit_all = cell(N_Sites,1); %tauex's for fit
sigma_tauex_fit_all = cell(N_Sites,1); %sigma tauex's for fit
tauratio_fit_all = cell(N_Sites,1); %tauratio's for fit
sigma_tauratio_fit_all = cell(N_Sites,1); %sigma tauratio's for fit
Qnorm_fit_all = cell(N_Sites,1); %Qnorm's for fit
sigma_Qnorm_fit_all = cell(N_Sites,1); %sigma_Qnorm's for fit
CQ_fit_all = cell(N_Sites,1); %CQ's for fit
sigma_CQ_fit_all = cell(N_Sites,1); %sigma_CQ's for fit

for i = 1:N_Sites
    ind_fit = find(fQ_bin_avg_all{i}>=fQ_Qtau_fit_min);
    tauex_fit_all{i} = tauex_bin_avg_all{i}(ind_fit);
    sigma_tauex_fit_all{i} = tauex_bin_sigma_all{i}(ind_fit);
    tauratio_fit_all{i} = tauratio_bin_avg_all{i}(ind_fit);
    sigma_tauratio_fit_all{i} = tauratio_bin_sigma_all{i}(ind_fit);
    Qnorm_fit_all{i} = Qnorm_bin_avg_all{i}(ind_fit); %Qnorm's for fit
    sigma_Qnorm_fit_all{i} = Qnorm_bin_sigma_all{i}(ind_fit); %sigma_Qnorm's for fit
    CQ_fit_all{i} = CQ_bin_avg_all{i}(ind_fit);
    sigma_CQ_fit_all{i} = CQ_bin_sigma_all{i}(ind_fit);
end

%% GET MEAN VALUE OF C_Q FOR ALL SITES
[CQ_bar_all, sigma_CQ_bar_all] = MeanUncertainty(CQ_all,sigma_CQ_all);

%% save data
save(SaveData_Path,'Site*','N_Sites','*Site','LitNames','N_Lit','*Lit','WindTunnel*','*WindTunnel','*all');