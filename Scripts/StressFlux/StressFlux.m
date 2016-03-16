%% ANALYZE SALTATION FLUX VERSUS SHEAR STRESS RELATIONSHIP

%%
%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%

%initialize
clearvars;
close all;

%parameters
rho_a = 1.18; %air density (kg/m^3)
g = 9.8; %gravity (m/s^2)
z0 = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m) at threshold
zW_limit = 4; %limit on number of unique Wenglor heights in profile

%fitting filters
fQ_Qtau_fit_min = 0.1; %minimum mean fQ for Q versus tau comparison (and determination of threshold)
N_sigma_tauit = 2; %std deviations from impact threshold for minimum tauex
fQ_TFEM_fit_min = [0.1,0.1,0.1]; %mininum for fitting to obtain impact/fluid thresholds
fQ_TFEM_fit_max = [0.95,0.95,0.95]; %maximum for fitting to obtain impact/fluid thresholds

%binning information
bin_N_min = 3; %minimum number of entries for bin
tau_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
fQ_bin_minrange = 0.1; %mininum range of fQ for binning
fQ_bin_maxrange = 0.2; %maximum range of fQ for binning

%information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for storing data output
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
folder_Plots = '../../PlotOutput/StressFlux/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%load grain size data
load(strcat(folder_GrainSizeData,'MeanGrainSize'));
d50_site = d50_surface_site;

%load flux stress window data
load(strcat(folder_AnalysisData,'StressFluxWindows'));
N_Sites = length(Sites);

%load nonlinear flux fit data
load(strcat(folder_AnalysisData,'NonlinearFluxFit'));

%load external data
load(strcat(folder_AnalysisData,'LitData')); %literature data
LitNames = {'Greeley et al. (1996)'; 'Namikas (2003)';'Farrell et al. (2012)'};
ust_lit = {ust_Greeley96, ust_Namikas03, ust_Farrell12};
sigma_ust_lit = {sigma_ust_Greeley96, sigma_ust_Namikas03, sigma_ust_Farrell12};
zq_lit = {zq_Greeley96, zq_Namikas03, zq_Farrell12};
sigma_zq_lit = {sigma_zq_Greeley96, sigma_zq_Namikas03, sigma_zq_Farrell12};
zqnorm_lit = {1e3*zq_Greeley96/d50_Greeley96, 1e3*zq_Namikas03/d50_Namikas03, 1e3*zq_Farrell12/d50_Farrell12};
sigma_zqnorm_lit = {1e3*sigma_zq_Greeley96/d50_Greeley96, 1e3*sigma_zq_Namikas03/d50_Namikas03, 1e3*sigma_zq_Farrell12/d50_Farrell12};
N_lit = 3;

%change zqnorm_lit so that Farrell doesn't show up 
zqnorm_lit{3} = zqnorm_lit{3}*NaN;

%set info for plotting
Markers_field = {'s','d','o','<','>'};
Colors_field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};
MarkerSize_field = 8;
LineWidth_field = 1;
Markers_lit = {'x','+','*'};
MarkerSize_lit = 8;
Colors_lit = {[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
LineWidth_lit = 1;
PlotFont = 14;

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

% remove points with too few Wenglor heights
for i = 1:N_Sites
    ind_zW = find(N_zW_unique_all{i}>=zW_limit);
    for j = 1:N_variables
        eval([variable_list{j},'{i}=',variable_list{j},'{i}(ind_zW);']);
    end
end

%%
%%%%%%%%%%%%%%%%%%%
% PRIMARY BINNING %
%%%%%%%%%%%%%%%%%%%

%% FOR EACH SITE, PERFORM BINNING BY TAU
%initialize lists of tau bins by site
tau_bin_values_all = cell(N_Sites,1); %tau's in each bin
sigma_tau_bin_values_all = cell(N_Sites,1); %sigma_tau's in each bin
tau_bin_avg_all = cell(N_Sites,1); %avg tau in bin
tau_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average tau for bin
tau_bin_std_all = cell(N_Sites,1); %std dev tau in bin
tau_bin_stdalt_all = cell(N_Sites,1); %compute std dev expected based on sigmaavg and N
tau_bin_SE_all = cell(N_Sites,1); %SE tau in bin
tau_bin_sigma_all = cell(N_Sites,1); %total estimated tau uncertainty for bin
tau_bin_min_all = cell(N_Sites,1); %get minimum tau for bin

%number of elements in bin
bin_N_all = cell(N_Sites,1);

%initialize lists of ust values for ust bins
ust_bin_values_all = cell(N_Sites,1); %ust's in each bin
sigma_ust_bin_values_all = cell(N_Sites,1); %sigma_ust's in each bin
ust_bin_avg_all = cell(N_Sites,1); %avg ust in bin
ust_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average ust for bin
ust_bin_std_all = cell(N_Sites,1); %std dev ust in bin
ust_bin_stdalt_all = cell(N_Sites,1); %compute std dev expected based on sigmaavg and N
ust_bin_SE_all = cell(N_Sites,1); %SE ust in bin
ust_bin_sigma_all = cell(N_Sites,1); %total estimated ust uncertainty for bin

%initialize lists of tauex values for tau bins
tauex_bin_avg_all = cell(N_Sites,1); %avg tauex in bin
tauex_bin_sigma_all = cell(N_Sites,1); %total estimated tauex uncertainty for bin (combination of tau and tauit uncertainty)
tauex_bin_min_all = cell(N_Sites,1); %minimum tauex in bin

%separate fQs into tau bins
fQ_bin_values_all = cell(N_Sites,1); %fQs into tau bins
fQ_bin_avg_all = cell(N_Sites,1); %get average fQ for tau bins
fQ_bin_std_all = cell(N_Sites,1); %get std fQ for tau bins
fQ_bin_SE_all = cell(N_Sites,1); %get SE fQ for tau bins
fQ_bin_sigma_all = cell(N_Sites,1); %get total uncertainty for fQ in tau bins

%separate Q's into tau bins
Q_bin_values_all = cell(N_Sites,1); %Q into tau bins
sigma_Q_bin_values_all = cell(N_Sites,1); %sigma_Q into tau bins
Q_bin_avg_all = cell(N_Sites,1); %get average Q for tau bins
Q_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average Q for tau bins
Q_bin_std_all = cell(N_Sites,1); %get std Q for tau bins
Q_bin_stdalt_all = cell(N_Sites,1); %compute std dev expected based on sigmaavg and N
Q_bin_SE_all = cell(N_Sites,1); %get SE Q for tau bins
Q_bin_sigma_all = cell(N_Sites,1); %total estimated Q uncertainty for bin
Q_bin_min_all = cell(N_Sites,1); %get minimum Q in bin

%separate zq's into tau bins
zq_bin_values_all = cell(N_Sites,1); %zq into tau bins
sigma_zq_bin_values_all = cell(N_Sites,1); %sigma_zq into tau bins
zq_bin_avg_all = cell(N_Sites,1); %get average zq for tau bins
zq_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average zq for tau bins
zq_bin_std_all = cell(N_Sites,1); %get std zq for tau bins
zq_bin_stdalt_all = cell(N_Sites,1); %compute std dev expected based on sigmaavg and N
zq_bin_SE_all = cell(N_Sites,1); %get SE zq for tau bins
zq_bin_sigma_all = cell(N_Sites,1); %total estimated zq uncertainty for bin
zq_bin_sigmaust_all = cell(N_Sites,1); %additional ust contribution to zq uncertainty
zq_bin_sigmatotal_all = cell(N_Sites,1);  %total uncertainty in zq from zq and ust contributions

%initialize lists of zqnorms
zqnorm_bin_avg_all = cell(N_Sites,1); %get average zqnorm for tau bins
zqnorm_bin_sigma_all = cell(N_Sites,1); %get uncertainty in zqnorm for tau bins
zqnorm_bin_sigmaust_all = cell(N_Sites,1); %additional ust contribution to zq uncertainty
zqnorm_bin_sigmatotal_all = cell(N_Sites,1); %total uncertainty in zqnorm from zqnorm and ust contributions


%get binned values for all sites
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
    
    %get other bin values
    N_bins = length(tau_bin_values); %number of bins
    sigma_tau_bin_values_all{i} = cell(N_bins,1); %initialize sigma_tau's
    ust_bin_values_all{i} = cell(N_bins,1); %initialize usts
    sigma_ust_bin_values_all{i} = cell(N_bins,1); %initialize sigma_usts
    fQ_bin_values_all{i} = cell(N_bins,1); %initialize fQ's
    Q_bin_values_all{i} = cell(N_bins,1); %initialize Q's
    sigma_Q_bin_values_all{i} = cell(N_bins,1); %initialize sigma_Q's
    Q_bin_min_all{i} = zeros(N_bins,1); %initialize minimum Q's in bins
    zq_bin_values_all{i} = cell(N_bins,1); %initialize zq's
    sigma_zq_bin_values_all{i} = cell(N_bins,1); %initialize sigma_zq's
    for j = 1:N_bins
        bin_ind = find(tau>=tau_bin_min(j)&tau<=tau_bin_max(j)); %get indices of values in bin
        sigma_tau_bin_values_all{i}{j} = sigma_tau(bin_ind); %sigma_tau's
        ust_bin_values_all{i}{j} = ust(bin_ind); %ust's
        sigma_ust_bin_values_all{i}{j} = sigma_ust(bin_ind); %sigma_ust's
        fQ_bin_values_all{i}{j} = fQ(bin_ind); %fQ's
        Q_bin_values_all{i}{j} = Q(bin_ind); %Q's
        sigma_Q_bin_values_all{i}{j} = sigma_Q(bin_ind); %sigma_Q's
        Q_bin_min_all{i}(j) = min(Q(bin_ind)); %minimum Q in bin
        zq_bin_values_all{i}{j} = zq(bin_ind); %zq's
        sigma_zq_bin_values_all{i}{j} = sigma_zq(bin_ind); %sigma_zq's
    end
    
    %get ust avg, std, and SE
    ust_bin_avg_all{i} = cellfun(@mean,ust_bin_values_all{i});
    ust_bin_std_all{i} = cellfun(@std,ust_bin_values_all{i});
    ust_bin_SE_all{i} = ust_bin_std_all{i}./sqrt(bin_N);    
    
    %get fQ avg, std, and SE
    fQ_bin_avg_all{i} = cellfun(@mean,fQ_bin_values_all{i});
    fQ_bin_std_all{i} = cellfun(@std,fQ_bin_values_all{i});
    fQ_bin_SE_all{i} = fQ_bin_std_all{i}./sqrt(bin_N);
    
    %get Q avg, std, and SE
    Q_bin_avg_all{i} = cellfun(@mean,Q_bin_values_all{i});
    Q_bin_std_all{i} = cellfun(@std,Q_bin_values_all{i});
    Q_bin_SE_all{i} = Q_bin_std_all{i}./sqrt(bin_N);
    
    %get zQ avg, std, and SE
    zq_bin_avg_all{i} = cellfun(@mean,zq_bin_values_all{i});
    zq_bin_std_all{i} = cellfun(@std,zq_bin_values_all{i});
    zq_bin_SE_all{i} = zq_bin_std_all{i}./sqrt(bin_N);
    
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
    
    %estimate stdalt for bins (std deviation of bin mean estimated from mean of random errors)
    tau_bin_stdalt_all{i} = tau_bin_sigmaavg_all{i}.*sqrt(bin_N); %alt estimate of tau std dev from sigmaavg and N
    ust_bin_stdalt_all{i} = ust_bin_sigmaavg_all{i}.*sqrt(bin_N); %alt estimate of ust std dev from sigmaavg and N
    Q_bin_stdalt_all{i} = Q_bin_sigmaavg_all{i}.*sqrt(bin_N); %alt estimate of Q std dev from sigmaavg and N
    zq_bin_stdalt_all{i} = zq_bin_sigmaavg_all{i}.*sqrt(bin_N); %alt estimate of zq std dev from sigmaavg and N

    %estimate total sigma for bins
    tau_bin_sigma_all{i} = max([tau_bin_SE_all{i}, tau_bin_sigmaavg_all{i}]')'; %total tau uncertainty (max of SE or uncertainty in mean)
    ust_bin_sigma_all{i} = max([ust_bin_SE_all{i}, ust_bin_sigmaavg_all{i}]')'; %total ust uncertainty (max of SE or uncertainty in mean)
    Q_bin_sigma_all{i} = max([Q_bin_SE_all{i}, Q_bin_sigmaavg_all{i}]')'; %total Q uncertainty (max of SE or uncertainty in mean)
    zq_bin_sigma_all{i} = max([zq_bin_SE_all{i}, zq_bin_sigmaavg_all{i}]')'; %total zq uncertainty (max of SE or uncertainty in mean)
    
    %for fQ, take total uncertainty simply as SE
    fQ_bin_sigma_all{i} = fQ_bin_SE_all{i}; 
end

%% RECOMPUTE SIGMA FOR BINS WITH FEW VALUES

%record original sigma values for comparison to later changed values
tau_bin_sigma_original_all = tau_bin_sigma_all;
ust_bin_sigma_original_all = ust_bin_sigma_all;
fQ_bin_sigma_original_all = fQ_bin_sigma_all;
Q_bin_sigma_original_all = Q_bin_sigma_all;
zq_bin_sigma_original_all = zq_bin_sigma_all;

%recompute sigma for bins with fewer than bin_N_min values based on nearest Oceano tau
i_Oceano = find(strncmp(Sites,'Oceano',6),1); %get index among site list for Oceano
ind_Oceano = find(bin_N_all{i_Oceano}>=bin_N_min); %indices of Oceano bins with bin_N >= bin_N_min
tau_avg_Oceano = tau_bin_avg_all{i_Oceano}(ind_Oceano); %avg tau values for these bins
tau_std_Oceano = max([tau_bin_std_all{i_Oceano}(ind_Oceano)'; tau_bin_stdalt_all{i_Oceano}(ind_Oceano)'])'; %std dev tau values from these bins (from max of std and stdalt)
tau_CV_Oceano = tau_std_Oceano./tau_avg_Oceano; %std/mean of tau values for these bins
ust_avg_Oceano = ust_bin_avg_all{i_Oceano}(ind_Oceano); %avg ust values for these bins
ust_std_Oceano = max([ust_bin_std_all{i_Oceano}(ind_Oceano)'; ust_bin_stdalt_all{i_Oceano}(ind_Oceano)'])'; %std dev ust values from these bins (from max of std and stdalt)
ust_CV_Oceano = ust_std_Oceano./ust_avg_Oceano; %std/mean of ust values for these bins
fQ_avg_Oceano = fQ_bin_avg_all{i_Oceano}(ind_Oceano); %avg fQ values for these bins
fQ_std_Oceano = fQ_bin_std_all{i_Oceano}(ind_Oceano); %std dev fQ values for these bins
fQ_CV_Oceano = fQ_std_Oceano./fQ_avg_Oceano; %std/mean of fQ values for these bins
Q_avg_Oceano = Q_bin_avg_all{i_Oceano}(ind_Oceano); %avg Q values for these bins
Q_std_Oceano = max([Q_bin_std_all{i_Oceano}(ind_Oceano)'; Q_bin_stdalt_all{i_Oceano}(ind_Oceano)'])'; %std dev Q values from these bins (from max of std and stdalt)
Q_CV_Oceano = Q_std_Oceano./Q_avg_Oceano; %std/mean of Q values for these bins
zq_avg_Oceano = zq_bin_avg_all{i_Oceano}(ind_Oceano); %avg zq values for these bins
zq_std_Oceano = max([zq_bin_std_all{i_Oceano}(ind_Oceano)'; zq_bin_stdalt_all{i_Oceano}(ind_Oceano)'])'; %std dev zq values from these bins (from max of std and stdalt)
zq_CV_Oceano = zq_std_Oceano./zq_avg_Oceano; %std/mean of zq values for these bins
for i = 1:N_Sites
    ind_recompute_sigma = find(bin_N_all{i}<bin_N_min);
    N_ind = length(ind_recompute_sigma);
    for j = 1:N_ind
        k = ind_recompute_sigma(j);
        tau_bin = tau_bin_avg_all{i}(k); %get avg tau for bin
        tau_diff_Oceano = abs(tau_bin-tau_avg_Oceano); %get distance of this tau from Oceano taus
        ind_Oceano = find(tau_diff_Oceano==min(tau_diff_Oceano),1); %get closest Oceano tau
        %compute sigma based on this nearest std dev and N values in bin
        tau_bin_sigma_pred = tau_CV_Oceano(ind_Oceano)*tau_bin_avg_all{i}(k)/sqrt(bin_N_all{i}(k));
        ust_bin_sigma_pred = ust_CV_Oceano(ind_Oceano)*ust_bin_avg_all{i}(k)/sqrt(bin_N_all{i}(k));
        fQ_bin_sigma_pred = fQ_CV_Oceano(ind_Oceano)*fQ_bin_avg_all{i}(k)/sqrt(bin_N_all{i}(k));
        Q_bin_sigma_pred = Q_CV_Oceano(ind_Oceano)*Q_bin_avg_all{i}(k)/sqrt(bin_N_all{i}(k));
        zq_bin_sigma_pred = zq_CV_Oceano(ind_Oceano)*zq_bin_avg_all{i}(k)/sqrt(bin_N_all{i}(k));
        %final sigma is max of current value in bin and this new predicted value
        tau_bin_sigma_all{i}(k) = max([tau_bin_sigma_all{i}(k) tau_bin_sigma_pred]);
        ust_bin_sigma_all{i}(k) = max([ust_bin_sigma_all{i}(k) ust_bin_sigma_pred]);
        fQ_bin_sigma_all{i}(k) = max([fQ_bin_sigma_all{i}(k) fQ_bin_sigma_pred]);
        Q_bin_sigma_all{i}(k) = max([Q_bin_sigma_all{i}(k) Q_bin_sigma_pred]);
        zq_bin_sigma_all{i}(k) = max([zq_bin_sigma_all{i}(k) zq_bin_sigma_pred]);
    end
end

%% COMPUTE ZQNORM
for i=1:N_Sites
    zqnorm_bin_avg_all{i} = 1000*zq_bin_avg_all{i}./d50_site(i); % determine zqnorm based on grain size by site
    zqnorm_bin_sigma_all{i} = 1000*zq_bin_sigma_all{i}./d50_site(i); %determine zqnorm based on grain size by site
end

%%
%%%%%%%%%%%
% FITTING %
%%%%%%%%%%%

%% FITTING ZQ VS U*
%values used in zq-ust linear fit
ust_zqustfit_all = cell(N_Sites,1);
ust_sigma_zqustfit_all = cell(N_Sites,1);
zq_zqustfit_all = cell(N_Sites,1);
zq_sigma_zqustfit_all = cell(N_Sites,1);
zqnorm_zqustfit_all = cell(N_Sites,1);
zqnorm_sigma_zqustfit_all = cell(N_Sites,1);
Chi2_zqustfit_all = zeros(N_Sites,1);
df_zqustfit_all = zeros(N_Sites,1);

%values resulting from zq-ust linear fit
intercept_zqustfit_all = zeros(N_Sites,1); %intercept of fit
sigma_intercept_zqustfit_all = zeros(N_Sites,1); %uncertainty in intercept of fit
slope_zqustfit_all = zeros(N_Sites,1); %slope of fit for zq versus ustar
sigma_slope_zqustfit_all = zeros(N_Sites,1); %uncertainty in slope for zq versus ustar
zq_pred_zqustfit_all = cell(N_Sites,1); %predicted zq
zq_pred_relativerange_all = zeros(N_Sites,1); %relative range of zq prediction
sigma_zq_pred_relativerange_all = zeros(N_Sites,1); %uncertainty in relative range of zq prediction
zqnorm_pred_zqustfit_all = cell(N_Sites,1); %predicted zqnorm

%mean values
zq_bar_all = zeros(N_Sites,1); %mean of zq
sigma_zq_bar_all = zeros(N_Sites,1); %uncertainty of zq
zqnorm_bar_all = zeros(N_Sites,1); %mean of zqnorm
sigma_zqnorm_bar_all = zeros(N_Sites,1); %uncertainty of zqnorm bar

for i = 1:N_Sites
    %get indices of values for fit (all fluxes in bin greater than 0)
    ind_fit = find(Q_bin_min_all{i}>0);

    %determine which values to use in fitting
    ust_zqustfit_all{i} = ust_bin_avg_all{i}(ind_fit);
    ust_sigma_zqustfit_all{i} = ust_bin_sigma_all{i}(ind_fit);
    zq_zqustfit_all{i} = zq_bin_avg_all{i}(ind_fit);
    zq_sigma_zqustfit_all{i} = zq_bin_sigma_all{i}(ind_fit);
    zqnorm_zqustfit_all{i} = zqnorm_bin_avg_all{i}(ind_fit);
    zqnorm_sigma_zqustfit_all{i} = zqnorm_bin_sigma_all{i}(ind_fit);

    %perform linear fit for each site
    [intercept, slope, sigma_intercept, sigma_slope, zq_pred, ~] = linearfit(ust_zqustfit_all{i}, zq_zqustfit_all{i}, zq_sigma_zqustfit_all{i});
    intercept_zqustfit_all(i) = intercept; %intercept of fit
    sigma_intercept_zqustfit_all(i) = sigma_intercept; %uncertainty in intercept of fit
    slope_zqustfit_all(i) = slope; %slope of fit
    sigma_slope_zqustfit_all(i) = sigma_slope; %uncertainty in slope of fit
    zq_pred_zqustfit_all{i} = zq_pred; %predicted zq
    zqnorm_pred_zqustfit_all{i} = 1000*zq_pred./d50_site(i); %predicted zq norm
    zq_pred_relativerange_all(i) = sign(slope)*range(zq_pred)./mean(zq_pred); %relative range of zq predictions
    sigma_zq_pred_relativerange_all(i) = abs(sigma_slope/slope)*range(zq_pred)/mean(zq_pred); %relative range of zq predictions
        
    %compute Chi2
    zq_residuals = zq_pred - zq_zqustfit_all{i};
    Chi2_zqustfit_all(i) = sum((zq_residuals./zq_sigma_zqustfit_all{i}).^2); %compute total Chi2
    df_zqustfit_all(i) = length(zq_residuals)-2;
    
    %get mean zq and zqnorm
    [zq_bar, sigma_zq_bar] = MeanUncertainty(zq_zqustfit_all{i}, zq_sigma_zqustfit_all{i});
    zq_bar_all(i) = zq_bar; %mean of zq
    sigma_zq_bar_all(i) = sigma_zq_bar; %uncertainty of zq bar
    
    [zqnorm_bar, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_zqustfit_all{i}, zqnorm_sigma_zqustfit_all{i});
    zqnorm_bar_all(i) = zqnorm_bar; %mean of zqnorm
    sigma_zqnorm_bar_all(i) = sigma_zqnorm_bar; %uncertainty of zqnorm bar
end

%fit values for zq versus ust - literature
intercept_zqust_lit_all = zeros(N_lit,1); %intercept of fit
sigma_intercept_zqust_lit_all = zeros(N_lit,1); %uncertainty in intercept of fit
slope_zqust_lit_all = zeros(N_lit,1); %slope of fit for zq versus ustar
sigma_slope_zqust_lit_all = zeros(N_lit,1); %uncertainty in slope for zq versus ustar
zq_pred_lit_all = cell(N_lit,1); %predicted zq for lit data
zq_pred_relativerange_lit_all = zeros(N_lit,1); %relative range of zq prediction
sigma_zq_pred_relativerange_lit_all = zeros(N_lit,1); %uncertainty in relative range of zq prediction

zqnorm_pred_lit_all = cell(N_lit,1); %predicted zqnorm
zqnorm_bar_lit_all = zeros(N_lit,1); %mean of zqnorm
sigma_zqnorm_bar_lit_all = zeros(N_lit,1); %uncertainty of zqnorm bar

for i = 1:N_lit
    [intercept, slope, sigma_intercept, sigma_slope, zq_pred, ~] = linearfit(ust_lit{i},zq_lit{i},sigma_zq_lit{i});
    intercept_zqust_lit_all(i) = intercept; %intercept of fit
    sigma_intercept_zqust_lit_all(i) = sigma_intercept; %uncertainty in intercept of fit
    slope_zqust_lit_all(i) = slope; %slope of fit
    sigma_slope_zqust_lit_all(i) = sigma_slope; %uncertainty in slope of fit
    zq_pred_lit_all{i} = zq_pred; %predicted zq for lit data
    zq_pred_relativerange_lit_all(i) = sign(slope)*range(zq_pred)/mean(zq_pred);
    sigma_zq_pred_relativerange_lit_all(i) = abs(sigma_slope/slope)*range(zq_pred)/mean(zq_pred); %relative range of zq predictions
    
    [intercept, slope, sigma_intercept, sigma_slope, zqnorm_pred, ~] = linearfit(ust_lit{i},zqnorm_lit{i},sigma_zqnorm_lit{i});
    zqnorm_pred_lit_all{i} = zqnorm_pred;
    [zqnorm_bar, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_lit{i}, sigma_zqnorm_lit{i});
    zqnorm_bar_lit_all(i) = zqnorm_bar;
    sigma_zqnorm_bar_lit_all(i) = sigma_zqnorm_bar;
end

%% FITTING Q VS TAU
%values used in Q-tau fit
ust_Qtaufit_all = cell(N_Sites,1);
ust_sigma_Qtaufit_all = cell(N_Sites,1);
tau_Qtaufit_all = cell(N_Sites,1);
tau_sigma_Qtaufit_all = cell(N_Sites,1);
fQ_Qtaufit_all = cell(N_Sites,1);
Q_Qtaufit_all = cell(N_Sites,1);
Q_sigma_Qtaufit_all = cell(N_Sites,1);

%values resulting from Q-tau linear fit
C_linearfit_all = zeros(N_Sites,1);
C_sigma_linearfit_all = zeros(N_Sites,1);
tauit_linearfit_all = zeros(N_Sites,1);
tauit_sigma_linearfit_all = zeros(N_Sites,1);
ustit_linearfit_all = zeros(N_Sites,1);
ustit_sigma_linearfit_all = zeros(N_Sites,1);
Q_pred_linearfit_all = cell(N_Sites,1); %predicted Q for linear fit
Q_sigmatau_linearfit_all = cell(N_Sites,1); %additional tau contribution to Q uncertainty for linear fit
Q_sigmatotal_linearfit_all = cell(N_Sites,1); %total uncertainty in Q from Q and tau contributions for linear fit
Q_residuals_linearfit_all = cell(N_Sites,1); %get residuals in Q for linear fit
Chi2_linearfit_all = zeros(N_Sites,1);
Chi2_contributions_linearfit_all = cell(N_Sites,1);
df_linearfit_all = zeros(N_Sites,1);
p_linearfit_all = zeros(N_Sites,1);

%values resulting from Q-tau nonlinear fit
n_nonlinearfit_all = zeros(N_Sites,1);
n_range_nonlinearfit_all = zeros(N_Sites,2);
tauit_nonlinearfit_all = zeros(N_Sites,1);
tauit_range_nonlinearfit_all = zeros(N_Sites,2);
C_nonlinearfit_all = zeros(N_Sites,1);
C_range_nonlinearfit_all = zeros(N_Sites,2);
Q_pred_nonlinearfit_all = cell(N_Sites,1);
Q_sigmatau_nonlinearfit_all = cell(N_Sites,1); %additional tau contribution to Q uncertainty for nonlinear fit
Q_sigmatotal_nonlinearfit_all = cell(N_Sites,1);  %total uncertainty in Q from Q and tau contributions for nonlinear fit
Q_residuals_nonlinearfit_all = cell(N_Sites,1); %get residuals in Q for nonlinear fit
Chi2_nonlinearfit_all = zeros(N_Sites,1);
Chi2_contributions_nonlinearfit_all = cell(N_Sites,1);
df_nonlinearfit_all = zeros(N_Sites,1);
p_nonlinearfit_all = zeros(N_Sites,1);
Q_sigmataualt_nonlinearfit_all = cell(N_Sites,1); %alternative computation of sigma tau from associated range in Q_fit for nonlinear fit

%values for plotting fits
tau_linearfit_plot = cell(N_Sites,1);
tau_nonlinearfit_plot = cell(N_Sites,1);
ust_nonlinearfit_plot = cell(N_Sites,1);
Q_linearfit_plot = cell(N_Sites,1);
Q_nonlinearfit_plot = cell(N_Sites,1);

for i = 1:N_Sites
    %get indices of values for fit
    ind_fit = find(fQ_bin_avg_all{i}>=fQ_Qtau_fit_min);
    
    %determine which values to use in fitting
    ust_Qtaufit_all{i} = ust_bin_avg_all{i}(ind_fit);
    ust_sigma_Qtaufit_all{i} = ust_bin_sigma_all{i}(ind_fit);
    tau_Qtaufit_all{i} = tau_bin_avg_all{i}(ind_fit);
    tau_sigma_Qtaufit_all{i} = tau_bin_sigma_all{i}(ind_fit);
    fQ_Qtaufit_all{i} = fQ_bin_avg_all{i}(ind_fit);
    Q_Qtaufit_all{i} = Q_bin_avg_all{i}(ind_fit);
    Q_sigma_Qtaufit_all{i} = Q_bin_sigma_all{i}(ind_fit);

    %perform linear fit for each site
    [intercept, slope, sigma_intercept, sigma_slope, Q_pred, ~] = ...
        linearfit(tau_Qtaufit_all{i}, Q_Qtaufit_all{i}, Q_sigma_Qtaufit_all{i});
    C_linearfit_all(i) = slope;
    C_sigma_linearfit_all(i) = sigma_slope;
    tauit_linearfit_all(i) = -intercept/slope;
    tauit_sigma_linearfit_all(i) = sigma_intercept/slope;
    ustit_linearfit_all(i) = sqrt(tauit_linearfit_all(i)/rho_a);
    ustit_sigma_linearfit_all(i) = tauit_sigma_linearfit_all(i)*(1/(2*rho_a))*(1/ustit_linearfit_all(i));
    Q_pred_linearfit_all{i} = Q_pred;
    Q_sigmatau_linearfit_all{i} = slope*tau_sigma_Qtaufit_all{i}; %uncertainty in Q due to tau
    Q_sigmatotal_linearfit_all{i} = sqrt(Q_sigma_Qtaufit_all{i}.^2+Q_sigmatau_linearfit_all{i}.^2); %total uncertainty in Q for linear fit
    Q_residuals_linearfit_all{i} = Q_Qtaufit_all{i}-Q_pred; %get residuals in Q for linear fit
    Chi2_contributions_linearfit_all{i} = (Q_residuals_linearfit_all{i}./Q_sigmatotal_linearfit_all{i}).^2; %compute individual contributions to Chi2 (Bevington and Robinson, Eq. 8.4)
    Chi2_linearfit_all(i) = sum(Chi2_contributions_linearfit_all{i}); %compute total Chi2
    df_linearfit_all(i) = length(ind_fit)-2;
    p_linearfit_all(i) = chi2cdf(Chi2_linearfit_all(i),df_linearfit_all(i));
    
    %perform nonlinear fit for each site
    [n, n_range, tauit, tauit_range, C, C_range, Chi2, Chi2_contributions, Q_pred] = ...
        NonlinearFluxFit(ust_Qtaufit_all{i}, ust_sigma_Qtaufit_all{i},...
        tau_Qtaufit_all{i}, tau_sigma_Qtaufit_all{i},...
        Q_Qtaufit_all{i}, Q_sigma_Qtaufit_all{i});
    n_nonlinearfit_all(i) = n;
    n_range_nonlinearfit_all(i,:) = n_range;
    tauit_nonlinearfit_all(i) = tauit;
    tauit_range_nonlinearfit_all(i,:) = tauit_range;
    C_nonlinearfit_all(i) = C;
    C_range_nonlinearfit_all(i,:) = C_range;
    Q_pred_nonlinearfit_all{i} = Q_pred;
    Q_sigmatau_nonlinearfit_all{i} = ust_sigma_Qtaufit_all{i}.*...
        abs(C*(rho_a*(n+2)*ust_Qtaufit_all{i}.^(n+1)-...
        n*ust_Qtaufit_all{i}.^(n-1)*tauit)); %uncertainty in Q due to tau
    Q_sigmatotal_nonlinearfit_all{i} = sqrt(Q_sigma_Qtaufit_all{i}.^2+Q_sigmatau_nonlinearfit_all{i}.^2); %total uncertainty in Q for linear fit
    Q_residuals_nonlinearfit_all{i} = Q_Qtaufit_all{i}-Q_pred; %get residuals in Q for nonlinear fit
    Chi2_nonlinearfit_all(i) = Chi2;
    Chi2_contributions_nonlinearfit_all{i} = Chi2_contributions;
    df_nonlinearfit_all(i) = length(ind_fit)-3;
    p_nonlinearfit_all(i) = chi2cdf(Chi2_nonlinearfit_all(i),df_nonlinearfit_all(i));
    
    %alternate calculation of tau contribution to Q
    tau_uncertainty_lower = tau_Qtaufit_all{i} - tau_sigma_Qtaufit_all{i}; %values of tau for lower limit of tau uncertainty range
    tau_uncertainty_upper = tau_Qtaufit_all{i} + tau_sigma_Qtaufit_all{i}; %values of tau for upper limit of tau uncertainty range
    ust_uncertainty_lower = ust_Qtaufit_all{i} - ust_sigma_Qtaufit_all{i}; %values of ust for lower limit of tau uncertainty range
    ust_uncertainty_upper = ust_Qtaufit_all{i} + ust_sigma_Qtaufit_all{i}; %values of ust for upper limit of tau uncertainty range
    Q_pred_uncertainty_lower = C*ust_uncertainty_lower.^n.*(tau_uncertainty_lower-tauit); %values of Q for lower limit of tau uncertainty range
    Q_pred_uncertainty_upper = C*ust_uncertainty_upper.^n.*(tau_uncertainty_upper-tauit); %values of Q for upper limit of tau uncertainty range
    Q_sigmataualt_nonlinearfit_all{i} = (Q_pred_uncertainty_upper-Q_pred_uncertainty_lower)/2; %uncertainty in Q due to tau based on this range

    %get values for plotting fits
    tau_linearfit_plot{i} = linspace(tauit_linearfit_all(i),max(tau_Qtaufit_all{i}),50);
    tau_nonlinearfit_plot{i} = linspace(tauit_nonlinearfit_all(i),max(tau_Qtaufit_all{i}),50);
    ust_nonlinearfit_plot{i} = sqrt(tau_nonlinearfit_plot{i}/rho_a);
    Q_linearfit_plot{i} = C_linearfit_all(i)*(tau_linearfit_plot{i}-tauit_linearfit_all(i));
    Q_nonlinearfit_plot{i} = C_nonlinearfit_all(i)*...
        ust_nonlinearfit_plot{i}.^n_nonlinearfit_all(i).*...
        (tau_nonlinearfit_plot{i}-tauit_nonlinearfit_all(i));
end


%% PERFORM BINNING AND FITS FOR FQ TO GET TFEM THRESHOLD

%initialize fQ and tauth bins
fQ_TFEM_avg_all = cell(N_Sites,1);
fQ_TFEM_sigma_all = cell(N_Sites,1);
tauth_TFEM_avg_all = cell(N_Sites,1);
tauth_TFEM_sigma_all = cell(N_Sites,1);

%initialize lists of best values for fitting
fQ_fQth_TFEM_all = cell(N_Sites,1);
sigma_fQ_fQth_TFEM_all = cell(N_Sites,1);
tauth_fQth_TFEM_all = cell(N_Sites,1);
sigma_tauth_fQth_TFEM_all = cell(N_Sites,1);

%initialize values for plotting best fit line
tauft_TFEM_all = zeros(N_Sites,1);
tauit_TFEM_all = zeros(N_Sites,1);
fQ_fit_TFEM_all = cell(N_Sites,1);
tauth_fit_TFEM_all = cell(N_Sites,1);

%go through each site
for i = 1:N_Sites;
    fQ = fQ_all{i};
    tauth = tauth_TFEM_all{i};

    %get fQ bins
    [fQ_bin_values, bin_N, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_sigma] = ...
        Binning(fQ, fQ_bin_minrange, fQ_bin_maxrange, bin_N_min);
    fQ_TFEM_avg_all{i} = fQ_bin_avg;
    fQ_TFEM_sigma_all{i} = fQ_bin_sigma;
    N_bins = length(fQ_bin_values); %number of bins
  
    %get tauth values
    tauth_TFEM_avg_all{i} = zeros(N_bins,1);
    tauth_TFEM_sigma_all{i} = zeros(N_bins,1);
    for j = 1:N_bins
        %ind_bin = intersect(find(fQ>=fQ_bin_min(j)),find(fQ<=fQ_bin_max(j)));
        ind_bin = intersect(find(fQ>=fQ_bin_min(j)),find(fQ<=fQ_bin_max(j)));
        tauth_bin = tauth(ind_bin);
        tauth_TFEM_avg_all{i}(j) = mean(tauth_bin);
        tauth_TFEM_sigma_all{i}(j) = std(tauth_bin)/sqrt(bin_N(j));
    end
    
    %get values for fitting
    fit_ind = intersect(find(fQ_bin_avg>=fQ_TFEM_fit_min(i)),find(fQ_bin_avg<=fQ_TFEM_fit_max(i)));
    fQ_fQth_TFEM_all{i} = fQ_bin_avg(fit_ind);
    sigma_fQ_fQth_TFEM_all{i} = fQ_bin_sigma(fit_ind);
    tauth_fQth_TFEM_all{i} = tauth_TFEM_avg_all{i}(fit_ind);
    sigma_tauth_fQth_TFEM_all{i} = tauth_TFEM_sigma_all{i}(fit_ind);
    
    %perform fits
    if i==3
        [a, b, ~, ~, ~, ~] = linearfit(fQ_fQth_TFEM_all{i}, tauth_fQth_TFEM_all{i}, sigma_tauth_fQth_TFEM_all{i});
    else %don't use confidence intervals for Jeri and RG, because there are insufficient data for these
        [a, b, ~, ~, ~, ~] = linearfit(fQ_fQth_TFEM_all{i}, tauth_fQth_TFEM_all{i});
    end
    
    %estimate fluid threshold and impact threshold
    tauft_TFEM_all(i) = a;
    tauit_TFEM_all(i) = a+b;
    
    %get fit lines
    fQ_fit_TFEM_all{i} = [0; 1];
    tauth_fit_TFEM_all{i} = [tauft_TFEM_all(i); tauit_TFEM_all(i)];
end


%%
%%%%%%%%%%%%%%%%%%
% DERIVED VALUES %
%%%%%%%%%%%%%%%%%%

%% COMPUTE DERIVED VARIABLES
%initialize Qhats 
Qhat_bin_avg_all = cell(N_Sites,1); %get average Qhat for tau bins
Qhat_bin_sigma_all = cell(N_Sites,1); %total estimated Qhat uncertainty for bin

%initialize alternative Qhats into tau bins
Qhatalt_bin_avg_all = cell(N_Sites,1); %get average alternative Qhat for tau bins
Qhatalt_bin_sigma_all = cell(N_Sites,1); %total estimated alternative Qhat uncertainty for bin

for i=1:N_Sites
    tauex_bin_avg_all{i} = tau_bin_avg_all{i} - tauit_linearfit_all(i); %determine tauex based on threshold by site
    tauex_bin_sigma_all{i} = sqrt(tau_bin_sigma_all{i}.^2 + tauit_sigma_linearfit_all(i).^2); %uncertainty in tauex combines tau uncertainty and threshold uncertainty
    tauex_bin_min_all{i} = tau_bin_min_all{i} - tauit_linearfit_all(i); %determine minimum tauex in bin based on threshold by site
    
    Qhat_bin_avg_all{i} = (Q_bin_avg_all{i}*1e-3)./((1/g)*ustit_linearfit_all(i).*tauex_bin_avg_all{i}); %compute Qhat
    Qhat_Q_sigma = Q_bin_sigma_all{i}.*(Qhat_bin_avg_all{i}./Q_bin_avg_all{i}); %contribution of Q to Qhat uncertainty
    Qhat_ustit_sigma = ustit_sigma_linearfit_all(i).*(Qhat_bin_avg_all{i}./ustit_linearfit_all(i)); %contribution of ustit to Qhat uncertainty
    Qhat_tauex_sigma = tauex_bin_sigma_all{i}.*(Qhat_bin_avg_all{i}./tauex_bin_avg_all{i}); %contribution of tauex to Qhat uncertainty
    Qhat_bin_sigma_all{i} = sqrt(Qhat_Q_sigma.^2+Qhat_ustit_sigma.^2+Qhat_tauex_sigma.^2); %total Qhat uncertainty
    
    Qhatalt_bin_avg_all{i} = (Q_bin_avg_all{i}*1e-3)./(sqrt(zq_bar_all(i)/g).*tauex_bin_avg_all{i}); %compute alternative Qhat
    Qhatalt_Q_sigma = Q_bin_sigma_all{i}.*(Qhatalt_bin_avg_all{i}./Q_bin_avg_all{i}); %contribution of Q to alternative Qhat uncertainty
    Qhatalt_zq_sigma = (1/2)*sigma_zq_bar_all(i).*(Qhatalt_bin_avg_all{i}./zq_bar_all(i)); %contribution of zq to alternative Qhat uncertainty
    Qhatalt_tauex_sigma = tauex_bin_sigma_all{i}.*(Qhatalt_bin_avg_all{i}./tauex_bin_avg_all{i}); %contribution of tauex to alternative Qhat uncertainty
    Qhatalt_bin_sigma_all{i} = sqrt(Qhatalt_Q_sigma.^2+Qhatalt_zq_sigma.^2+Qhatalt_tauex_sigma.^2); %total alternative Qhat uncertainty
end

%% COMPUTE PARAMETERS ASSOCIATED WITH DERIVED VALUES
%compute fit values for Q versus tauex
tauex_Qtauexfit_all = cell(N_Sites,1); %tauex's for fit
sigma_tauex_Qtauexfit_all = cell(N_Sites,1); %sigma tauex's for fit
Q_Qtauexfit_all = cell(N_Sites,1); %Q's for fit
sigma_Q_Qtauexfit_all = cell(N_Sites,1); %sigma_Q's for fit
Qtauratio_bar_all = zeros(N_Sites,1); %scaling coefficient for Q/tauex versus tauex 

Qhat_Qtauexfit_all = cell(N_Sites,1); %Qhat's for fit
sigma_Qhat_Qtauexfit_all = cell(N_Sites,1); %sigma_Qhat's for fit
Qhat_bar_all = zeros(N_Sites,1); %scaling coefficient for Qhat versus tauex
sigma_Qhat_bar_all = zeros(N_Sites,1); %scaling coefficient for Qhat versus tauex

Qhatalt_Qtauexfit_all = cell(N_Sites,1); %alternative Qhat's for fit
sigma_Qhatalt_Qtauexfit_all = cell(N_Sites,1); %alternative Qhat's for fit
Qhatalt_bar_all = zeros(N_Sites,1); %scaling coefficient for alternative Qhat versus tauex
sigma_Qhatalt_bar_all = zeros(N_Sites,1); %scaling coefficient for alternative Qhat versus tauex

%go through sites
for i = 1:N_Sites
    %get values for fitting
    ind_fit = find(tauex_bin_min_all{i}>=N_sigma_tauit*tauit_sigma_linearfit_all(i)); %must exceed threshold by at least 2 sigma
    tauex_Qtauexfit_all{i} = tauex_bin_avg_all{i}(ind_fit);
    sigma_tauex_Qtauexfit_all{i} = tauex_bin_sigma_all{i}(ind_fit);
    Q_Qtauexfit_all{i} = Q_bin_avg_all{i}(ind_fit);
    sigma_Q_Qtauexfit_all{i} = Q_bin_sigma_all{i}(ind_fit);
    Qtauratio_fit = Q_Qtauexfit_all{i}./tauex_Qtauexfit_all{i};
    sigma_Qtauratio_fit = sigma_Q_Qtauexfit_all{i}./tauex_Qtauexfit_all{i};
    Qhat_Qtauexfit_all{i} = Qhat_bin_avg_all{i}(ind_fit);
    sigma_Qhat_Qtauexfit_all{i} = Qhat_bin_sigma_all{i}(ind_fit);
    Qhatalt_Qtauexfit_all{i} = Qhatalt_bin_avg_all{i}(ind_fit);
    sigma_Qhatalt_Qtauexfit_all{i} = Qhatalt_bin_sigma_all{i}(ind_fit);
    
    %get mean values and uncertainties
    [Qtauratio_bar, ~] = MeanUncertainty(Qtauratio_fit, sigma_Qtauratio_fit);
    Qtauratio_bar_all(i) = Qtauratio_bar;
    [Qhat_bar, sigma_Qhat_bar] = MeanUncertainty(Qhat_Qtauexfit_all{i}, sigma_Qhat_Qtauexfit_all{i});
    Qhat_bar_all(i) = Qhat_bar;
    sigma_Qhat_bar_all(i) = sigma_Qhat_bar;
    [Qhatalt_bar, sigma_Qhatalt_bar] = MeanUncertainty(Qhatalt_Qtauexfit_all{i}, sigma_Qhatalt_Qtauexfit_all{i});
    Qhatalt_bar_all(i) = Qhatalt_bar;
    sigma_Qhatalt_bar_all(i) = sigma_Qhatalt_bar;
end

%%
%%%%%%%%%%%%%%%%%
% PRIMARY PLOTS %
%%%%%%%%%%%%%%%%%

%% PLOT ZQ VERSUS UST

%PANEL A - dimensional
figure; subplot(2,1,1); hold on;
legend_items = cell(N_Sites,1);
%plot binned field data
for i = 1:N_Sites
    errorbar(ust_zqustfit_all{i},zq_zqustfit_all{i},zq_sigma_zqustfit_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
    %legend_items{i} = ['\chi^2_{\nu} = ', num2str(Chi2_zqustfit_all(i)/df_zqustfit_all(i),'%.2f')];
end
%plot literature data
for i = 1:N_lit
    errorbar(ust_lit{i},zq_lit{i},sigma_zq_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end
%plot fit for field data
for i = 1:N_Sites;
    plot(ust_zqustfit_all{i},zq_pred_zqustfit_all{i},'Color',Colors_field{i},'LineWidth',LineWidth_field);
end
%plot fit for lit data
for i = 1:N_lit
    plot(ust_lit{i},zq_pred_lit_all{i},'Color',Colors_lit{i},'LineWidth',LineWidth_lit);
end
%organize plot
xlim([0.25 0.6]);
ylim([0 0.15]);
%legend(legend_items,'Location','SouthEast');
text(0.255, 0.14,'(a)','FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
ylabel('\textbf{Saltation height, $$z_q$$ (m)}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%PANEL B - dimensionless
subplot(2,1,2); hold on; %alternative subplots if legend is inside
%plot binned field data
for i = 1:N_Sites
    errorbar(ust_zqustfit_all{i},zqnorm_zqustfit_all{i},zqnorm_sigma_zqustfit_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end
%plot literature data
for i = 1:N_lit        
    errorbar(ust_lit{i},zqnorm_lit{i},sigma_zqnorm_lit{i},Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_lit,'LineWidth',LineWidth_lit);
end
%plot fit for field data
for i = 1:N_Sites
    plot(ust_zqustfit_all{i},zqnorm_pred_zqustfit_all{i},'Color',Colors_field{i},'LineWidth',LineWidth_field);
end
for i = 1:2; %only plot first two entries, ignorning Farrell (2012) with no d50
    plot(ust_lit{i},zqnorm_pred_lit_all{i},'Color',Colors_lit{i},'LineWidth',LineWidth_lit);
end
%organize plot
xlim([0.25 0.6]);
ylim([0 250]);
text(0.255, 235,'(b)','FontSize',PlotFont);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
ylabel('\textbf{Dimensionless salt. ht., $$z_q/d_{50}$$}','Interpreter','Latex');
legend_items = [Sites;LitNames];
legend(legend_items,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9]);
print([folder_Plots,'zq_ust.png'],'-dpng');

%% PLOT Q VERSUS TAUEX BINNED DATA
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    errorbar(tauex_Qtauexfit_all{i},Q_Qtauexfit_all{i},sigma_Q_Qtauexfit_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit, Qtauratio_bar_all(i)*tauex_fit,'Color',Colors_field{i});
end

xlim([0 0.3]);
ylim([0 50]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
legend(Sites,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'Q_tauex.png'],'-dpng');


%% PLOT QHAT VERSUS TAUEX
figure; clf; hold on;

%plot binned values
for i = 1:N_Sites
    errorbar(tauex_Qtauexfit_all{i},Qhat_Qtauexfit_all{i},sigma_Qhat_Qtauexfit_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field,'LineWidth',LineWidth_field);
end

%plot fit values
tauex_fit = [0 0.3];
for i = 1:N_Sites
    plot(tauex_fit,Qhat_bar_all(i)*ones(2,1),'Color',Colors_field{i});
end
xlim([0 0.3]);
ylim([0 12]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Dimensionless saltation flux, $$\hat{Q}$$}','Interpreter','Latex');
legend(Sites,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'Qhat_tauex.png'],'-dpng');

%%
%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT Q VERSUS TAU
figure; hold on;
legend_items = cell(N_Sites*3,1);
for i = 1:N_Sites
    errorbar(tau_Qtaufit_all{i},Q_Qtaufit_all{i},Q_sigmatotal_linearfit_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i}); % include uncertainties from tau
    plot(tau_linearfit_plot{i},Q_linearfit_plot{i},'Color',Colors_field{i}); %plot fit for these values
    plot(tau_nonlinearfit_plot{i},Q_nonlinearfit_plot{i},'--','Color',Colors_field{i}); %plot nonlinear fit
    legend_items{3*i-2} = Sites{i}; %add to list of legend items
    legend_items{3*i-1} = ['linear, \chi^2_{\nu} = ',...
        num2str(Chi2_linearfit_all(i)/df_linearfit_all(i),'%.2f')];
    legend_items{3*i} = ['nonlinear, n=[',...
        num2str(n_range_nonlinearfit_all(i,1),'%.2f'),', ',...
        num2str(n_range_nonlinearfit_all(i,2),'%.2f'),'], \chi^2_{\nu} = ',...
        num2str(Chi2_nonlinearfit_all(i)/df_nonlinearfit_all(i),'%.2f')];
end
xlim([0.05 0.35]);
ylim([0 50]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{$$Q$$ (g/m/s)}','Interpreter','Latex');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition',[0 0 8 6]);
print([folder_Plots,'Q_tau.png'],'-dpng');


%% PLOT RAW ZQ / ZQNORM DATA
figure;
for i = 1:N_Sites

    %zq plot
    subplot(1,N_Sites,i);
    errorbar(ustRe_all{i},zq_all{i},sigma_zq_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    xlim([0.25 0.55]);
    ylim([0 0.14]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Shear velocity, $$u_{*}$$ (m/s)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Saltation height, $$z_q$$ (m)}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'zq_ust_raw.png'],'-dpng');

%% PLOT ALL Q VERSUS ALL TAU
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i);
    errorbar(tauRe_all{i},Q_all{i},sigma_Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{i});
    xlim([0 0.35]);
    ylim([0 50]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$Q$$ (g/m/s)}','Interpreter','Latex');
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 11 4]);
print([folder_Plots,'Q_tau_raw.png'],'-dpng');


%% PLOT Q VERSUS TRANSPORT FREQUENCY (RAW)
figure; hold on;
for i = 1:N_Sites
    errorbar(fQ_all{i},Q_all{i},sigma_Q_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{i});

end
set(gca,'yscale','log');
xlim([0 1]);
ylim([1e-2 1e2]);
legend(Sites,'Location','SouthEast');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Transport frequency, $$f_{q}$$}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 5.5 4]);
print([folder_Plots,'Q_fQ_raw.png'],'-dpng');


%% PLOT Q VERSUS TRANSPORT FREQUENCY (BINNED)
figure; hold on;
for i = 1:N_Sites
    errorbar(fQ_bin_avg_all{i},Q_bin_avg_all{i},Q_bin_sigma_all{i},Markers_field{i},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{i});
end
set(gca,'yscale','log');
xlim([0 1]);
ylim([1e-2 1e2]);
legend(Sites,'Location','SouthEast');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Transport frequency, $$f_{q}$$}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 11 4]);
print([folder_Plots,'Q_fQ.png'],'-dpng');


%% PLOT TRANSPORT FREQUENCY VERSUS SHEAR STRESS (RAW)
figure; hold on;
for i = 1:N_Sites
    plot(tauRe_all{i},fQ_all{i},Markers_field{i},'LineWidth',LineWidth_field,'MarkerSize',MarkerSize_field,'Color',Colors_field{i});
end
xlim([0 0.35]);
ylim([0 1]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear stress, $$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Transport frequency, $$f_{q}$$}','Interpreter','Latex');
legend(Sites,'Location','SouthEast');
set(gca,'FontSize',PlotFont);
    
%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'fQ_tau_raw.png'],'-dpng');


%% PLOT TRANSPORT FREQUENCY VERSUS SHEAR STRESS (BINNED)
figure; hold on;
for i = 1:N_Sites
    errorbar(tau_bin_avg_all{i},fQ_bin_avg_all{i},fQ_bin_sigma_all{i},Markers_field{i},'LineWidth',LineWidth_field,'MarkerSize',MarkerSize_field,'Color',Colors_field{i});
end
xlim([0 0.35]);
ylim([0 1]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear stress, $$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Transport frequency, $$f_{q}$$}','Interpreter','Latex');
legend(Sites,'Location','SouthEast');
    set(gca,'FontSize',PlotFont);
    
%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'fQ_tau.png'],'-dpng');


%% PLOT TAUTH VERSUS FQ
figure; clf; hold on; %initialize plot
for i = 1:N_Sites
    if i~=3
        plot(fQ_fQth_TFEM_all{i},tauth_fQth_TFEM_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i});
    else
        errorbar(fQ_fQth_TFEM_all{i},tauth_fQth_TFEM_all{i},sigma_tauth_fQth_TFEM_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i});
        plot(fQ_fit_TFEM_all{i},tauth_fit_TFEM_all{i},'Color',Colors_field{i});
    end
end
legend_items = Sites;
legend_items{4} = 'Fit to Oceano';
xlabel('frequency of transport');
ylabel('TFEM inferred threshold stress, \tau_{th}');
legend(legend_items,'Location','SouthWest');
set(gca,'FontSize',PlotFont);

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6]);
print([folder_Plots,'tauth_fQ_all.png'],'-dpng');


%%%%%%%%%%%%%%%%%%%%
% DIAGNOSTIC PLOTS %
%%%%%%%%%%%%%%%%%%%%

%% PLOT TAU RELATIVE UNCERTAINTIES VERSUS TAU
figure; hold on;
for i = 1:N_Sites
    plot(tauRe_all{i},sigma_tauRe_all{i}./tauRe_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i});
end
xlim([0 0.35]);
ylim([0 0.12]);
legend(Sites,'Location','NorthEast');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{$$\sigma_{\tau}/\tau$$}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 8 6]);
print([folder_Plots,'sigma_tau.png'],'-dpng');

%% PLOT RELATIVE UNCERTAINTY IN Q VERSUS TOTAL UNCERTAINTY
figure; hold on;
for i = 1:N_Sites
    plot(Q_all{i},sigma_Q_all{i}./Q_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field);
end
xlim([0 55]);
ylim([0 0.3]);
xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
ylabel('\textbf{Wenglor relative uncertainty, $$\sigma_Q/Q$$}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'sigma_Q_Wenglor_profilefit.png'],'-dpng');

%% PLOT RELATIVE UNCERTAINTY IN ZQ VERSUS TOTAL UNCERTAINTY
figure; hold on;
for i = 1:N_Sites
    plot(Q_all{i},sigma_zq_all{i}./zq_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field);
end
xlim([0 55]);
ylim([0 0.15]);
xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
ylabel('\textbf{Wenglor relative uncertainty, $$\sigma_{zq}/{z_q}$$}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'sigma_zq_Wenglor_profilefit.png'],'-dpng');


%% PLOT VARIANTS ON TAU BIN UNCERTAINTY
figure;
for i = 1:N_Sites
    %full plot
    subplot(2,3,i); hold on;
    plot(tau_bin_avg_all{i},tau_bin_sigma_all{i}./tau_bin_avg_all{i},Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(tau_bin_avg_all{i},tau_bin_sigmaavg_all{i}./tau_bin_avg_all{i},Markers_field{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{2});
    plot(tau_bin_avg_all{i},tau_bin_SE_all{i}./tau_bin_avg_all{i},Markers_field{3},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    legend('total uncertainty','random error','standard error','Location','NorthEast');
    title(Sites{i});
    xlim([0 0.4]);
    ylim([0 0.4]);
    plot(xlim,[0.05 0.05],'k-.'); %plot line showing cutoff for lower plot
    xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$\sigma_{\tau}/\tau$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
    
    %replot with zoomed in range
    subplot(2,3,i+3); hold on;
    plot(tau_bin_avg_all{i},tau_bin_sigma_all{i}./tau_bin_avg_all{i},Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(tau_bin_avg_all{i},tau_bin_sigmaavg_all{i}./tau_bin_avg_all{i},Markers_field{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{2});
    plot(tau_bin_avg_all{i},tau_bin_SE_all{i}./tau_bin_avg_all{i},Markers_field{3},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    xlim([0 0.4]);
    ylim([0 0.05]);
    xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$\sigma_{\tau}/\tau$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
print([folder_Plots,'tau_bin_uncertainty.png'],'-dpng');


%% PLOT VARIANTS ON Q BIN UNCERTAINTY
figure;
for i = 1:N_Sites
    %full plot
    subplot(2,3,i); hold on;
    plot(Q_bin_avg_all{i},Q_bin_sigma_all{i}./Q_bin_avg_all{i},Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(Q_bin_avg_all{i},Q_bin_sigmaavg_all{i}./Q_bin_avg_all{i},Markers_field{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{2});
    plot(Q_bin_avg_all{i},Q_bin_SE_all{i}./Q_bin_avg_all{i},Markers_field{3},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    legend('total uncertainty','random error','standard error','Location','NorthEast');
    title(Sites{i});
    xlim([0 50]);
    ylim([0 1]);
    plot(xlim,[0.3 0.3],'k-.'); %plot line showing cutoff for lower plot
    xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
    
    %replot with zoomed in range
    subplot(2,3,i+3); hold on;
    plot(Q_bin_avg_all{i},Q_bin_sigma_all{i}./Q_bin_avg_all{i},Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(Q_bin_avg_all{i},Q_bin_sigmaavg_all{i}./Q_bin_avg_all{i},Markers_field{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{2});
    plot(Q_bin_avg_all{i},Q_bin_SE_all{i}./Q_bin_avg_all{i},Markers_field{3},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    xlim([0 50]);
    ylim([0 0.3]);
    xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
print([folder_Plots,'Q_bin_uncertainty.png'],'-dpng');


%% PLOT CONTRIBUTIONS TO Q BIN UNCERTAINTY
figure;
for i = 1:N_Sites
    %for linear fit
    subplot(2,3,i); hold on;
    plot(Q_Qtaufit_all{i},Q_sigmatotal_linearfit_all{i}./Q_Qtaufit_all{i},Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(Q_Qtaufit_all{i},Q_sigma_Qtaufit_all{i}./Q_Qtaufit_all{i},Markers_field{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{2});
    plot(Q_Qtaufit_all{i},Q_sigmatau_linearfit_all{i}./Q_Qtaufit_all{i},Markers_lit{1},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    legend('total uncertainty','Q contribution','\tau contribution','Location','NorthEast');
    title([Sites{i}, ' - ','Linear Fit']);
    xlim([0 50]);
    ylim([0 0.8]);
    xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
    
    %for nonlinear fit
    subplot(2,3,i+3); hold on;
    plot(Q_Qtaufit_all{i},Q_sigmatotal_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(Q_Qtaufit_all{i},Q_sigma_Qtaufit_all{i}./Q_Qtaufit_all{i},Markers_field{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{2});
    plot(Q_Qtaufit_all{i},Q_sigmatau_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_lit{1},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    plot(Q_Qtaufit_all{i},Q_sigmataualt_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_lit{2},'MarkerSize',MarkerSize_field/2,'Color',Colors_field{3});
    legend('total uncertainty','Q contribution','\tau contribution','\tau contr. (alt calc)','Location','NorthEast');
    title([Sites{i}, ' - ','Nonlinear Fit']);
    xlim([0 50]);
    ylim([0 0.8]);
    xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
    ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
end
%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
print([folder_Plots,'Q_uncertainty_contributions.png'],'-dpng');


%% PLOT CONTRIBUTIONS TO CHI2 UNCERTAINTY
figure;
for i = 1:N_Sites
    %for linear fit
    subplot(1,3,i); hold on;
    plot(tau_Qtaufit_all{i},Chi2_contributions_linearfit_all{i}/df_linearfit_all(i),Markers_field{1},'MarkerSize',MarkerSize_field,'Color',Colors_field{1});
    plot(tau_Qtaufit_all{i},Chi2_contributions_nonlinearfit_all{i}/df_nonlinearfit_all(i),Markers_field{2},'MarkerSize',MarkerSize_field,'Color',Colors_field{2});
    legend('Linear fit','Nonlinear fit','Location','NorthEast');
    title(Sites{i});
    xlim([0 0.4]);
    ylim([0 0.5]);
    xlabel('\textbf{$$\tau$$}','Interpreter','Latex');
    ylabel('\textbf{$$\chi^2_{\nu}$$}','Interpreter','Latex');
    set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'Chi2_contributions_Qtaufit.png'],'-dpng');


%% PLOT THETA MINUS MEAN THETA VERSUS TAU
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(tauRe_all{i}(Q_all{i}>0),theta_all{i}(Q_all{i}>0)-mean(theta_all{i}(Q_all{i}>0)),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    plot(tauRe_all{i}(Q_all{i}==0),theta_all{i}(Q_all{i}==0)-mean(theta_all{i}(Q_all{i}>0)),Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    legend('transport','no transport');
    xlim([0 0.4]);
    ylim([-90 90]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Shear stress, $$\tau$$ (Pa)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Wind angle, $$\theta (^{\circ})$$}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'theta_tau.png'],'-dpng');


%% PLOT z/L VERSUS TAU
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(tauRe_all{i}(Q_all{i}>0),zL_all{i}(Q_all{i}>0),Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    plot(tauRe_all{i}(Q_all{i}==0),zL_all{i}(Q_all{i}==0),Markers_lit{i},'Color',Colors_lit{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    legend('transport','no transport','Location','SouthEast');
    xlim([0 0.4]);
    ylim([-0.5 0]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Shear stress, $$\tau$$ (Pa)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Stability parameter, $$z/L$$}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'zL_tau.png'],'-dpng');


%% PLOT CHI2 FOR WENGLOR PROFILE VERSUS Q
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(Q_all{i},Chi2_Qfit_all{i}./df_Qfit_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    xlim([0 50]);
    ylim([1e-3 1e2]);
    set(gca,'yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
%    xlabel('\textbf{Saltation flux, $$\tau$$ (Pa)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Quality of fit for Wenglor profile, $$\chi^2_{\nu}$$}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'Chi2_qz_profilefit.png'],'-dpng');


%% PLOT QCAL VERSUS Z
ID_unique = unique(horzcat(W_ID_all{3}{Q_all{3}>0})); %get Wenglor IDs from Oceano
N_ID = length(ID_unique);
qcal_plot_markers = {'r+','go','b*','cx','ms','kd','r^','gv','b>','c>','mp','kh','b+','co'};

figure;
for i = 1:N_Sites
    %identify time intervals with transport
    ind_trans = find(Q_all{i}>0);
    N_trans = length(ind_trans);
    
    %initialize list of calibration coefficients and heights
    qcal_ID = cell(N_ID,1);
    zW_ID = cell(N_ID,1);
    
    %sort calibration values by Wenglor ID
    for j = 1:N_ID;
        qcal_list = []; %initialize list of qcal
        zW_list = []; %initialize list of zW
        for k = 1:N_trans; %go through each interval to get values for list
            zW_interval = zW_all{i}{ind_trans(k)};
            qcal_interval = qcal_all{i}{ind_trans(k)};
            ID_interval = W_ID_all{i}{ind_trans(k)};
            ind_ID = find(strcmp(ID_interval,ID_unique(j)));
            if length(ind_ID)==1
                qcal_list = [qcal_list; qcal_interval(ind_ID)];
                zW_list = [zW_list; zW_interval(ind_ID)];
            end
        end
        qcal_ID{j} = qcal_list;
        zW_ID{j} = zW_list;
    end
    
    %plot
    if i~=3
        subplot(1,4,i); hold on;
    else
        subplot(1,4,3:4); hold on;
    end
    for j = 1:N_ID
        plot(zW_ID{j},qcal_ID{j},qcal_plot_markers{j});
    end
    if i==3
        legend(ID_unique,'Location','EastOutside');
    end
    
    title(Sites{i});
    ylim([1 1000]);
    set(gca,'yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Wenglor ht., $$z$$ (m)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Calibration factor, $$q/n$$ (g m$$^{-1}$$ s$$^{-1}$$ count$$^{-1}$$)}','Interpreter','Latex');
    end
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 5]);
print([folder_Plots,'qcal_z.png'],'-dpng');


%% PLOT NORMALIZED RESIDUALS FOR FIT

figure;
subplot(2,1,1); hold on; %residuals for linear fit
for i = 1:N_Sites
    errorbar(tau_Qtaufit_all{i},Q_residuals_linearfit_all{i}./Q_Qtaufit_all{i},Q_sigmatotal_linearfit_all{i}./Q_Qtaufit_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i}); % include uncertainties from tau
end
xlim([0 0.4]);
ylim([-1 1]);
plot(xlim,[0 0]); %plot zero line
title('Linear fit');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{$$(Q - Q_{pred})/Q$$}','Interpreter','Latex');
legend(Sites,'Location','NorthEast');
set(gca,'FontSize',PlotFont);
    
subplot(2,1,2); hold on; %residuals for nonlinear fit
for i = 1:N_Sites
    errorbar(tau_Qtaufit_all{i},Q_residuals_nonlinearfit_all{i}./Q_Qtaufit_all{i},Q_sigmatotal_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_field{i},'MarkerSize',MarkerSize_field,'Color',Colors_field{i}); % include uncertainties from tau
end
xlim([0 0.4]);
ylim([-1 1]);
plot(xlim,[0 0]); %plot zero line
title('Nonlinear fit');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{$$(Q - Q_{pred})/Q$$}','Interpreter','Latex');
legend(Sites,'Location','NorthEast');
set(gca,'FontSize',PlotFont);

set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 8 10]);
print([folder_Plots,'Qresiduals_tau.png'],'-dpng');

%%
%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%

path_SaveData = strcat(folder_AnalysisData,'StressFluxWindows_Analysis');
save(path_SaveData,'*all');