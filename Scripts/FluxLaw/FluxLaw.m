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

%load data
folder_LoadData = '../../AnalysisData/FluxLaw/'; %folder for outputs of this analysis
LoadData_Path = strcat(folder_LoadData,'FluxLawWindows_30min'); %path for 30 minute data
load(LoadData_Path);
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions
N_Sites = length(Sites);

%binning information
bin_N_min = 3; %minimum number of entries for bin
tau_bin_minrange = 0.01; %minimum difference between upper and lower value in bin
tau_bin_maxrange = 0.025; %maximum difference between upper and lower value in bin
fQ_bin_minrange = 0.1; %mininum range of fQ for binning
fQ_bin_maxrange = 0.2; %maximum range of fQ for binning

%load grain size data
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
load(strcat(folder_GrainSizeData,'MeanGrainSize'));
d50_Site = d50_surface_site;
sigma_d50_Site = sigma_d50_surface_site;

%load external data
folder_LitData = '../../AnalysisData/Literature/'; %folder for loading lit data
load(strcat(folder_LitData,'LitData')); %Literature data
LitNames = {'Greeley et al. (1996)'; 'Namikas (2003)';'Farrell et al. (2012)'};
ust_Lit = {ust_Greeley96, ust_Namikas03, ust_Farrell12};
sigma_ust_Lit = {sigma_ust_Greeley96, sigma_ust_Namikas03, sigma_ust_Farrell12};
zq_Lit = {zq_Greeley96, zq_Namikas03, zq_Farrell12};
sigma_zq_Lit = {sigma_zq_Greeley96, sigma_zq_Namikas03, sigma_zq_Farrell12};
zqnorm_Lit = {1e3*zq_Greeley96/d50_Greeley96, 1e3*zq_Namikas03/d50_Namikas03, 1e3*zq_Farrell12/d50_Farrell12};
sigma_zqnorm_Lit = {1e3*sigma_zq_Greeley96/d50_Greeley96, 1e3*sigma_zq_Namikas03/d50_Namikas03, 1e3*sigma_zq_Farrell12/d50_Farrell12};
d50_Lit = [d50_Greeley96 d50_Namikas03 d50_Farrell12];
N_Lit = 3;
zqnorm_Lit{3} = zqnorm_Lit{3}*NaN; %change zqnorm_Lit so that Farrell doesn't show up

%information about where to save data and plots
folder_SaveData = '../../AnalysisData/FluxLaw/'; %folder for storing data output
path_SaveData = strcat(folder_SaveData,'FluxLawData_Binned');
folder_Plots = '../../PlotOutput/FluxLaw/'; %folder for plots

%set info for plotting
Markers_Field = {'s','d','o','<','>'};
Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};
MarkerSize_Field = 5;
LineWidth_Field = 1;
Markers_Lit = {'<','>','^'};
MarkerSize_Lit = 5;
Colors_Lit = {[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
LineWidth_Lit = 1;
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
%initialize lists of tau bins by Site
tau_bin_values_all = cell(N_Sites,1); %tau's in each bin
sigma_tau_bin_values_all = cell(N_Sites,1); %sigma_tau's in each bin
tau_bin_avg_all = cell(N_Sites,1); %avg tau in bin
tau_bin_sigmaavg_all = cell(N_Sites,1); %get uncertainty in average tau for bin
tau_bin_std_all = cell(N_Sites,1); %std dev tau in bin
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

%separate Q's into tau bins
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
    
    %get other bin values
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
    ind_transport = find(fQ_bin_max_all{i}>=fQ_Qtau_fit_min); %indices for bins with transport
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

%mean values - zq
zq_bar_all = zeros(N_Sites,1); %mean of zq
zq_sigmaavg_all = zeros(N_Sites,1); %uncertainty of zq mean
zq_std_all = zeros(N_Sites,1); %std dev of zq
zq_SE_all = zeros(N_Sites,1); %standard error for zq

%mean values - zqnorm
zqnorm_bar_all = zeros(N_Sites,1); %mean of zqnorm
zqnorm_sigmaavg_all = zeros(N_Sites,1); %uncertainty of zq mean
zqnorm_std_all = zeros(N_Sites,1); %std dev of zqnorm
zqnorm_SE_all = zeros(N_Sites,1); %standard error for zqnorm

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
    
    %perform linear fit for each Site
    [intercept, slope, sigma_intercept, sigma_slope, zq_pred, ~] = linearfit(ust_zqustfit_all{i}, zq_zqustfit_all{i}, zq_sigma_zqustfit_all{i});
    intercept_zqustfit_all(i) = intercept; %intercept of fit
    sigma_intercept_zqustfit_all(i) = sigma_intercept; %uncertainty in intercept of fit
    slope_zqustfit_all(i) = slope; %slope of fit
    sigma_slope_zqustfit_all(i) = sigma_slope; %uncertainty in slope of fit
    zq_pred_zqustfit_all{i} = zq_pred; %predicted zq
        
    %compute Chi2
    zq_residuals = zq_pred - zq_zqustfit_all{i};
    Chi2_zqustfit_all(i) = sum((zq_residuals./zq_sigma_zqustfit_all{i}).^2); %compute total Chi2
    df_zqustfit_all(i) = length(zq_residuals)-2;
    
    %get mean, uncertainty of mean, std dev, and std error for zq
    zq_bar_all(i) = mean(zq_zqustfit_all{i});
    [~, sigma_zq_bar] = MeanUncertainty(zq_zqustfit_all{i}, zq_sigma_zqustfit_all{i});
    zq_sigmaavg_all(i) = sigma_zq_bar; %uncertainty of zq bar
    zq_std_all(i) = std(zq_zqustfit_all{i});
    zq_SE_all(i) = zq_std_all(i)/sqrt(length(zq_zqustfit_all{i}));
    
    %get mean, uncertainty of mean, std dev, and std error for zqnorm
    zqnorm_bar_all(i) = mean(zqnorm_zqustfit_all{i});
    [~, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_zqustfit_all{i}, zqnorm_sigma_zqustfit_all{i});
    zqnorm_sigmaavg_all(i) = sigma_zqnorm_bar; %uncertainty of zqnorm bar
    zqnorm_std_all(i) = std(zqnorm_zqustfit_all{i});
    zqnorm_SE_all(i) = zqnorm_std_all(i)/sqrt(length(zq_zqustfit_all{i}));
end

%fit values for zq versus ust - Literature
intercept_zqust_Lit_all = zeros(N_Lit,1); %intercept of fit
sigma_intercept_zqust_Lit_all = zeros(N_Lit,1); %uncertainty in intercept of fit
slope_zqust_Lit_all = zeros(N_Lit,1); %slope of fit for zq versus ustar
sigma_slope_zqust_Lit_all = zeros(N_Lit,1); %uncertainty in slope for zq versus ustar

zq_bar_Lit_all = zeros(N_Lit,1); %mean of zq (unweighted)
zq_sigmaavg_Lit_all = zeros(N_Lit,1); %uncertainty of mean of zq
zq_SE_Lit_all = zeros(N_Lit,1); %standard error of zq
zq_std_Lit_all = zeros(N_Lit,1); %std dev of zq

zqnorm_bar_Lit_all = zeros(N_Lit,1); %mean of zqnorm (unweighted)
zqnorm_sigmaavg_Lit_all = zeros(N_Lit,1); %uncertainty of zqnorm bar
zqnorm_SE_Lit_all = zeros(N_Lit,1); %standard error of zqnorm
zqnorm_std_Lit_all = zeros(N_Lit,1); %std dev of zqnorm

for i = 1:N_Lit
    [intercept, slope, sigma_intercept, sigma_slope, zq_pred, ~] = linearfit(ust_Lit{i},zq_Lit{i},sigma_zq_Lit{i});
    intercept_zqust_Lit_all(i) = intercept; %intercept of fit
    sigma_intercept_zqust_Lit_all(i) = sigma_intercept; %uncertainty in intercept of fit
    slope_zqust_Lit_all(i) = slope; %slope of fit
    sigma_slope_zqust_Lit_all(i) = sigma_slope; %uncertainty in slope of fit
   
    %get mean zq and zqnorm
    zq_bar_Lit_all(i) = mean(zq_Lit{i});
    [~, sigma_zq_bar] = MeanUncertainty(zq_Lit{i}, sigma_zq_Lit{i});
    zq_sigmaavg_Lit_all(i) = sigma_zq_bar; %uncertainty of zq bar
    zq_std_Lit_all(i) = std(zq_Lit{i});
    zq_SE_Lit_all(i) = zq_std_Lit_all(i)/sqrt(length(zq_Lit{i}));
    
    %get mean zqnorm
    zqnorm_bar_Lit_all(i) = mean(zqnorm_Lit{i});
    [~, sigma_zqnorm_bar] = MeanUncertainty(zqnorm_Lit{i}, sigma_zqnorm_Lit{i});
    zqnorm_sigmaavg_Lit_all(i) = sigma_zqnorm_bar;
    zqnorm_std_Lit_all(i) = std(zqnorm_Lit{i});
    zqnorm_SE_Lit_all(i) = zqnorm_std_Lit_all(i)/sqrt(length(zqnorm_Lit{i}));
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
tau_linearfit_plot = cell(N_Sites,1); %values for plotting fits
Q_linearfit_plot = cell(N_Sites,1); %values for plotting fits

%values resulting from Q-tau three halves fit
tauit_threehalvesfit_all = zeros(N_Sites,1);
tauit_sigma_threehalvesfit_all = zeros(N_Sites,2);
C_threehalvesfit_all = zeros(N_Sites,1);
C_sigma_threehalvesfit_all = zeros(N_Sites,2);
Q_pred_threehalvesfit_all = cell(N_Sites,1);
Q_sigmatau_threehalvesfit_all = cell(N_Sites,1); %additional tau contribution to Q uncertainty for threehalves fit
Q_sigmatotal_threehalvesfit_all = cell(N_Sites,1);  %total uncertainty in Q from Q and tau contributions for threehalves fit
Q_residuals_threehalvesfit_all = cell(N_Sites,1); %get residuals in Q for threehalves fit
Chi2_threehalvesfit_all = zeros(N_Sites,1);
Chi2_contributions_threehalvesfit_all = cell(N_Sites,1);
df_threehalvesfit_all = zeros(N_Sites,1);
tau_threehalvesfit_plot = cell(N_Sites,1); %values for plotting fits
ust_threehalvesfit_plot = cell(N_Sites,1); %values for plotting fits
Q_threehalvesfit_plot = cell(N_Sites,1); %values for plotting fits

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

    %perform linear fit for each Site
    [intercept, slope, sigma_intercept, sigma_slope, Q_pred, ~] = ...
        linearfit(tau_Qtaufit_all{i}, Q_Qtaufit_all{i}, Q_sigma_Qtaufit_all{i});
    C_linearfit_all(i) = slope;
    C_sigma_linearfit_all(i) = sigma_slope;
    tauit_linearfit_all(i) = -intercept/slope;
    tauit_sigma_linearfit_all(i) = sigma_intercept/slope;
    ustit_linearfit_all(i) = sqrt(tauit_linearfit_all(i)/rho_a(i));
    ustit_sigma_linearfit_all(i) = tauit_sigma_linearfit_all(i)*(1/(2*rho_a(i)))*(1/ustit_linearfit_all(i));
    Q_pred_linearfit_all{i} = Q_pred;
    Q_sigmatau_linearfit_all{i} = slope*tau_sigma_Qtaufit_all{i}; %uncertainty in Q due to tau
    Q_sigmatotal_linearfit_all{i} = sqrt(Q_sigma_Qtaufit_all{i}.^2+Q_sigmatau_linearfit_all{i}.^2); %total uncertainty in Q for linear fit
    Q_residuals_linearfit_all{i} = Q_Qtaufit_all{i}-Q_pred; %get residuals in Q for linear fit
    Chi2_contributions_linearfit_all{i} = (Q_residuals_linearfit_all{i}./Q_sigmatotal_linearfit_all{i}).^2; %compute individual contributions to Chi2 (Bevington and Robinson, Eq. 8.4)
    Chi2_linearfit_all(i) = sum(Chi2_contributions_linearfit_all{i}); %compute total Chi2
    df_linearfit_all(i) = length(ind_fit)-2;

    %perform three-halves fit for each Site
    [tauit, tauit_range, C, C_range, Chi2, Chi2_contributions, Q_pred] = ...
        ThreeHalvesFluxFit(ust_Qtaufit_all{i}, ust_sigma_Qtaufit_all{i},...
        tau_Qtaufit_all{i}, tau_sigma_Qtaufit_all{i},...
        Q_Qtaufit_all{i}, Q_sigma_Qtaufit_all{i});
    tauit_threehalvesfit_all(i) = tauit;
    tauit_sigma_threehalvesfit_all(i,:) = tauit_range;
    C_threehalvesfit_all(i) = C;
    C_sigma_threehalvesfit_all(i,:) = C_range;
    Q_pred_threehalvesfit_all{i} = Q_pred;
    Q_sigmatau_threehalvesfit_all{i} = ust_sigma_Qtaufit_all{i}.*...
        C.*abs(3*tau_Qtaufit_all{i}-tauit); %uncertainty in Q due to tau
    Q_sigmatotal_threehalvesfit_all{i} = sqrt(Q_sigma_Qtaufit_all{i}.^2+Q_sigmatau_threehalvesfit_all{i}.^2); %total uncertainty in Q for linear fit
    Q_residuals_threehalvesfit_all{i} = Q_Qtaufit_all{i}-Q_pred; %get residuals in Q for threehalves fit
    Chi2_threehalvesfit_all(i) = Chi2;
    Chi2_contributions_threehalvesfit_all{i} = Chi2_contributions;
    df_threehalvesfit_all(i) = length(ind_fit)-2;
      
    %get values for plotting fits
    tau_linearfit_plot{i} = linspace(tauit_linearfit_all(i),max(tau_Qtaufit_all{i}),50);
    Q_linearfit_plot{i} = C_linearfit_all(i)*(tau_linearfit_plot{i}-tauit_linearfit_all(i));
    tau_threehalvesfit_plot{i} = linspace(tauit_threehalvesfit_all(i),max(tau_Qtaufit_all{i}),50);
    ust_threehalvesfit_plot{i} = sqrt(tau_threehalvesfit_plot{i}/rho_a(i));
    Q_threehalvesfit_plot{i} = C_threehalvesfit_all(i)*...
        ust_threehalvesfit_plot{i}.*...
        (tau_threehalvesfit_plot{i}-tauit_threehalvesfit_all(i));
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

%initialize tauratio
tauratio_bin_avg_all = cell(N_Sites,1); %avg tau/tauit in bin
tauratio_bin_sigma_all = cell(N_Sites,1); %uncertainty in tau/tauit for bin
tauratio_bin_min_all = cell(N_Sites,1); %determine minimum tau/tauit in bin based on threshold by Site

%initialize CQs 
CQ_bin_avg_all = cell(N_Sites,1); %get average CQ for tau bins
CQ_bin_sigma_all = cell(N_Sites,1); %total estimated CQ uncertainty for bin

%initialize alternative CQs into tau bins
CQalt_bin_avg_all = cell(N_Sites,1); %get average alternative CQ for tau bins
CQalt_bin_sigma_all = cell(N_Sites,1); %total estimated alternative CQ uncertainty for bin

for i=1:N_Sites
    tauex_bin_avg_all{i} = tau_bin_avg_all{i} - tauit_linearfit_all(i); %determine tauex based on threshold by Site
    tauex_bin_sigma_all{i} = sqrt(tau_bin_sigma_all{i}.^2 + tauit_sigma_linearfit_all(i).^2); %uncertainty in tauex combines tau uncertainty and threshold uncertainty
    tauex_bin_min_all{i} = tau_bin_min_all{i} - tauit_linearfit_all(i); %determine minimum tauex in bin based on threshold by Site

    tauratio_bin_avg_all{i} = tau_bin_avg_all{i}/tauit_linearfit_all(i); %determine tauratio based on threshold by Site
    tauratio_bin_sigma_all{i} = sqrt(tau_bin_sigma_all{i}.^2 + (tauit_sigma_linearfit_all(i).*tauratio_bin_avg_all{i}).^2)/tauit_linearfit_all(i); %uncertainty in tauratio combines tau uncertainty and threshold uncertainty
    tauratio_bin_min_all{i} = tau_bin_min_all{i}/tauit_linearfit_all(i); %determine minimum tauratio in bin based on threshold by Site
    
    CQ_bin_avg_all{i} = (Q_bin_avg_all{i}*1e-3)./((1/g)*ustit_linearfit_all(i).*tauex_bin_avg_all{i}); %compute CQ
    CQ_Q_sigma = Q_bin_sigma_all{i}.*(CQ_bin_avg_all{i}./Q_bin_avg_all{i}); %contribution of Q to CQ uncertainty
    CQ_ustit_sigma = ustit_sigma_linearfit_all(i).*(CQ_bin_avg_all{i}./ustit_linearfit_all(i)); %contribution of ustit to CQ uncertainty
    CQ_tauex_sigma = tauex_bin_sigma_all{i}.*(CQ_bin_avg_all{i}./tauex_bin_avg_all{i}); %contribution of tauex to CQ uncertainty
    CQ_bin_sigma_all{i} = sqrt(CQ_Q_sigma.^2+CQ_ustit_sigma.^2+CQ_tauex_sigma.^2); %total CQ uncertainty
    
    CQalt_bin_avg_all{i} = (Q_bin_avg_all{i}*1e-3)./(sqrt(zq_bar_all(i)/g).*tauex_bin_avg_all{i}); %compute alternative CQ
    CQalt_Q_sigma = Q_bin_sigma_all{i}.*(CQalt_bin_avg_all{i}./Q_bin_avg_all{i}); %contribution of Q to alternative CQ uncertainty
    CQalt_zq_sigma = (1/2)*zq_sigmaavg_all(i).*(CQalt_bin_avg_all{i}./zq_bar_all(i)); %contribution of zq to alternative CQ uncertainty
    CQalt_tauex_sigma = tauex_bin_sigma_all{i}.*(CQalt_bin_avg_all{i}./tauex_bin_avg_all{i}); %contribution of tauex to alternative CQ uncertainty
    CQalt_bin_sigma_all{i} = sqrt(CQalt_Q_sigma.^2+CQalt_zq_sigma.^2+CQalt_tauex_sigma.^2); %total alternative CQ uncertainty
end

%% COMPUTE PARAMETERS ASSOCIATED WITH DERIVED VALUES

%initialize fit values for tauex and tauratio
tauex_fit_all = cell(N_Sites,1); %tauex's for fit
sigma_tauex_fit_all = cell(N_Sites,1); %sigma tauex's for fit
tauratio_fit_all = cell(N_Sites,1); %tauratio's for fit
sigma_tauratio_fit_all = cell(N_Sites,1); %sigma tauratio's for fit

%initialize values for Q, CQ, and CQalt
Q_fit_all = cell(N_Sites,1); %Q's for fit
sigma_Q_fit_all = cell(N_Sites,1); %sigma_Q's for fit
CQ_fit_all = cell(N_Sites,1); %CQ's for fit
sigma_CQ_fit_all = cell(N_Sites,1); %sigma_CQ's for fit
CQalt_fit_all = cell(N_Sites,1); %alternative CQ's for fit
sigma_CQalt_fit_all = cell(N_Sites,1); %alternative CQ's for fit

%initialize Q/tauex ratio, CQ, and CQalt calculations
Qtauexratio_bar_all = zeros(N_Sites,1); %mean Q/tauex ratio
CQ_bar_all = zeros(N_Sites,1); %mean CQ
CQ_sigmaavg_all = zeros(N_Sites,1); %uncertainty for avg CQ
CQ_std_all = zeros(N_Sites,1); %standard deviation of CQs
CQ_SE_all = zeros(N_Sites,1); %standard error of CQs
CQalt_bar_all = zeros(N_Sites,1); %mean CQalt
CQalt_sigmaavg_all = zeros(N_Sites,1); %uncertainty for avg CQalt
CQalt_std_all = zeros(N_Sites,1); %standard deviation of CQalts
CQalt_SE_all = zeros(N_Sites,1); %standard error of CQalts

%go through Sites
for i = 1:N_Sites
    %get values for fitting - stress
    ind_fit = find(tauex_bin_min_all{i}>=N_sigma_tauit*tauit_sigma_linearfit_all(i)); %min tau must exceed threshold by at least 2 sigma
    tauex_fit_all{i} = tauex_bin_avg_all{i}(ind_fit);
    sigma_tauex_fit_all{i} = tauex_bin_sigma_all{i}(ind_fit);
    tauratio_fit_all{i} = tauratio_bin_avg_all{i}(ind_fit);
    sigma_tauratio_fit_all{i} = tauratio_bin_sigma_all{i}(ind_fit);
    
    %get values for fitting - flux
    Q_fit_all{i} = Q_bin_avg_all{i}(ind_fit);
    sigma_Q_fit_all{i} = Q_bin_sigma_all{i}(ind_fit);
    CQ_fit_all{i} = CQ_bin_avg_all{i}(ind_fit);
    sigma_CQ_fit_all{i} = CQ_bin_sigma_all{i}(ind_fit);
    CQalt_fit_all{i} = CQalt_bin_avg_all{i}(ind_fit);
    sigma_CQalt_fit_all{i} = CQalt_bin_sigma_all{i}(ind_fit);
    
    %computations for flux
    Qtauexratio_bar_all(i) = mean(Q_fit_all{i}./tauex_fit_all{i}); %mean Q/tauex ratio
    CQ_bar_all(i) = mean(CQ_fit_all{i}); %mean CQ
    [~, sigma_CQ_bar] = MeanUncertainty(CQ_fit_all{i}, sigma_CQ_fit_all{i}); %uncertainty in mean CQ
    CQ_sigmaavg_all(i) = sigma_CQ_bar; %uncertainty in mean CQ
    CQ_std_all(i) = std(CQ_fit_all{i}); %standard deviation of CQs
    CQ_SE_all(i) = CQ_std_all(i)/sqrt(length(CQ_fit_all{i})); %standard error of CQs
    CQalt_bar_all(i) = mean(CQalt_fit_all{i}); %mean CQalt
    [~, sigma_CQalt_bar] = MeanUncertainty(CQalt_fit_all{i}, sigma_CQalt_fit_all{i}); %uncertainty in mean CQalt
    CQalt_sigmaavg_all(i) = sigma_CQalt_bar; %uncertainty in mean CQalt
    CQalt_std_all(i) = std(CQalt_fit_all{i}); %standard deviation of CQalts
    CQalt_SE_all(i) = CQalt_std_all(i)/sqrt(length(CQalt_fit_all{i})); %standard error of CQalts
end

%%
%%%%%%%%%%%%%%%%%
% PRIMARY PLOTS %
%%%%%%%%%%%%%%%%%

%% PLOT ZQ VERSUS UST
legend_items = [SiteNames;LitNames];

%PANEL A - dimensional
figure(1); subplot(2,7,1:7); hold on;
%plot binned Field data
for i = 1:N_Sites
    plot(ust_zqustfit_all{i},zq_zqustfit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);
end
%plot Literature data
for i = 1:N_Lit
    plot(ust_Lit{i},zq_Lit{i},Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_Lit);
end

%plot error bars separately
for i = 1:N_Sites %plot binned Field data
    for j = 1:length(ust_zqustfit_all{i})
        plot(ones(2,1)*ust_zqustfit_all{i}(j), zq_zqustfit_all{i}(j)+[-1 1]*zq_sigma_zqustfit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %y error-bars
        %plot(ust_zqustfit_all{i}(j)+[-1 1]*ust_sigma_zqustfit_all{i}(j), ones(2,1)*zq_zqustfit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %x error-bars
    end
end
for i = 1:N_Sites %plot Literature data
    for j = 1:length(ust_Lit{i})
        plot(ones(2,1)*ust_Lit{i}(j), zq_Lit{i}(j)+[-1 1]*sigma_zq_Lit{i}(j),'Color',Colors_Lit{i},'LineWidth',LineWidth_Field); %y error-bars
    end
end

%organize plot
xlim([0.25 0.6]);
ylim([0 0.15]);
legend(legend_items,'Location','EastOutside');
text(0.255, 0.14,'(a)','FontSize',PlotFont);
set(gca,'XMinorTick','On','XScale','log','YMinorTick','On','Box','On');
xlabel('\textbf{Shear velocity, $$u_{*}$$ (ms$$^{-1}$$)}','Interpreter','Latex');
ylabel('\textbf{Saltation height, $$z_q$$ (m)}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%PANEL B - fit values
subplot(2,7,8:9); hold on;

%plot Field data
for i = 1:N_Sites
    plot(d50_Site(i),slope_zqustfit_all(i),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field*1.5); %values
    plot(ones(2,1)*d50_Site(i),slope_zqustfit_all(i)+[-1 1]*sigma_slope_zqustfit_all(i),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %error bars
end
%plot Literature data
for i = 1:2; %neglect Farrell (2012) data
    plot(d50_Lit(i),slope_zqust_Lit_all(i),Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_Lit*1.5); %values
    plot(ones(2,1)*d50_Lit(i),slope_zqust_Lit_all(i)+[-1 1]*sigma_slope_zqust_Lit_all(i),'Color',Colors_Lit{i},'LineWidth',LineWidth_Lit); %error bars
end

%plot 0 line
plot([0.2 0.6],[0 0],'k--','LineWidth',1);

%organize plot
xlim([0.2 0.6]);
ylim([-0.2 0.2]);
text(0.205, 0.17,'(b)','FontSize',PlotFont);
set(gca,'XMinorTick','On','XScale','log','YMinorTick','On','Box','On');
xlabel('\textbf{Particle diameter, $$d_{50}$$ (mm)}','Interpreter','Latex');
ylabel('\textbf{$$z_q$$-$$u_{*}$$ fit slope, $$b$$ (m/ms$$^{-1}$$)}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%PANEL B2 - fit values for Farrell et al. (2012)
subplot(2,7,10); hold on;

%plot best fit values
plot(d50_Lit(3),slope_zqust_Lit_all(3),Markers_Lit{3},'Color',Colors_Lit{3},'MarkerSize',MarkerSize_Lit*1.5); %values
plot(ones(2,1)*d50_Lit(3),slope_zqust_Lit_all(3)+[-1 1]*sigma_slope_zqust_Lit_all(3),'Color',Colors_Lit{3},'LineWidth',LineWidth_Lit); %error bars

%plot 0 line
plot([-1 1],[0 0],'k--','LineWidth',1);

%organize plot
xlim([-1 1]);
ylim([-0.2 0.2]);
set(gca,'XTick',[],'XTickLabel',{''},'YTickLabel',{''},'YMinorTick','On','Box','On');
set(gca,'FontSize',PlotFont);

%PANEL C - dimensionless saltation heights
subplot(2,7,12:14); hold on;

%plot dimensionless heights versus d50 for Field data
for i = 1:N_Sites
    plot(d50_Site(i), zqnorm_bar_all(i), Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field*1.5); %values
    plot(ones(2,1)*d50_Site(i), zqnorm_bar_all(i)+[-1 1]*zqnorm_std_all(i),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %error bars
end
%plot dimensionless heights versus d50 for lit data
for i = 1:2 %only plot first two entries, ignorning Farrell (2012) with no d50
    plot(d50_Lit(i), zqnorm_bar_Lit_all(i), Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_Lit*1.5); %values
    plot(ones(2,1)*d50_Lit(i), zqnorm_bar_Lit_all(i)+[-1 1]*zqnorm_std_Lit_all(i),'Color',Colors_Lit{i},'MarkerSize',MarkerSize_Lit*1.5); %error bars
end
%organize plot
xlim([0.2 0.6]);
ylim([0 250]);
text(0.203, 235,'(c)','FontSize',PlotFont);
set(gca,'XMinorTick','On','XScale','log','YMinorTick','On','Box','On');
xlabel('\textbf{Particle diameter, $$d_{50}$$ (mm)}','Interpreter','Latex');
ylabel('\textbf{Dimensionless salt. ht., $$\langle z_q \rangle/d_{50}$$}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8.5 8]);
print([folder_Plots,'zq_ust.png'],'-dpng');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 18 12.7]);
print([folder_Plots,'zq_ust.eps'],'-depsc');


%% PLOT Q VERSUS TAUEX BINNED DATA
figure(2); clf; hold on;

%plot binned values
for i = 1:N_Sites
    plot(tauex_fit_all{i},Q_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);
end

%plot error bars
for i = 1:N_Sites
    for j = 1:length(tauex_fit_all{i})
        plot(ones(2,1)*tauex_fit_all{i}(j),Q_fit_all{i}(j)+[-1 1]*sigma_Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %y error
        plot(tauex_fit_all{i}(j)+[1 -1]*sigma_tauex_fit_all{i}(j),ones(2,1)*Q_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %x error
    end
end

%plot fit values
tauex_fit = [0 0.35];
for i = 1:N_Sites
    plot(tauex_fit, Qtauexratio_bar_all(i)*tauex_fit,'Color',Colors_Field{i});
end

xlim([0 0.35]);
ylim([0 65]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Excess shear stress, $$\tau_{ex}$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)}','Interpreter','Latex');
legend(SiteNames,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'Q_tauex.png'],'-dpng');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8.7 7]);
print([folder_Plots,'Q_tauex.eps'],'-depsc');


%% PLOT CQ VERSUS TAUEX
figure(3); clf; hold on;

%plot binned values
for i = 1:N_Sites
    plot(tauratio_fit_all{i},CQ_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);
    %errorbar(tauex_fit_all{i},CQ_fit_all{i},sigma_CQ_fit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
end

%plot error bars
for i = 1:N_Sites
    for j = 1:length(tauratio_fit_all{i})
        plot(ones(2,1)*tauratio_fit_all{i}(j),CQ_fit_all{i}(j)+[-1, 1]*sigma_CQ_fit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field);
    end
end

%plot fit values
tauratio_fit = [1 4];
for i = 1:N_Sites
    plot(tauratio_fit, ones(2,1)*CQ_bar_all(i),'Color',Colors_Field{i});
end

xlim([1 4]);
ylim([0 10]);
set(gca,'XMinorTick','On','YMinorTick','On','XScale','Log','Box','On');
xlabel('\textbf{Dimensionless shear stress, $$\tau/\tau_{it}$$}','Interpreter','Latex');
ylabel('\textbf{Dimensionless saltation flux, $$\hat{C_{Q}}$$}','Interpreter','Latex');
legend(SiteNames,'Location','SouthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gca,'LooseInset',get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 5]);
print([folder_Plots,'CQ_tauratio.png'],'-dpng');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8.7 6.4]);
print([folder_Plots,'CQ_tauratio.eps'],'-depsc');

%%
%%%%%%%%%%%%%%%%%%%%%%%
% SUPPLEMENTARY PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%

%% PLOT THETA MINUS MEAN THETA VERSUS TAU
figure(13); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(tauRe_all{i}(Q_all{i}>0),theta_all{i}(Q_all{i}>0)-mean(theta_all{i}(Q_all{i}>0)),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
    plot(tauRe_all{i}(Q_all{i}==0),theta_all{i}(Q_all{i}==0)-mean(theta_all{i}(Q_all{i}>0)),Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
    legend('transport','no transport');
    xlim([0 0.4]);
    ylim([-90 90]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{30 min. shear stress, $$\tilde{\tau}$$ (Pa)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Wind angle, $$\tilde{\theta} (^{\circ})$$}','Interpreter','Latex');
        text(0.01, 83,'(a)','FontSize',PlotFont);
    elseif i==2
        text(0.01, 83,'(b)','FontSize',PlotFont);
    elseif i==3
        text(0.01, 83,'(c)','FontSize',PlotFont);
    end
    title(SiteNames{i});
    set(gca,'FontSize',13);
end

%calculate adjusted theta
theta_adjusted_all = cell(N_Sites,1);
for i = 1:N_Sites
    theta_adjusted_all{i}=theta_all{i}-mean(theta_all{i}(Q_all{i}>0));
end

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'theta_tau.png'],'-dpng');


%% PLOT z/L VERSUS TAU
figure(14); clf;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(tauRe_all{i}(Q_all{i}>0),zL_all{i}(Q_all{i}>0),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
    plot(tauRe_all{i}(Q_all{i}==0),zL_all{i}(Q_all{i}==0),Markers_Lit{i},'Color',Colors_Lit{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
    legend('transport','no transport','Location','SouthEast');
    xlim([0 0.4]);
    ylim([-0.5 0]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{30 min. shear stress, $$\tilde{\tau}$$ (Pa)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Stability parameter, $$z/L$$}','Interpreter','Latex');
        text(0.01, -0.015,'(a)','FontSize',PlotFont);
    elseif i==2
        text(0.01, -0.015,'(b)','FontSize',PlotFont);
    elseif i==3
        text(0.01, -0.015,'(c)','FontSize',PlotFont);
    end
    title(SiteNames{i});
    set(gca,'FontSize',13);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'zL_tau.png'],'-dpng');

%% PLOT RAW ZQ DATA
figure(17);
for i = 1:N_Sites
    %zq plot
    subplot(1,N_Sites,i); hold on;
    
    %plot data
    plot(ustRe_all{i},zq_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);

    %plot error bars
    for j = 1:length(ustRe_all{i})
        plot(ones(2,1)*ustRe_all{i}(j),zq_all{i}(j)+[-1 1]*sigma_zq_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %y-error
        plot(ustRe_all{i}(j)+[-1 1]*sigma_ustRe_all{i}(j),ones(2,1)*zq_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %x-error
    end
    
    xlim([0.25 0.6]);
    ylim([0 0.14]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','log');
    xlabel('\textbf{30 min. shear vel., $$\tilde{u}_{*}$$ (ms$$^{-1}$$)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{30 min. saltation height, $$\tilde{z}_q$$ (m)}','Interpreter','Latex');
        text(0.26, 0.135,'(a)','FontSize',PlotFont);
    elseif i==2
        text(0.26, 0.135,'(b)','FontSize',PlotFont);
    elseif i==3
        text(0.26, 0.135,'(c)','FontSize',PlotFont);
    end
    title(SiteNames{i});
    set(gca,'FontSize',13);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'zq_ust_raw.png'],'-dpng');

%% PLOT ALL Q VERSUS ALL TAU
figure(18);
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;

    %plot data
    plot(tauRe_all{i},Q_all{i},Markers_Field{i},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{i});

    %plot error bars
    for j = 1:length(tauRe_all{i})
        plot(ones(2,1)*tauRe_all{i}(j),Q_all{i}(j)+[-1 1]*sigma_Q_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %y error
        plot(tauRe_all{i}(j)+[1 -1]*sigma_tauRe_all{i}(j),ones(2,1)*Q_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %x error
    end
    xlim([0 0.45]);
    ylim([0 65]);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{30 min. wind stress, $$\tilde{\tau}$$ (Pa)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{30 min. saltation flux, $$\tilde{Q}$$ (gm$$^{-1}$$s$$^{-1}$$)}','Interpreter','Latex');
        text(0.01, 62,'(a)','FontSize',PlotFont);
    elseif i==2
        text(0.01, 62,'(b)','FontSize',PlotFont);
    elseif i==3
        text(0.01, 62,'(c)','FontSize',PlotFont);
    end
    title(SiteNames{i});
    set(gca,'FontSize',13);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf, 'PaperPosition',[0 0 11 5]);
print([folder_Plots,'Q_tau_raw.png'],'-dpng');

%% PLOT TYPICAL STANDARD DEVIATIONS FOR Q's IN BINS

figure(19); hold on;
%plot standard deviations
for i = 1:N_Sites
    ind_full = find(bin_N_all{i}>=3); %indices for full bins
    ind_transport = find(fQ_bin_max_all{i}>=fQ_Qtau_fit_min); %indices for bins with transport
    ind_stdmed_Q = intersect(ind_full,ind_transport); %indices for computing minimum standard deviation
    plot(Q_bin_avg_all{i}(ind_stdmed_Q),Q_bin_std_all{i}(ind_stdmed_Q),Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
end

%plot median standard deviations
for i = 1:N_Sites
    plot(xlim,[Q_bin_stdmed_all(i) Q_bin_stdmed_all(i)],'Color',Colors_Field{i},'LineWidth',LineWidth_Field);
end

set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Saltation mass flux, $$Q_i$$ (gm$$^{-1}$$s$$^{-1}$$)}','Interpreter','Latex');
ylabel('\textbf{Salt. mass flux std. dev., $$SD_{Q_i}$$ (gm$$^{-1}$$s$$^{-1}$$)}','Interpreter','Latex');
legend(SiteNames,'Location','NorthEast');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition',[0 0 8 5]);
print([folder_Plots,'Q_std.png'],'-dpng');

%% PLOT Q VERSUS TAU
figure(20); hold on;
legend_items = cell(N_Sites*3,1);

%plot binned values
for i = 1:N_Sites
    plot(tau_Qtaufit_all{i},Q_Qtaufit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);
    plot(tau_linearfit_plot{i},Q_linearfit_plot{i},'Color',Colors_Field{i}); %plot linear fit
    plot(tau_threehalvesfit_plot{i},Q_threehalvesfit_plot{i},'--','Color',Colors_Field{i}); %plot 3/2 fit
    legend_items{3*i-2} = SiteNames{i}; %add to list of legend items
    legend_items{3*i-1} = 'linear fit';
    legend_items{3*i} = 'nonlinear 3/2 fit';

    %print out fit values
    linear_output = ['linear, \chi^2_{\nu} = ',...
        num2str(Chi2_linearfit_all(i)/df_linearfit_all(i),'%.2f')]
    threehalves_output = ['3/2, \chi^2_{\nu} = ',...
        num2str(Chi2_threehalvesfit_all(i)/df_threehalvesfit_all(i),'%.2f')]
end

%plot error bars
for i = 1:N_Sites
    for j = 1:length(tau_Qtaufit_all{i})
        plot(ones(2,1)*tau_Qtaufit_all{i}(j),Q_Qtaufit_all{i}(j)+[-1 1]*Q_sigma_Qtaufit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %y error
        plot(tau_Qtaufit_all{i}(j)+[1 -1]*tau_sigma_Qtaufit_all{i}(j),ones(2,1)*Q_Qtaufit_all{i}(j),'Color',Colors_Field{i},'LineWidth',LineWidth_Field); %x error
    end
end

xlim([0 0.45]);
ylim([0 65]);
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('\textbf{Shear stress, $$\tau$$ (Pa)}','Interpreter','Latex');
ylabel('\textbf{Saltation mass flux, $$Q$$ (gm$$^{-1}$$s$$^{-1}$$)}','Interpreter','Latex');
legend(legend_items,'Location','NorthWest');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition',[0 0 8 6]);
print([folder_Plots,'Q_tau.png'],'-dpng');

%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%

save(path_SaveData,'*all');

%%%%%%%%%%%%%%%%%%%%
% DIAGNOSTIC PLOTS %
%%%%%%%%%%%%%%%%%%%%

% %% PLOT TAU RELATIVE UNCERTAINTIES VERSUS TAU
% figure; hold on;
% for i = 1:N_Sites
%     plot(tauRe_all{i},sigma_tauRe_all{i}./tauRe_all{i},Markers_Field{i},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{i});
% end
% xlim([0 0.35]);
% ylim([0 0.12]);
% legend(Sites,'Location','NorthEast');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
% ylabel('\textbf{$$\sigma_{\tau}/\tau$$}','Interpreter','Latex');
% set(gca,'FontSize',PlotFont);
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf, 'PaperPosition',[0 0 8 6]);
% print([folder_Plots,'sigma_tau.png'],'-dpng');

% %% PLOT RELATIVE UNCERTAINTY IN Q VERSUS TOTAL UNCERTAINTY
% figure; hold on;
% for i = 1:N_Sites
%     plot(Q_all{i},sigma_Q_all{i}./Q_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);
% end
% xlim([0 55]);
% ylim([0 0.3]);
% xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
% ylabel('\textbf{Wenglor relative uncertainty, $$\sigma_Q/Q$$}','Interpreter','Latex');
% set(gca,'FontSize',PlotFont);
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
% print([folder_Plots,'sigma_Q_Wenglor_profilefit.png'],'-dpng');
% 
% %% PLOT RELATIVE UNCERTAINTY IN ZQ VERSUS TOTAL UNCERTAINTY
% figure; hold on;
% for i = 1:N_Sites
%     plot(Q_all{i},sigma_zq_all{i}./zq_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field);
% end
% xlim([0 55]);
% ylim([0 0.15]);
% xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
% ylabel('\textbf{Wenglor relative uncertainty, $$\sigma_{zq}/{z_q}$$}','Interpreter','Latex');
% set(gca,'FontSize',PlotFont);
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
% print([folder_Plots,'sigma_zq_Wenglor_profilefit.png'],'-dpng');


% %% PLOT VARIANTS ON TAU BIN UNCERTAINTY
% figure;
% for i = 1:N_Sites
%     %full plot
%     subplot(2,3,i); hold on;
%     plot(tau_bin_avg_all{i},tau_bin_sigma_all{i}./tau_bin_avg_all{i},Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(tau_bin_avg_all{i},tau_bin_sigmaavg_all{i}./tau_bin_avg_all{i},Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     plot(tau_bin_avg_all{i},tau_bin_SE_all{i}./tau_bin_avg_all{i},Markers_Field{3},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{3});
%     legend('total uncertainty','random error','standard error','Location','NorthEast');
%     title(Sites{i});
%     xlim([0 0.4]);
%     ylim([0 0.4]);
%     plot(xlim,[0.05 0.05],'k-.'); %plot line showing cutoff for lower plot
%     xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
%     ylabel('\textbf{$$\sigma_{\tau}/\tau$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
%     
%     %replot with zoomed in range
%     subplot(2,3,i+3); hold on;
%     plot(tau_bin_avg_all{i},tau_bin_sigma_all{i}./tau_bin_avg_all{i},Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(tau_bin_avg_all{i},tau_bin_sigmaavg_all{i}./tau_bin_avg_all{i},Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     plot(tau_bin_avg_all{i},tau_bin_SE_all{i}./tau_bin_avg_all{i},Markers_Field{3},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{3});
%     xlim([0 0.4]);
%     ylim([0 0.05]);
%     xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
%     ylabel('\textbf{$$\sigma_{\tau}/\tau$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
% print([folder_Plots,'tau_bin_uncertainty.png'],'-dpng');
% 
% 
% %% PLOT VARIANTS ON Q BIN UNCERTAINTY
% figure;
% for i = 1:N_Sites
%     %full plot
%     subplot(2,3,i); hold on;
%     plot(Q_bin_avg_all{i},Q_bin_sigma_all{i}./Q_bin_avg_all{i},Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(Q_bin_avg_all{i},Q_bin_sigmaavg_all{i}./Q_bin_avg_all{i},Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     plot(Q_bin_avg_all{i},Q_bin_SE_all{i}./Q_bin_avg_all{i},Markers_Field{3},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{3});
%     legend('total uncertainty','random error','standard error','Location','NorthEast');
%     title(Sites{i});
%     xlim([0 50]);
%     ylim([0 1]);
%     plot(xlim,[0.3 0.3],'k-.'); %plot line showing cutoff for lower plot
%     xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
%     ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
%     
%     %replot with zoomed in range
%     subplot(2,3,i+3); hold on;
%     plot(Q_bin_avg_all{i},Q_bin_sigma_all{i}./Q_bin_avg_all{i},Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(Q_bin_avg_all{i},Q_bin_sigmaavg_all{i}./Q_bin_avg_all{i},Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     plot(Q_bin_avg_all{i},Q_bin_SE_all{i}./Q_bin_avg_all{i},Markers_Field{3},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{3});
%     xlim([0 50]);
%     ylim([0 0.3]);
%     xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
%     ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
% print([folder_Plots,'Q_bin_uncertainty.png'],'-dpng');


% %% PLOT CONTRIBUTIONS TO Q BIN UNCERTAINTY
% figure;
% for i = 1:N_Sites
%     %for linear fit
%     subplot(2,3,i); hold on;
%     plot(Q_Qtaufit_all{i},Q_sigmatotal_linearfit_all{i}./Q_Qtaufit_all{i},Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(Q_Qtaufit_all{i},Q_sigma_Qtaufit_all{i}./Q_Qtaufit_all{i},Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     plot(Q_Qtaufit_all{i},Q_sigmatau_linearfit_all{i}./Q_Qtaufit_all{i},Markers_Lit{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{3});
%     legend('total uncertainty','Q contribution','\tau contribution','Location','NorthEast');
%     title([Sites{i}, ' - ','Linear Fit']);
%     xlim([0 50]);
%     ylim([0 0.8]);
%     xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
%     ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
%     
%     %for nonlinear fit
%     subplot(2,3,i+3); hold on;
%     plot(Q_Qtaufit_all{i},Q_sigmatotal_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(Q_Qtaufit_all{i},Q_sigma_Qtaufit_all{i}./Q_Qtaufit_all{i},Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     plot(Q_Qtaufit_all{i},Q_sigmatau_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_Lit{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{3});
%     legend('total uncertainty','Q contribution','\tau contribution','\tau contr. (alt calc)','Location','NorthEast');
%     title([Sites{i}, ' - ','Nonlinear Fit']);
%     xlim([0 50]);
%     ylim([0 0.8]);
%     xlabel('\textbf{$$Q$$ (Pa)}','Interpreter','Latex');
%     ylabel('\textbf{$$\sigma_{Q}/Q$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 8]);
% print([folder_Plots,'Q_uncertainty_contributions.png'],'-dpng');


% %% PLOT CONTRIBUTIONS TO CHI2 UNCERTAINTY
% figure;
% for i = 1:N_Sites
%     %for linear fit
%     subplot(1,3,i); hold on;
%     plot(tau_Qtaufit_all{i},Chi2_contributions_linearfit_all{i}/df_linearfit_all(i),Markers_Field{1},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{1});
%     plot(tau_Qtaufit_all{i},Chi2_contributions_nonlinearfit_all{i}/df_nonlinearfit_all(i),Markers_Field{2},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{2});
%     legend('Linear fit','Nonlinear fit','Location','NorthEast');
%     title(Sites{i});
%     xlim([0 0.4]);
%     ylim([0 0.5]);
%     xlabel('\textbf{$$\tau$$}','Interpreter','Latex');
%     ylabel('\textbf{$$\chi^2_{\nu}$$}','Interpreter','Latex');
%     set(gca,'FontSize',PlotFont,'XMinorTick','On','YMinorTick','On','Box','On');
% end
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
% print([folder_Plots,'Chi2_contributions_Qtaufit.png'],'-dpng');
% 
% 
% %% PLOT CHI2 FOR WENGLOR PROFILE VERSUS Q
% figure;
% for i = 1:N_Sites
%     subplot(1,N_Sites,i); hold on;
%     plot(Q_all{i},Chi2_Qfit_all{i}./df_Qfit_all{i},Markers_Field{i},'Color',Colors_Field{i},'MarkerSize',MarkerSize_Field,'LineWidth',LineWidth_Field);
%     xlim([0 50]);
%     ylim([1e-3 1e2]);
%     set(gca,'yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
%     xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-1}$$ s$$^{-1}$$)}','Interpreter','Latex');
%     if i==1
%         ylabel('\textbf{Quality of fit for Wenglor profile, $$\chi^2_{\nu}$$}','Interpreter','Latex');
%     end
%     title(Sites{i});
%     set(gca,'FontSize',PlotFont);
% end
% 
% %print plot for draft
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
% print([folder_Plots,'Chi2_qz_profilefit.png'],'-dpng');
% 
% 
% %% PLOT CQN VERSUS Z
% ID_unique = unique(horzcat(W_ID_all{3}{Q_all{3}>0})); %get Wenglor IDs from Oceano
% N_ID = length(ID_unique);
% Cqn_plot_markers = {'r+','go','b*','cx','ms','kd','r^','gv','b>','c>','mp','kh','b+','co'};
% 
% figure;
% for i = 1:N_Sites
%     %identify time intervals with transport
%     ind_trans = find(Q_all{i}>0);
%     N_trans = length(ind_trans);
%     
%     %initialize list of calibration coefficients and heights
%     Cqn_ID = cell(N_ID,1);
%     zW_ID = cell(N_ID,1);
%     
%     %sort calibration values by Wenglor ID
%     for j = 1:N_ID
%         Cqn_list = []; %initialize list of Cqn
%         zW_list = []; %initialize list of zW
%         for k = 1:N_trans; %go through each interval to get values for list
%             zW_interval = zW_all{i}{ind_trans(k)};
%             Cqn_interval = Cqn_all{i}{ind_trans(k)};
%             ID_interval = W_ID_all{i}{ind_trans(k)};
%             ind_ID = find(strcmp(ID_interval,ID_unique(j)));
%             if length(ind_ID)==1
%                 Cqn_list = [Cqn_list; Cqn_interval(ind_ID)];
%                 zW_list = [zW_list; zW_interval(ind_ID)];
%             end
%         end
%         Cqn_ID{j} = Cqn_list;
%         zW_ID{j} = zW_list;
%     end
%     
%     %plot
%     if i~=3
%         subplot(1,4,i); hold on;
%     else
%         subplot(1,4,3:4); hold on;
%     end
%     for j = 1:N_ID
%         plot(zW_ID{j},Cqn_ID{j},Cqn_plot_markers{j});
%     end
%     if i==3
%         legend(ID_unique,'Location','EastOutside');
%     end
%     
%     title(Sites{i});
%     ylim([0.1 100]);
%     set(gca,'yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
%     xlabel('\textbf{Wenglor ht., $$z$$ (m)}','Interpreter','Latex');
%     if i==1
%         ylabel('\textbf{Calibration factor, $$C_{qn}$$ (g m$$^{-1}$$ s$$^{-1}$$ count$$^{-1}$$)}','Interpreter','Latex');
%     end
%     set(gca,'FontSize',PlotFont);
% end
% 
% %print plot
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 5]);
% print([folder_Plots,'Cqn_z.png'],'-dpng');
% 
%
% %% PLOT NORMALIZED RESIDUALS FOR FIT
% 
% figure;
% subplot(2,1,1); hold on; %residuals for linear fit
% for i = 1:N_Sites
%     errorbar(tau_Qtaufit_all{i},Q_residuals_linearfit_all{i}./Q_Qtaufit_all{i},Q_sigmatotal_linearfit_all{i}./Q_Qtaufit_all{i},Markers_Field{i},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{i}); % include uncertainties from tau
% end
% xlim([0 0.4]);
% ylim([-1 1]);
% plot(xlim,[0 0]); %plot zero line
% title('Linear fit');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
% ylabel('\textbf{$$(Q - Q_{pred})/Q$$}','Interpreter','Latex');
% legend(Sites,'Location','NorthEast');
% set(gca,'FontSize',PlotFont);
%     
% subplot(2,1,2); hold on; %residuals for nonlinear fit
% for i = 1:N_Sites
%     errorbar(tau_Qtaufit_all{i},Q_residuals_nonlinearfit_all{i}./Q_Qtaufit_all{i},Q_sigmatotal_nonlinearfit_all{i}./Q_Qtaufit_all{i},Markers_Field{i},'MarkerSize',MarkerSize_Field,'Color',Colors_Field{i}); % include uncertainties from tau
% end
% xlim([0 0.4]);
% ylim([-1 1]);
% plot(xlim,[0 0]); %plot zero line
% title('Nonlinear fit');
% set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
% xlabel('\textbf{$$\tau$$ (Pa)}','Interpreter','Latex');
% ylabel('\textbf{$$(Q - Q_{pred})/Q$$}','Interpreter','Latex');
% legend(Sites,'Location','NorthEast');
% set(gca,'FontSize',PlotFont);
% 
% set(gca, 'LooseInset', get(gca,'TightInset'));
% set(gcf, 'PaperPosition',[0 0 8 10]);
% print([folder_Plots,'Qresiduals_tau.png'],'-dpng');

%%
