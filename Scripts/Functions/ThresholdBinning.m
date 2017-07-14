%% FUNCTION TO ANALYZE THRESHOLD USING TFEM METHOD

function [fQ_bin_avg,fQ_bin_SE,...
    uth_bin_avg,uth_bin_SE,...
    ustth_bin_avg,ustth_bin_SE,...
    tauth_bin_avg,tauth_bin_SE,...
    tauft, sigma_tauft,...
    tauit, sigma_tauit,...
    ustitftratio, sigma_ustitftratio,...
    Chi2nu, sigma_b] = ...
ThresholdBinning(rho,z0,sigma_z0,...
    Deltat_all,deltat_all,...
    Deltat_analysis, deltat_analysis,...
    Sites,Site_analysis,...
    fQ_all, uth_all, zU_all,...
    theta_all,theta_max,zL_all,zL_max,...
    starttime_all,timeofday_all,...
    starthour_analysis,endhour_analysis,...
    startdate_analysis,enddate_analysis)

%set time range if none given 
if nargin<=18
    starthour_analysis = 0;
    endhour_analysis = 24;
end

%set date range if none given
if nargin<=20
    startdate_analysis = datetime(0,0,0);
    enddate_analysis = datetime(10000,0,0);
end

%% parameter values
kappa = 0.4; %von Karman parameter

%% binning info
fQ_min = 0.05;
fQ_max = 0.95;
fQ_bin_minrange = 0.1;
fQ_bin_maxrange = 0.2;
fQ_bin_N_min = 3;

%% get index of site for analysis
ind_Site = find(strcmp(Site_analysis,Sites));

%% get indices of timescales for analysis
ind_Deltat = find(Deltat_all==Deltat_analysis); %get measurement interval
ind_deltat = find(deltat_all==deltat_analysis); %get sampling interval

%get values
starttime = starttime_all{ind_Site}{ind_Deltat,ind_deltat};
timeofday = timeofday_all{ind_Site}{ind_Deltat,ind_deltat};
fQ = fQ_all{ind_Site}{ind_Deltat,ind_deltat};
uth = uth_all{ind_Site}{ind_Deltat,ind_deltat};
zU = zU_all{ind_Site}{ind_Deltat,ind_deltat};
theta = theta_all{ind_Site}{ind_Deltat};
zL = zL_all{ind_Site}{ind_Deltat,ind_deltat};

%get indices for binning
ind_fQ = intersect(find(fQ>=fQ_min),find(fQ<=fQ_max)); %get indices for fQ range
ind_time = intersect(find(timeofday>=starthour_analysis),find(timeofday<=endhour_analysis)); %get indices of time of day range
ind_date = intersect(find(days(starttime-startdate_analysis)>0),find(days(starttime-enddate_analysis)<0)); %get indices of date range if necessary
ind_datetime = intersect(ind_time,ind_date);
ind_theta = find(abs(theta)<=theta_max);
ind_zL = find(abs(zL)<=zL_max);
ind_wind = intersect(ind_theta,ind_zL);
ind_binning = intersect(ind_fQ,intersect(ind_datetime,ind_wind));

%keep only values for binning
fQ_binning = fQ(ind_binning);
uth_binning = uth(ind_binning);
zU_binning = zU(ind_binning);

%% PERFORM BINNING - TOTAL
if ~isempty(fQ_binning)

    %get binned values for fQ
    [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
    N_fQ_bins = length(fQ_bin_min);

    %get binned values for uth
    uth_bin_avg = zeros(N_fQ_bins,1);
    uth_bin_SE = zeros(N_fQ_bins,1);
    for k=1:N_fQ_bins   
        fQ_bin_ind = find(fQ_binning>=fQ_bin_min(k)&fQ_binning<=fQ_bin_max(k));
        uth_bin_avg(k) = mean(uth_binning(fQ_bin_ind));
        uth_bin_SE(k) = std(uth_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
    end
    
    %get binned values for zU
    zU_bin_avg = zeros(N_fQ_bins,1);
    zU_bin_SE = zeros(N_fQ_bins,1);
    for k=1:N_fQ_bins   
        fQ_bin_ind = find(fQ_binning>=fQ_bin_min(k)&fQ_binning<=fQ_bin_max(k));
        zU_bin_avg(k) = mean(zU_binning(fQ_bin_ind));
        zU_bin_SE(k) = std(zU_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
    end
end

%% CONVERT TO SHEAR VELOCITY AND STRESS
ustth_bin_avg = (kappa*uth_bin_avg)./log(zU_bin_avg./z0); %threshold shear velocity
%ustth_bin_SE = (uth_bin_SE./uth_bin_avg).*ustth_bin_avg; %threshold shear velocity - uncertainty
ustth_bin_SE = (kappa./log(zU_bin_avg./z0)).*sqrt(uth_bin_SE.^2 + (log(sigma_z0).^2).*(ustth_bin_avg./log(zU_bin_avg./z0).^2)); %threshold shear velocity - uncertainty - including z0 uncertainty
tauth_bin_avg = rho*ustth_bin_avg.^2; %threshold shear stress
tauth_bin_SE = 2*rho*ustth_bin_avg.*ustth_bin_SE; %threshold shear stress - uncertainty

%get best fit values
[a, b, sigma_a, sigma_b, yfit, ~, sigma2_ab, ~, ~] = linearfit(fQ_bin_avg, tauth_bin_avg, tauth_bin_SE);
tauft = a;
sigma_tauft = sigma_a;
tauit = a+b;
sigma_tauit = sqrt(sigma_a^2+sigma_b^2+2*sigma2_ab);

%estimate chi2
tauth_residuals = yfit - tauth_bin_avg; %get residuals for fit
Chi2 = sum((tauth_residuals./tauth_bin_SE).^2); %compute total Chi2
df = length(tauth_residuals)-2; %compute degrees of freedom
Chi2nu = Chi2/df;

%get ratio of ustit / ustft
ustitftratio = sqrt(tauit/tauft);
sigma_ustitftratio = (1/2)*sqrt((sigma_tauft^2/(tauft*tauit))+...
    (sigma_tauit^2*tauft/tauit^3));