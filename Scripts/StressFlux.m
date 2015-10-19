%% ANALYZE FLUX VERSUS STRESS RELATIONSHIP

%initialize
clearvars;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots

%load data
%load(strcat(folder_ProcessedData,'StressFluxContinuousWindows')); %stress/flux continuous windows
load(strcat(folder_ProcessedData,'StressFluxWindows')); %stress/flux all windows
load(strcat(folder_ProcessedData,'C_ustThr_sqrtd50')); %threshold scaling
load(strcat(folder_ProcessedData,'GreeleyNamikasData')); %literature data

%get info about sites, set markers for plotting
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%set parameter values
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

%% GET FLUX VS U* IN A VARIETY OF WAYS, WITH U* THRESHOLD COMPUTED FROM SQRT(D50) FROM FITTING PARAMETER

%initialize cell arrays of lists by site
Q_Site = cell(N_Sites,1); %absolute flux (g/m/s)
Q_Einstein_Site = cell(N_Sites,1); %flux with Einstein flux non-dimensionalization
Q_norm2_Site = cell(N_Sites,1); %Q normalized assuming squared law, has units of (m/s)^2
Q_norm3_Site = cell(N_Sites,1); %Q normalized assuming cubed law, has units of (m/s)^3
Q_nondim2_Site = cell(N_Sites,1); %Q2 normalization divided by excess stress times air density
Q_nondim3_Site = cell(N_Sites,1); %Q3 normalization divided by excess stress times shear velocity times air density
ustThr_Site = cell(N_Sites,1); %u* threshold
ust_ratio_Site = cell(N_Sites,1); %ratio of u*/u*th
ust_ex_Site = cell(N_Sites,1); %excess shear velocity, u*-u*th
tauThr_Site = cell(N_Sites,1); %stress threshold
tau_ratio_Site = cell(N_Sites,1); %ratio of tau/tau_th
tau_ex_Site = cell(N_Sites,1); %excess shear stress, tau-tau_th
tau_Shields_Site = cell(N_Sites,1); %Shields stress
ust_tau_ex_Site = cell(N_Sites,1); %product of u* and excess stress
d50_Site = cell(N_Sites,1); %median grain diameter (mm)

%initialize full combined lists
Q_all = []; %absolute flux (g/m/s)
Q_Einstein_all = []; %with Einstein flux non-dimensionalization
Q_norm2_all = []; %Q normalized assuming squared law, has units of (m/s)^2
Q_norm3_all = []; %Q normalized assuming cubed law, has units of (m/s)^3
Q_nondim2_all = []; %Q2 normalization divided by excess stress times air density
Q_nondim3_all = []; %Q3 normalization divided by excess stress times shear velocity times air density
ust_all = []; %u*
tau_all = []; %tau
ust_ratio_all = []; %ratio of u*/u*th
tau_ratio_all = []; %ratio of tau/tau_th
ust_ex_all = []; %excess shear velocity, u*-u*th
tau_ex_all = []; %excess shear stress, tau-tau_th
tau_Shields_all = []; %shields stress
ust_tau_ex_all = []; %product of u* and excess stress

%initialize binning - ust_bins
ust_bins_min = 0.3:.025:0.55;
ust_bins_max = 0.325:.025:0.575;
ust_bins_mid = mean([ust_bins_min;ust_bins_max]);
N_ust_bins = length(ust_bins_mid);
Q_ust_bin_avg = cell(N_Sites,1);
Q_ust_bin_SE = cell(N_Sites,1);

%ust_ratio bins
ust_ratio_bins_min = 1:.05:1.85;
ust_ratio_bins_max = 1.05:.05:1.9;
ust_ratio_bins_mid = mean([ust_ratio_bins_min;ust_ratio_bins_max]);
N_ust_ratio_bins = length(ust_ratio_bins_mid);
Q_Einstein_ust_ratio_bin_avg = cell(N_Sites,1);
Q_Einstein_ust_ratio_bin_SE = cell(N_Sites,1);
C_Q_nondim2_ust_ratio_bin_avg = cell(N_Sites,1);
C_Q_nondim2_ust_ratio_bin_SE = cell(N_Sites,1);
C_Q_nondim3_ust_ratio_bin_avg = cell(N_Sites,1);
C_Q_nondim3_ust_ratio_bin_SE = cell(N_Sites,1);

%ust_ex bins
ust_ex_bins_min = 0:.025:0.275;
ust_ex_bins_max = 0.025:.025:0.3;
ust_ex_bins_mid = mean([ust_ex_bins_min;ust_ex_bins_max]);
N_ust_ex_bins = length(ust_ex_bins_mid);
Q_ust_ex_bin_avg = cell(N_Sites,1);
Q_ust_ex_bin_SE = cell(N_Sites,1);
Q_Einstein_ust_ex_bin_avg = cell(N_Sites,1);
Q_Einstein_ust_ex_bin_SE = cell(N_Sites,1);

%tau_bins
tau_bins_min = 0:.025:0.375;
tau_bins_max = 0.025:.025:0.4;
tau_bins_mid = mean([tau_bins_min;tau_bins_max]);
N_tau_bins = length(tau_bins_mid);
Q_tau_bin_avg = cell(N_Sites,1);
Q_tau_bin_SE = cell(N_Sites,1);

%tau_ratio bins
tau_ratio_bins_min = 1:.5:4;
tau_ratio_bins_max = 1.5:.5:4.5;
tau_ratio_bins_mid = mean([tau_ratio_bins_min;tau_ratio_bins_max]);
N_tau_ratio_bins = length(tau_ratio_bins_mid);
Q_Einstein_tau_ratio_bin_avg = cell(N_Sites,1);
Q_Einstein_tau_ratio_bin_SE = cell(N_Sites,1);

%tau_ex bins
tau_ex_bins_min = 0:.025:0.275;
tau_ex_bins_max = 0.025:.025:0.3;
tau_ex_bins_mid = mean([tau_ex_bins_min;tau_ex_bins_max]);
N_tau_ex_bins = length(tau_ex_bins_mid);
Q_tau_ex_bin_avg = cell(N_Sites,1);
Q_tau_ex_bin_SE = cell(N_Sites,1);
Q_Einstein_tau_ex_bin_avg = cell(N_Sites,1);
Q_Einstein_tau_ex_bin_SE = cell(N_Sites,1);
Q_norm2_tau_ex_bin_avg = cell(N_Sites,1);
Q_norm2_tau_ex_bin_SE = cell(N_Sites,1);

%ust_tau_ex bins
ust_tau_ex_bins_min = 0.03:.01:0.12;
ust_tau_ex_bins_max = 0.04:.01:0.13;
ust_tau_ex_bins_mid = mean([ust_tau_ex_bins_min;ust_tau_ex_bins_max]);
N_ust_tau_ex_bins = length(ust_tau_ex_bins_mid);
Q_norm3_ust_tau_ex_bin_avg = cell(N_Sites,1);
Q_norm3_ust_tau_ex_bin_SE = cell(N_Sites,1);

%tau_Shields bins
tau_Shields_bins_min = 0.0075:.0025:0.03;
tau_Shields_bins_max = 0.01:.0025:0.0325;
tau_Shields_bins_mid = mean([tau_Shields_bins_min;tau_Shields_bins_max]);
N_tau_Shields_bins = length(tau_Shields_bins_mid);
Q_Einstein_tau_Shields_bin_avg = cell(N_Sites,1);
Q_Einstein_tau_Shields_bin_SE = cell(N_Sites,1);

for i = 1:N_Sites
%     %get values from lists - continuous windows only
%     Q_Site{i} = Q_continuous{i};
%     d50_Site{i} = d50_continuous{i};
%     ust_Site{i} = ustRe_cal_continuous{i};
%     tau_Site{i} = tauRe_cal_continuous{i};
    
    %get values from lists - all windows
    Q_Site{i} = Q_list{i};
    d50_Site{i} = d50_list{i};
    ust_Site{i} = ustRe_cal_list{i};
    tau_Site{i} = tauRe_cal_list{i};

    %compute thresholds
    ustThr_Site{i} = C_ustThr_sqrtd50*sqrt(d50_Site{i});
    tauThr_Site{i} = rho_a*ustThr_Site{i}.^2;
    
    %perform non-dimensionalization and excess calculations
    ust_ratio_Site{i} = ust_Site{i}./ustThr_Site{i};
    ust_ex_Site{i} = ust_Site{i} - ustThr_Site{i};
    tau_ratio_Site{i} = ust_ratio_Site{i}.^2;
    tau_ex_Site{i} = tau_Site{i} - tauThr_Site{i};
    tau_Shields_Site{i} = tau_Site{i}./(rho_s*g*(d50_Site{i}/1000));
    ust_tau_ex_Site{i} = ust_Site{i}.*tau_ex_Site{i}; %product of u* and excess stress

    Q_Einstein_Site{i} = (Q_Site{i}/1000)./(rho_s*g^(1/2).*(d50_Site{i}/1000).^(3/2));
    Q_norm2_Site{i} = Q_Site{i}*g./(1000*rho_a*ustThr_Site{i}); 
    Q_norm3_Site{i} = Q_Site{i}*g/(1000*rho_a);
    Q_nondim2_Site{i} = Q_norm2_Site{i}./(rho_a*tau_ex_Site{i});
    Q_nondim3_Site{i} = Q_norm3_Site{i}./(rho_a*tau_ex_Site{i}.*ust_Site{i});
    
    %add to lists
    ust_all = [ust_all; ust_Site{i}];
    ust_ratio_all = [ust_ratio_all; ust_ratio_Site{i}];
    tau_ratio_all = [tau_ratio_all; tau_ratio_Site{i}];
    tau_all = [tau_all; tau_Site{i}];
    ust_ex_all = [ust_ex_all; ust_ex_Site{i}];
    tau_ex_all = [tau_ex_all; tau_ex_Site{i}];
    tau_Shields_all = [tau_Shields_all; tau_Shields_Site{i}];
    ust_tau_ex_all = [ust_tau_ex_all; ust_tau_ex_Site{i}];
    Q_all = [Q_all; Q_Site{i}];
    Q_Einstein_all = [Q_Einstein_all; Q_Einstein_Site{i}];
    Q_norm2_all = [Q_norm2_all; Q_norm2_Site{i}];
    Q_norm3_all = [Q_norm3_all; Q_norm3_Site{i}];
    Q_nondim2_all = [Q_nondim2_all; Q_nondim2_Site{i}];
    Q_nondim3_all = [Q_nondim3_all; Q_nondim3_Site{i}];
    
    %PERFORM BINNING, COMPUTE STANDARD ERRORS
    for j=1:N_ust_bins
        bin_ind = find(ust_Site{i}>=ust_bins_min(j)&ust_Site{i}<=ust_bins_max(j));
        Q_ust_bin_avg{i}(j) = mean(Q_Site{i}(bin_ind));
        Q_ust_bin_SE{i}(j) = std(Q_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_ust_ratio_bins
        bin_ind = find(ust_ratio_Site{i}>=ust_ratio_bins_min(j)&ust_ratio_Site{i}<=ust_ratio_bins_max(j));
        Q_Einstein_ust_ratio_bin_avg{i}(j) = mean(Q_Einstein_Site{i}(bin_ind));
        Q_Einstein_ust_ratio_bin_SE{i}(j) = std(Q_Einstein_Site{i}(bin_ind))/sqrt(length(bin_ind));
        C_Q_nondim2_ust_ratio_bin_avg{i}(j) = mean(Q_nondim2_Site{i}(bin_ind));
        C_Q_nondim2_ust_ratio_bin_SE{i}(j) = std(Q_nondim2_Site{i}(bin_ind))/sqrt(length(bin_ind));
        C_Q_nondim3_ust_ratio_bin_avg{i}(j) = mean(Q_nondim3_Site{i}(bin_ind));
        C_Q_nondim3_ust_ratio_bin_SE{i}(j) = std(Q_nondim3_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_ust_ex_bins
        bin_ind = find(ust_ex_Site{i}>=ust_ex_bins_min(j)&ust_ex_Site{i}<=ust_ex_bins_max(j));
        Q_ust_ex_bin_avg{i}(j) = mean(Q_Site{i}(bin_ind));
        Q_ust_ex_bin_SE{i}(j) = std(Q_Site{i}(bin_ind))/sqrt(length(bin_ind));
        Q_Einstein_ust_ex_bin_avg{i}(j) = mean(Q_Einstein_Site{i}(bin_ind));
        Q_Einstein_ust_ex_bin_SE{i}(j) = std(Q_Einstein_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_tau_bins
        bin_ind = find(tau_Site{i}>=tau_bins_min(j)&tau_Site{i}<=tau_bins_max(j));
        Q_tau_bin_avg{i}(j) = mean(Q_Site{i}(bin_ind));
        Q_tau_bin_SE{i}(j) = std(Q_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_tau_ratio_bins
        bin_ind = find(tau_ratio_Site{i}>=tau_ratio_bins_min(j)&tau_ratio_Site{i}<=tau_ratio_bins_max(j));
        Q_Einstein_tau_ratio_bin_avg{i}(j) = mean(Q_Einstein_Site{i}(bin_ind));
        Q_Einstein_tau_ratio_bin_SE{i}(j) = std(Q_Einstein_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_tau_ex_bins
        bin_ind = find(tau_ex_Site{i}>=tau_ex_bins_min(j)&tau_ex_Site{i}<=tau_ex_bins_max(j));
        Q_tau_ex_bin_avg{i}(j) = mean(Q_Site{i}(bin_ind));
        Q_tau_ex_bin_SE{i}(j) = std(Q_Site{i}(bin_ind));
        Q_Einstein_tau_ex_bin_avg{i}(j) = mean(Q_Einstein_Site{i}(bin_ind));
        Q_Einstein_tau_ex_bin_SE{i}(j) = std(Q_Einstein_Site{i}(bin_ind));
        Q_norm2_tau_ex_bin_avg{i}(j) = mean(Q_norm2_Site{i}(bin_ind));
        Q_norm2_tau_ex_bin_SE{i}(j) = std(Q_norm2_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_tau_Shields_bins
        bin_ind = find(tau_Shields_Site{i}>=tau_Shields_bins_min(j)&tau_Shields_Site{i}<=tau_Shields_bins_max(j));
        Q_Einstein_tau_Shields_bin_avg{i}(j) = mean(Q_Einstein_Site{i}(bin_ind));
        Q_Einstein_tau_Shields_bin_SE{i}(j) = std(Q_Einstein_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
    
    for j=1:N_ust_tau_ex_bins
        bin_ind = find(ust_tau_ex_Site{i}>=ust_tau_ex_bins_min(j)&ust_tau_ex_Site{i}<=ust_tau_ex_bins_max(j));
        Q_norm3_ust_tau_ex_bin_avg{i}(j) = mean(Q_norm3_Site{i}(bin_ind));
        Q_norm3_ust_tau_ex_bin_SE{i}(j) = std(Q_norm3_Site{i}(bin_ind))/sqrt(length(bin_ind));
    end
end

%% GET STRESS AND SHEAR VELOCITY VALUES FOR FITTING

%get indices of ust>ustThr
ind_positive = find(ust_ratio_all>1);

%get tau values for plotting fits
tau_fit = linspace(min(tau_all), max(tau_all),50);
tau_ratio_fit = linspace(1, max(tau_ratio_all),50);
tau_ex_fit = linspace(0, max(tau_ex_all),50);
tau_Shields_fit = linspace(0, max(tau_Shields_all), 50);

%get u* values for plotting fits
ust_fit = sqrt(tau_fit/rho_a);
ust_ratio_fit = sqrt(tau_ratio_fit);
ust_ex_fit = linspace(0, max(ust_ex_all),50);
ust_tau_ex_fit = linspace(0, max(ust_tau_ex_all), 50);


%% PERFORM COMPARISONS

%fit Q versus u* and get correlation coefficient
P_Q_ust = polyfit(ust_all(ind_positive),Q_all(ind_positive),1);
Q_ust_fit = polyval(P_Q_ust,ust_fit);
R_Q_ust = corrcoef(ust_all(ind_positive),Q_all(ind_positive));
R_Q_ust = R_Q_ust(2);
chi2_Q_ust = sum((Q_all(ind_positive)-polyval(P_Q_ust,ust_all(ind_positive))).^2)/(var(Q_all(ind_positive))*length(ind_positive))

%fit Q versus tau and get correlation coefficient
P_Q_tau = polyfit(tau_all(ind_positive),Q_all(ind_positive),1);
Q_tau_fit = polyval(P_Q_tau,tau_fit);
R_Q_tau = corrcoef(tau_all(ind_positive),Q_all(ind_positive));
R_Q_tau = R_Q_tau(2);
chi2_Q_tau = sum((Q_all(ind_positive)-polyval(P_Q_tau,tau_all(ind_positive))).^2)/(var(Q_all(ind_positive))*length(ind_positive))

%fit Q versus excess ust and get correlation coefficient
P_Q_ust_ex = polyfit(ust_ex_all(ind_positive),Q_all(ind_positive),1);
Q_ust_ex_fit = polyval(P_Q_ust_ex, ust_ex_fit);
R_Q_ust_ex = corrcoef(ust_ex_all(ind_positive),Q_all(ind_positive));
R_Q_ust_ex = R_Q_ust_ex(2);
chi2_Q_ust_ex = sum((Q_all(ind_positive)-polyval(P_Q_ust_ex,ust_ex_all(ind_positive))).^2)/(var(Q_all(ind_positive))*length(ind_positive))

%fit Q versus excess tau and get correlation coefficient
P_Q_tau_ex = polyfit(tau_ex_all(ind_positive),Q_all(ind_positive),1);
Q_tau_ex_fit = polyval(P_Q_tau_ex, tau_ex_fit);
R_Q_tau_ex = corrcoef((tau_ex_all(ind_positive)),Q_all(ind_positive));
R_Q_tau_ex = R_Q_tau_ex(2);
chi2_Q_tau_ex = sum((Q_all(ind_positive)-polyval(P_Q_tau_ex,tau_ex_all(ind_positive))).^2)/(var(Q_all(ind_positive))*length(ind_positive))

%fit Einstein Q versus ust excess, get correlation coefficient
P_Q_Einstein_ust_ex = polyfit(ust_ex_all(ind_positive),Q_Einstein_all(ind_positive),1);
Q_Einstein_ust_ex_fit = polyval(P_Q_Einstein_ust_ex, ust_ex_fit);
R_Q_Einstein_ust_ex = corrcoef(ust_ex_all(ind_positive),Q_Einstein_all(ind_positive));
R_Q_Einstein_ust_ex = R_Q_Einstein_ust_ex(2);
chi2_Q_Einstein_ust_ex = sum((Q_Einstein_all(ind_positive)-polyval(P_Q_Einstein_ust_ex,ust_ex_all(ind_positive))).^2)/(var(Q_Einstein_all(ind_positive))*length(ind_positive))

%fit Einstein Q versus tau excess, get correlation coefficient
P_Q_Einstein_tau_ex = polyfit(tau_ex_all(ind_positive),Q_Einstein_all(ind_positive),1);
Q_Einstein_tau_ex_fit = polyval(P_Q_Einstein_tau_ex, tau_ex_fit);
R_Q_Einstein_tau_ex = corrcoef(tau_ex_all(ind_positive),Q_Einstein_all(ind_positive));
R_Q_Einstein_tau_ex = R_Q_Einstein_tau_ex(2);
chi2_Q_Einstein_tau_ex = sum((Q_Einstein_all(ind_positive)-polyval(P_Q_Einstein_tau_ex,tau_ex_all(ind_positive))).^2)/(var(Q_Einstein_all(ind_positive))*length(ind_positive))

%fit Einstein Q versus ust ratio and get correlation coefficient
P_Q_Einstein_ust_ratio = polyfit(ust_ratio_all(ind_positive),Q_Einstein_all(ind_positive),1);
Q_Einstein_ust_ratio_fit = polyval(P_Q_Einstein_ust_ratio,ust_ratio_fit);
R_Q_Einstein_ust_ratio = corrcoef(ust_ratio_all(ind_positive),Q_Einstein_all(ind_positive));
R_Q_Einstein_ust_ratio = R_Q_Einstein_ust_ratio(2);
chi2_Q_Einstein_ust_ratio = sum((Q_Einstein_all(ind_positive)-polyval(P_Q_Einstein_ust_ratio,ust_ratio_all(ind_positive))).^2)/(var(Q_Einstein_all(ind_positive))*length(ind_positive))

%fit Einstein Q versus stress ratio and get correlation coefficient
P_Q_Einstein_tau_ratio = polyfit(tau_ratio_all(ind_positive),Q_Einstein_all(ind_positive),1);
Q_Einstein_tau_ratio_fit = polyval(P_Q_Einstein_tau_ratio,tau_ratio_fit);
R_Q_Einstein_tau_ratio = corrcoef(tau_ratio_all(ind_positive),Q_Einstein_all(ind_positive));
R_Q_Einstein_tau_ratio = R_Q_Einstein_tau_ratio(2);
chi2_Q_Einstein_tau_ratio = sum((Q_Einstein_all(ind_positive)-polyval(P_Q_Einstein_tau_ratio,tau_ratio_all(ind_positive))).^2)/(var(Q_Einstein_all(ind_positive))*length(ind_positive))

%fit Einstein Q versus Shields stress and get correlation coefficient
P_Q_Einstein_tau_Shields = polyfit(tau_Shields_all(ind_positive),Q_Einstein_all(ind_positive),1);
Q_Einstein_tau_Shields_fit = polyval(P_Q_Einstein_tau_Shields,tau_Shields_fit);
R_Q_Einstein_tau_Shields = corrcoef(tau_Shields_all(ind_positive),Q_Einstein_all(ind_positive));
R_Q_Einstein_tau_Shields = R_Q_Einstein_tau_Shields(2);
chi2_Q_Einstein_tau_Shields = sum((Q_Einstein_all(ind_positive)-polyval(P_Q_Einstein_tau_Shields,tau_Shields_all(ind_positive))).^2)/(var(Q_Einstein_all(ind_positive))*length(ind_positive))

%fit Q_norm2 versus tau excess, get correlation coefficient
P_Q_norm2_tau_ex = polyfit(tau_ex_all(ind_positive),Q_norm2_all(ind_positive),1);
Q_norm2_tau_ex_fit = polyval(P_Q_norm2_tau_ex, tau_ex_fit);
R_Q_norm2_tau_ex = corrcoef(tau_ex_all(ind_positive),Q_norm2_all(ind_positive));
R_Q_norm2_tau_ex = R_Q_norm2_tau_ex(2);
chi2_Q_norm2_tau_ex = sum((Q_norm2_all(ind_positive)-polyval(P_Q_norm2_tau_ex,tau_ex_all(ind_positive))).^2)/(var(Q_norm2_all(ind_positive))*length(ind_positive))

%fit Q_norm3 versus ust times tau excess, get correlation coefficient
P_Q_norm3_ust_tau_ex = polyfit(ust_tau_ex_all(ind_positive),Q_norm3_all(ind_positive),1);
Q_norm3_ust_tau_ex_fit = polyval(P_Q_norm3_ust_tau_ex, ust_tau_ex_fit);
R_Q_norm3_ust_tau_ex = corrcoef(ust_tau_ex_all(ind_positive),Q_norm3_all(ind_positive));
R_Q_norm3_ust_tau_ex = R_Q_norm3_ust_tau_ex(2);
chi2_Q_norm3_ust_tau_ex = sum((Q_norm3_all(ind_positive)-polyval(P_Q_norm3_ust_tau_ex,ust_tau_ex_all(ind_positive))).^2)/(var(Q_norm3_all(ind_positive))*length(ind_positive))

%get median values of Q_nondim2 and Q_nondim3
C_Q_nondim2_ust_ratio_median = median(Q_nondim2_all(ind_positive));
C_Q_nondim3_ust_ratio_median = median(Q_nondim3_all(ind_positive));


%% MAKE CALCULATIONS FOR LITERATURE VALUES

%thresholds
ustThr_Greeley96 = C_ustThr_sqrtd50*sqrt(d50_Greeley96); %get Greeley ustThr value
ustThr_Namikas03 = C_ustThr_sqrtd50*sqrt(d50_Namikas03); %get Namikas ustThr value

%shear velocities
ust_ratio_Greeley96 = ust_Greeley96/ustThr_Greeley96; %get ust/ustThr for Greeley
ust_ratio_Namikas03 = ust_Namikas03/ustThr_Namikas03; %get ust/ustThr for Namikas
ust_ex_Greeley96 = ust_Greeley96-ustThr_Greeley96; %get ust_ex for Greeley data
ust_ex_Namikas03 = ust_Namikas03-ustThr_Namikas03; %get ust_ex for Namikas data

%stresses
tauThr_Greeley96 = rho_a*ustThr_Greeley96.^2; %get Greeley tauThr value
tauThr_Namikas03 = rho_a*ustThr_Namikas03.^2; %get Namikas tauThr value
tau_Greeley96 = rho_a*ust_Greeley96.^2; %get tau for Greeley data
tau_Namikas03 = rho_a*ust_Namikas03.^2; %get tau for Namikas data
tau_ratio_Greeley96 = ust_ratio_Greeley96.^2; %get tau/tauThr for Greeley
tau_ratio_Namikas03 = ust_ratio_Namikas03.^2; %get tau/tauThr for Namikas
tau_ex_Greeley96 = tau_Greeley96-tauThr_Greeley96; %get tau_ex for Greeley data
tau_ex_Namikas03 = tau_Namikas03-tauThr_Namikas03; %get tau_ex for Namikas data
tau_Shields_Greeley96 = tau_Greeley96./(rho_s*g*(d50_Greeley96/1000)); %get Shields stress for Greeley data
tau_Shields_Namikas03 = tau_Namikas03./(rho_s*g*(d50_Namikas03/1000)); %get Shields stress for Namikas data
ust_tau_ex_Greeley96 = ust_Greeley96.*tau_ex_Greeley96; %u* tau_ex product
ust_tau_ex_Namikas03 = ust_Namikas03.*tau_ex_Namikas03; %u* tau_ex product

%fluxes
Q_Einstein_Greeley96 = (Q_Greeley96/1000)./(rho_s*g^(1/2).*(d50_Greeley96/1000).^(3/2)); %get Einstein flux for Greeley
Q_norm2_Greeley96 = Q_Greeley96*g./(1000*rho_a*ustThr_Greeley96); 
Q_norm3_Greeley96 = Q_Greeley96*g/(1000*rho_a);
Q_nondim2_Greeley96 = Q_norm2_Greeley96./(rho_a*tau_ex_Greeley96);
Q_nondim3_Greeley96 = Q_norm3_Greeley96./(rho_a*tau_ex_Greeley96.*ust_Greeley96);
Q_Einstein_Namikas03 = (Q_Namikas03/1000)./(rho_s*g^(1/2).*(d50_Greeley96/1000).^(3/2)); %get Einstein flux for Namikas
Q_norm2_Namikas03 = Q_Namikas03*g./(1000*rho_a*ustThr_Namikas03); 
Q_norm3_Namikas03 = Q_Namikas03*g/(1000*rho_a);
Q_nondim2_Namikas03 = Q_norm2_Namikas03./(rho_a*tau_ex_Namikas03);
Q_nondim3_Namikas03 = Q_norm3_Namikas03./(rho_a*tau_ex_Namikas03.*ust_Namikas03);

%% GENERATE PLOTS
close all; %close all windows if they were open

%FLUX VERSUS SHEAR VELOCITY - NO NORMALIZATION
figure(1); clf; hold on; %initialize plot
for i = 1:N_Sites
    %errorbar(ust_bins_mid,Q_ust_bin_avg{i},Q_ust_bin_SE{i},Markers{i},'MarkerSize',5);
    plot(ust_Site{i},Q_Site{i},Markers{i},'MarkerSize',2);
end
plot(ust_fit, Q_ust_fit,'k-','LineWidth',2); %fit
plot(ust_Greeley96,Q_Greeley96,'k^'); %add Greeley data to plot
plot(ust_Namikas03,Q_Namikas03,'kd'); %add Namikas data to plot
ylim([0 ceil(max(cell2mat(Q_Site)))]);
xlabel('u_{*} (m/s)');
ylabel('Q (g m^{-1} s^{-1})');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['linear fit, R=',num2str(R_Q_ust)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
title('Flux versus shear velocity, no normalization');
print([folder_Plots,'Flux_Ust_all.png'],'-dpng');

%FLUX VERSUS STRESS - NO NORMALIZATION
figure(2); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(tau_bins_mid,Q_tau_bin_avg{i},Q_tau_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(tau_Site{i},Q_Site{i},Markers{i},'MarkerSize',2);
end
plot(tau_fit, Q_tau_fit,'k-','LineWidth',2); %fit
plot(tau_Greeley96,Q_Greeley96,'k^'); %add Greeley data to plot
plot(tau_Namikas03,Q_Namikas03,'kd'); %add Namikas data to plot
ylim([0 ceil(max(cell2mat(Q_Site)))]);
xlabel('\tau (Pa)');
ylabel('Q (g m^{-1} s^{-1})');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['linear fit, R=',num2str(R_Q_tau)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
title('Flux versus stress, no normalization');
print([folder_Plots,'Flux_Tau_all.png'],'-dpng');

%FLUX VERSUS EXCESS SHEAR VELOCITY
figure(3); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(ust_ex_bins_mid,Q_ust_ex_bin_avg{i},Q_ust_ex_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(ust_ex_Site{i},Q_Site{i},Markers{i},'MarkerSize',2);
end
plot(ust_ex_fit,Q_ust_ex_fit,'k','LineWidth',2); %add fit to plot
plot(ust_ex_Greeley96,Q_Greeley96,'k^'); %add Greeley data to plot
plot(ust_ex_Namikas03,Q_Namikas03,'kd'); %add Namikas data to plot
xlabel('u_{*} - u_{*,th} (m/s)');
ylabel('Q (g m^{1} s^{-1})');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_ust_ex)];
legend_values{length(Sites)+2} = ['Greeley et al. (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
title('Flux versus excess shear velocity');
print([folder_Plots,'Flux_UstEx_all.png'],'-dpng');

%FLUX VERSUS EXCESS STRESS
figure(4); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(tau_ex_bins_mid,Q_tau_ex_bin_avg{i},Q_tau_ex_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(tau_ex_Site{i},Q_Site{i},Markers{i},'MarkerSize',2);
end
plot(tau_ex_fit,Q_tau_ex_fit,'k','LineWidth',2); %add fit to plot
plot(tau_ex_Greeley96,Q_Greeley96,'k^'); %add Greeley data to plot
plot(tau_ex_Namikas03,Q_Namikas03,'kd'); %add Namikas data to plot
xlabel('\tau - \tau_{th} (Pa)');
ylabel('Q (g m^{-1} s^{-1})');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_tau_ex)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
title('Flux versus excess stress');
print([folder_Plots,'Flux_TauEx_all.png'],'-dpng');

%EINSTEIN FLUX VERSUS EXCESS SHEAR VELOCITY
figure(5); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(ust_ex_bins_mid,Q_Einstein_ust_ex_bin_avg{i},Q_Einstein_ust_ex_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(ust_ex_Site{i},Q_Einstein_Site{i},Markers{i},'MarkerSize',2);
end
plot(ust_ex_fit,Q_Einstein_ust_ex_fit,'k','LineWidth',2); %add fit to plot
plot(ust_ex_Greeley96,Q_Einstein_Greeley96,'k^'); %add Greeley data to plot
plot(ust_ex_Namikas03,Q_Einstein_Namikas03,'kd'); %add Namikas data to plot
xlabel('u_{*} - u_{*,th} (m/s)','FontSize',16);
ylabel('Q_{*} = Q \rho_{s}^{-1} g^{-1/2} D^{-3/2}','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_Einstein_ust_ex)];
legend_values{length(Sites)+2} = ['Greeley et al. (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
title('Einstein flux versus excess shear velocity');
print([folder_Plots,'FluxEinstein_UstEx_all.png'],'-dpng');

%EINSTEIN FLUX VERSUS EXCESS STRESS
figure(6); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(tau_ex_bins_mid,Q_Einstein_tau_ex_bin_avg{i},Q_Einstein_tau_ex_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(tau_ex_Site{i},Q_Einstein_Site{i},Markers{i},'MarkerSize',2);
end
plot(tau_ex_fit,Q_Einstein_tau_ex_fit,'k','LineWidth',2); %add fit to plot
plot(tau_ex_Greeley96,Q_Einstein_Greeley96,'k^'); %add Greeley data to plot
plot(tau_ex_Namikas03,Q_Einstein_Namikas03,'kd'); %add Namikas data to plot
xlabel('\tau - \tau_{th} (Pa)');
ylabel('Q_{*} = Q \rho_{s}^{-1} g^{-1/2} D^{-3/2}');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_Einstein_tau_ex)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
title('Einstein flux versus excess stress');
print([folder_Plots,'FluxEinstein_TauEx_all.png'],'-dpng');

%EINSTEIN FLUX VERSUS SHEAR VELOCITY RATIO
figure(7); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(ust_ratio_bins_mid,Q_Einstein_ust_ratio_bin_avg{i},Q_Einstein_ust_ratio_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(ust_ratio_Site{i},Q_Einstein_Site{i},Markers{i},'MarkerSize',2);
end
plot(ust_ratio_fit,Q_Einstein_ust_ratio_fit,'k','LineWidth',2); %add fit to plot
plot(ust_ratio_Greeley96,Q_Einstein_Greeley96,'k^'); %add Greeley data to plot
plot(ust_ratio_Namikas03,Q_Einstein_Namikas03,'kd'); %add Namikas data to plot
xlabel('u_{*}/u_{*,th}');
ylabel('Q_{*} = Q \rho_{s}^{-1} g^{-1/2} D^{-3/2}');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_Einstein_ust_ratio)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
title('Einstein flux versus shear velocity ratio');
print([folder_Plots,'FluxEinstein_UstRatio_all.png'],'-dpng');

%EINSTEIN FLUX VERSUS STRESS RATIO
figure(8); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(tau_ratio_bins_mid,Q_Einstein_tau_ratio_bin_avg{i},Q_Einstein_tau_ratio_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(tau_ratio_Site{i},Q_Einstein_Site{i},Markers{i},'MarkerSize',2);
end
plot(tau_ratio_fit,Q_Einstein_tau_ratio_fit,'k','LineWidth',2); %add fit to plot
plot(tau_ratio_Greeley96,Q_Einstein_Greeley96,'k^'); %add Greeley data to plot
plot(tau_ratio_Namikas03,Q_Einstein_Namikas03,'kd'); %add Namikas data to plot
xlabel('\tau/\tau_{th}');
ylabel('Q_{*} = Q \rho_{s}^{-1} g^{-1/2} D^{-3/2}');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_Einstein_tau_ratio)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
title('Einstein flux versus stress ratio');
print([folder_Plots,'FluxEinstein_TauRatio_all.png'],'-dpng');

%EINSTEIN FLUX VERSUS SHIELDS STRESS
figure(9); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(tau_Shields_bins_mid,Q_Einstein_tau_Shields_bin_avg{i},Q_Einstein_tau_Shields_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(tau_Shields_Site{i},Q_Einstein_Site{i},Markers{i},'MarkerSize',2);
end
plot(tau_Shields_fit,Q_Einstein_tau_Shields_fit,'k','LineWidth',2); %add fit to plot
plot(tau_Shields_Greeley96,Q_Einstein_Greeley96,'k^'); %add Greeley data to plot
plot(tau_Shields_Namikas03,Q_Einstein_Namikas03,'kd'); %add Namikas data to plot
xlabel('\tau_{*} = \tau \rho_{s}^{-1} g^{-1} D^{-1}');
ylabel('Q_{*} = Q \rho_{s}^{-1} g^{-1/2} D^{-3/2}');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_Einstein_tau_Shields)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
title('Einstein flux versus Shields stress');
print([folder_Plots,'FluxEinstein_TauShields_all.png'],'-dpng');

%SQUARE NORMALIZED FLUX VERSUS EXCESS SHEAR STRESS
figure(10); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(tau_ex_bins_mid,Q_norm2_tau_ex_bin_avg{i},Q_norm2_tau_ex_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(tau_ex_Site{i},Q_norm2_Site{i},Markers{i},'MarkerSize',2);
end
plot(tau_ex_fit,Q_norm2_tau_ex_fit,'k','LineWidth',2); %add fit to plot
plot(tau_ex_Greeley96,Q_norm2_Greeley96,'k^'); %add Greeley data to plot
plot(tau_ex_Namikas03,Q_norm2_Namikas03,'kd'); %add Namikas data to plot
xlabel('\tau_{ex} (Pa)','FontSize',16);
ylabel('gQ/{(\rho_{a}u_{*,th})} (m^2/s^2)','FontSize',16);
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_norm2_tau_ex)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','NorthWest');
set(h_legend,'FontSize',16);
title('Squared normalization');
print([folder_Plots,'FluxNorm2_TauEx_all.png'],'-dpng');

%CUBIC NORMALIZED FLUX VERSUS SHEAR VELOCITY CUBED
figure(11); clf; hold on; %initialize plot
for i = 1:N_Sites
    errorbar(ust_tau_ex_bins_mid,Q_norm3_ust_tau_ex_bin_avg{i},Q_norm3_ust_tau_ex_bin_SE{i},Markers{i},'MarkerSize',5);
%    plot(ust_tau_ex_Site{i}/rho_a,Q_norm3_Site{i},Markers{i},'MarkerSize',2);
end
plot(ust_tau_ex_fit,Q_norm3_ust_tau_ex_fit,'k','LineWidth',2); %add fit to plot
plot(ust_tau_ex_Greeley96,Q_norm3_Greeley96,'k^'); %add Greeley data to plot
plot(ust_tau_ex_Namikas03,Q_norm3_Namikas03,'kd'); %add Namikas data to plot
xlabel('\tau_{ex}u_{*} (Pa-m/s)');
ylabel('gQ/{\rho_{a}} (m^3/s^3)');
set(gca,'FontSize',16);
legend_values = Sites;
legend_values{length(Sites)+1} = ['fit, R= ',num2str(R_Q_norm3_ust_tau_ex)];
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthEast');
set(h_legend,'FontSize',16);
title('Cubic normalization');
print([folder_Plots,'FluxNorm3_Ust3_all.png'],'-dpng');

%SQUARED DIMENSIONLESS FLUX VERSUS SHEAR VELOCITY RATIO
figure(12); clf; hold on; %initialize plot
for i = 1:N_Sites
    notoutlierind = find(C_Q_nondim2_ust_ratio_bin_avg{i}<10);
    errorbar(ust_ratio_bins_mid(notoutlierind),C_Q_nondim2_ust_ratio_bin_avg{i}(notoutlierind),C_Q_nondim2_ust_ratio_bin_SE{i}(notoutlierind),Markers{i},'MarkerSize',5);
%    plot(ust_ratio_Site{i},Q_nondim2_Site{i},Markers{i},'MarkerSize',2);
end
plot(xlim,C_Q_nondim2_ust_ratio_median*[1 1],'k','LineWidth',2); %add fit to plot
plot(ust_ratio_Greeley96,Q_nondim2_Greeley96,'k^'); %add Greeley data to plot
plot(ust_ratio_Namikas03,Q_nondim2_Namikas03,'kd'); %add Namikas data to plot
xlabel('u_{*}/u_{*,th}','FontSize',16);
ylabel('gQ/{(\tau_{ex}u_{*,th})}','FontSize',16);
xlim([1 max(ust_ratio_all)]);
title('Squared law');
set(gca,'FontSize',16);
set(gca,'xscale','log');
legend_values = Sites;
legend_values{length(Sites)+1} = 'median';
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthWest');
set(h_legend,'FontSize',16);
title('Squared dimensionless flux versus shear velocity ratio');
print([folder_Plots,'FluxNondim2_UstRatio_all.png'],'-dpng');

%CUBIC DIMENSIONLESS FLUX VERSUS SHEAR VELOCITY RATIO
figure(13); clf; hold on; %initialize plot
for i = 1:N_Sites
    notoutlierind = find(C_Q_nondim3_ust_ratio_bin_avg{i}<10);
    errorbar(ust_ratio_bins_mid(notoutlierind),C_Q_nondim3_ust_ratio_bin_avg{i}(notoutlierind),C_Q_nondim3_ust_ratio_bin_SE{i}(notoutlierind),Markers{i},'MarkerSize',5);
%    plot(ust_ratio_Site{i},Q_nondim3_Site{i},Markers{i},'MarkerSize',2);
end
plot(xlim,C_Q_nondim3_ust_ratio_median*[1 1],'k','LineWidth',2); %add fit to plot
plot(ust_ratio_Greeley96,Q_nondim3_Greeley96,'k^'); %add Greeley data to plot
plot(ust_ratio_Namikas03,Q_nondim3_Namikas03,'kd'); %add Namikas data to plot
xlabel('u_{*}/u_{*,th}');
ylabel('gQ/{(\tau_{ex}u_{*})}');
xlim([1 max(ust_ratio_all)]);
title('Cubed law');
set(gca,'FontSize',16);
set(gca,'xscale','log');
legend_values = Sites;
legend_values{length(Sites)+1} = 'median';
legend_values{length(Sites)+2} = ['Greeley (1996)'];
legend_values{length(Sites)+3} = ['Namikas (2003)'];
h_legend = legend(legend_values,'Location','SouthWest');
set(h_legend,'FontSize',16);
title('Cubic dimensionless flux versus shear velocity ratio');
print([folder_Plots,'FluxNondim3_UstRatio_all.png'],'-dpng');