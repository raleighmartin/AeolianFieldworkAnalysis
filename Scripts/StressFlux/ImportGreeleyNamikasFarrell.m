%% SCRIPT TO IMPORT GREELEY (1996), NAMIKAS (2003), AND FARRELL (2012) FLUX VALUES FROM TEXT FILE
clearvars; %initialize


%% constants
ust_relerr = 0.1; %10% relative error for u*


%% filepath information
folder_LiteratureData = '../../../../Google Drive/Data/Literature/'; %location where Namikas and Greeley data are stored
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for storing outputs of this analysis
folder_Plots = '../../PlotOutput/StressFlux/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions
SaveData_Path = strcat(folder_AnalysisData,'LitData'); %folder for saving data

%% processing variables
delimiter = '\t';
startRow = 2;


%% manually enter basic info

% Greeley et al. 1996
filename_Greeley96 = [folder_LiteratureData,'Greeley96.dat'];
formatSpec_Greeley96 = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
N_Greeley96 = 10;
ust_Greeley96 = [0.47; 0.51; 0.49; 0.31; 0.35; 0.31; 0.54; 0.49; 0.41; 0.42]; %u* values (m/s) used in Kok and Renno (2008)
sigma_ust_Greeley96 = ust_Greeley96*ust_relerr; %calculate u* uncertainty based on assumed relative error
Q_pub_Greeley96 = [NaN; NaN; NaN; NaN; 16.1; 2.7; 19; 19; 10.4; 10.4]; %values from Table 1 in paper
d50_Greeley96 = 0.230; %median grain size (mm) value from paper
%Greeley96_OutlierInd = [4];
Greeley96_OutlierInd = [];

% Namikas 2003
filename_Namikas03 = [folder_LiteratureData,'Namikas03verMF.dat'];
formatSpec_Namikas03 = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
N_Namikas03 = 9;
ust_Namikas03 = [0.32; 0.27; 0.32; 0.30; 0.38; 0.37; 0.38; 0.47; 0.63]; %u* values (m/s) from paper
sigma_ust_Namikas03 = ust_Namikas03*ust_relerr; %calculate u* uncertainty based on assumed relative error
Q_pub_Namikas03 = [0.35; 0.76; 2.73; 2.12; 5.61; 5.86; 5.57; 11.80; 29.11]; %values from Table 1 in paper
d50_Namikas03 = 0.250; %median grain size (mm) value from paper
%Namikas03_OutlierInd = [9];
Namikas03_OutlierInd = [];

% Farrell 2012
filename_Farrell12 = [folder_LiteratureData,'FarrellEtAl2012.dat'];
formatSpec_Farrell12 = '%f%f%f%f%f%f%f%[^\n\r]';
N_Farrell12 = 14;
ust_Farrell12 = [0.54; 0.47; 0.53; 0.49; 0.50; 0.50; 0.47; 0.45; 0.51; 0.48; 0.41; 0.50; 0.50; 0.49]; %u* values (m/s) from paper
sigma_ust_Farrell12 = ust_Farrell12*ust_relerr; %calculate u* uncertainty based on assumed relative error
d50_Farrell12 = 0.545659491010631; %use median grain size (mm) from Jeri deployment since they don't report it
Farrell12_OutlierInd = [1,12:14];
%exclude first sample because of shorter time interval
%exclude last three because they are from coarser CSFC (as opposed to the rest of the samples from finer CSFF)


%% Load and process Greeley et al. 1996 data
Greeley96_GoodInd = setdiff(1:N_Greeley96,Greeley96_OutlierInd); %information about which values to include in analysis
fileID_Greeley96 = fopen(filename_Greeley96,'r');
dataArray_Greeley96 = textscan(fileID_Greeley96, formatSpec_Greeley96, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID_Greeley96);

% Initialize cell arrays with all values
z_Greeley96 = cell(N_Greeley96,1); %measurement heights (m)
q_Greeley96 = cell(N_Greeley96,1); %fluxes (g/m^2/s)
RunName_Greeley96 = cell(N_Greeley96,1);

% Allocate imported array to column variable names and cell arrays
z_cm_1a = dataArray_Greeley96{:, 1};
q_gpercm2s_1a = dataArray_Greeley96{:, 2};
z_Greeley96{1} = 1e-2*z_cm_1a(~isnan(z_cm_1a));
q_Greeley96{1} = 1e4*q_gpercm2s_1a(~isnan(q_gpercm2s_1a));
RunName_Greeley96{1} = '1a';

z_cm_1b = dataArray_Greeley96{:, 3};
q_gpercm2s_1b = dataArray_Greeley96{:, 4};
z_Greeley96{2} = 1e-2*z_cm_1b(~isnan(z_cm_1b));
q_Greeley96{2} = 1e4*q_gpercm2s_1b(~isnan(q_gpercm2s_1b));
RunName_Greeley96{2} = '1b';

z_cm_1c = dataArray_Greeley96{:, 5};
q_gpercm2s_1c = dataArray_Greeley96{:, 6};
z_Greeley96{3} = 1e-2*z_cm_1c(~isnan(z_cm_1c));
q_Greeley96{3} = 1e4*q_gpercm2s_1c(~isnan(q_gpercm2s_1c));
RunName_Greeley96{3} = '1c';

z_cm_2 = dataArray_Greeley96{:, 7};
q_gpercm2s_2 = dataArray_Greeley96{:, 8};
z_Greeley96{4} = 1e-2*z_cm_2(~isnan(z_cm_2));
q_Greeley96{4} = 1e4*q_gpercm2s_2(~isnan(q_gpercm2s_2));
RunName_Greeley96{4} = '2';

z_cm_3 = dataArray_Greeley96{:, 9};
q_gpercm2s_3 = dataArray_Greeley96{:, 10};
z_Greeley96{5} = 1e-2*z_cm_3(~isnan(z_cm_3));
q_Greeley96{5} = 1e4*q_gpercm2s_3(~isnan(q_gpercm2s_3));
RunName_Greeley96{5} = '3';

z_cm_4 = dataArray_Greeley96{:, 11};
q_gpercm2s_4 = dataArray_Greeley96{:, 12};
z_Greeley96{6} = 1e-2*z_cm_4(~isnan(z_cm_4));
q_Greeley96{6} = 1e4*q_gpercm2s_4(~isnan(q_gpercm2s_4));
RunName_Greeley96{6} = '4';

z_cm_5a = dataArray_Greeley96{:, 13};
q_gpercm2s_5a = dataArray_Greeley96{:, 14};
z_Greeley96{7} = 1e-2*z_cm_5a(~isnan(z_cm_5a));
q_Greeley96{7} = 1e4*q_gpercm2s_5a(~isnan(q_gpercm2s_5a));
RunName_Greeley96{7} = '5a';

z_cm_5b = dataArray_Greeley96{:, 15};
q_gpercm2s_5b = dataArray_Greeley96{:, 16};
z_Greeley96{8} = 1e-2*z_cm_5b(~isnan(z_cm_5b));
q_Greeley96{8} = 1e4*q_gpercm2s_5b(~isnan(q_gpercm2s_5b));
RunName_Greeley96{8} = '5b';

z_cm_6a = dataArray_Greeley96{:, 17};
q_gpercm2s_6a = dataArray_Greeley96{:, 18};
z_Greeley96{9} = 1e-2*z_cm_6a(~isnan(z_cm_6a));
q_Greeley96{9} = 1e4*q_gpercm2s_6a(~isnan(q_gpercm2s_6a));
RunName_Greeley96{9} = '6a';

z_cm_6b = dataArray_Greeley96{:, 19};
q_gpercm2s_6b = dataArray_Greeley96{:, 20};
z_Greeley96{10} = 1e-2*z_cm_6b(~isnan(z_cm_6b));
q_Greeley96{10} = 1e4*q_gpercm2s_6b(~isnan(q_gpercm2s_6b));
RunName_Greeley96{10} = '6b';

% Calculate saltation flux and height
zq_Greeley96 = zeros(N_Greeley96,1); %(m)
sigma_zq_Greeley96 = zeros(N_Greeley96,1); %(m)
Q_fit_Greeley96 = zeros(N_Greeley96,1); %(g/m/s)

for i=1:N_Greeley96
    [q0,zq,~,sigma_q0,sigma_zq,~] = qz_profilefit(q_Greeley96{i}, z_Greeley96{i}); %fit will assume equal uncertainties
    Q = q0*zq;
    zq_Greeley96(i) = zq;
    sigma_zq_Greeley96(i) = sigma_zq; %flux height uncertainty
    Q_fit_Greeley96(i) = Q;
end

%% Remove outlier value(s) from data
N_Greeley96 = length(Greeley96_GoodInd);
q_Greeley96 = q_Greeley96(Greeley96_GoodInd);
Q_pub_Greeley96 = Q_pub_Greeley96(Greeley96_GoodInd);
Q_fit_Greeley96 =Q_fit_Greeley96(Greeley96_GoodInd);
RunName_Greeley96 = RunName_Greeley96(Greeley96_GoodInd);
ust_Greeley96 = ust_Greeley96(Greeley96_GoodInd);
z_Greeley96 = z_Greeley96(Greeley96_GoodInd);
zq_Greeley96 = zq_Greeley96(Greeley96_GoodInd);
sigma_zq_Greeley96 = sigma_zq_Greeley96(Greeley96_GoodInd);


%% Load and process Namikas 2003 data
Namikas03_GoodInd = setdiff(1:N_Namikas03,Namikas03_OutlierInd); %information about which values to include in analysis
fileID_Namikas03 = fopen(filename_Namikas03,'r');
dataArray_Namikas03 = textscan(fileID_Namikas03, formatSpec_Namikas03, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID_Namikas03);

% Initialize cell arrays with all values
z_Namikas03 = cell(N_Namikas03,1); %trap height (m)
q_Namikas03 = cell(N_Namikas03,1); %flux (g/m^2/s)
RunName_Namikas03 = cell(N_Namikas03,1);

% Allocate imported array to column variable names and cell arrays
z_cm_5 = dataArray_Namikas03{:, 1};
q_kgperm2s_5 = dataArray_Namikas03{:, 2};
z_Namikas03{1} = z_cm_5(~isnan(q_kgperm2s_5));
q_Namikas03{1} = 1e3*q_kgperm2s_5(~isnan(q_kgperm2s_5));
RunName_Namikas03{1} = '5';

z_cm_3 = dataArray_Namikas03{:, 1};
q_kgperm2s_3 = dataArray_Namikas03{:, 3};
z_Namikas03{2} = z_cm_3(~isnan(q_kgperm2s_3));
q_Namikas03{2} = 1e3*q_kgperm2s_3(~isnan(q_kgperm2s_3));
RunName_Namikas03{2} = '3';

z_cm_4 = dataArray_Namikas03{:, 1};
q_kgperm2s_4 = dataArray_Namikas03{:, 4};
z_Namikas03{3} = z_cm_4(~isnan(q_kgperm2s_4));
q_Namikas03{3} = 1e3*q_kgperm2s_4(~isnan(q_kgperm2s_4));
RunName_Namikas03{3} = '4';

z_cm_8 = dataArray_Namikas03{:, 1};
q_kgperm2s_8 = dataArray_Namikas03{:, 5};
z_Namikas03{4} = z_cm_8(~isnan(q_kgperm2s_8));
q_Namikas03{4} = 1e3*q_kgperm2s_8(~isnan(q_kgperm2s_8));
RunName_Namikas03{4} = '8';

z_cm_10 = dataArray_Namikas03{:, 1};
q_kgperm2s_10 = dataArray_Namikas03{:, 6};
z_Namikas03{5} = z_cm_10(~isnan(q_kgperm2s_10));
q_Namikas03{5} = 1e3*q_kgperm2s_10(~isnan(q_kgperm2s_10));
RunName_Namikas03{5} = '10';

z_cm_6 = dataArray_Namikas03{:, 1};
q_kgperm2s_6 = dataArray_Namikas03{:, 7};
z_Namikas03{6} = z_cm_6(~isnan(q_kgperm2s_6));
q_Namikas03{6} = 1e3*q_kgperm2s_6(~isnan(q_kgperm2s_6));
RunName_Namikas03{6} = '6';

z_cm_9 = dataArray_Namikas03{:, 1};
q_kgperm2s_9 = dataArray_Namikas03{:, 8};
z_Namikas03{7} = z_cm_9(~isnan(q_kgperm2s_9));
q_Namikas03{7} = 1e3*q_kgperm2s_9(~isnan(q_kgperm2s_9));
RunName_Namikas03{7} = '9';

z_cm_13 = dataArray_Namikas03{:, 1};
q_kgperm2s_13 = dataArray_Namikas03{:, 9};
z_Namikas03{8} = z_cm_13(~isnan(q_kgperm2s_13));
q_Namikas03{8} = 1e3*q_kgperm2s_13(~isnan(q_kgperm2s_13));
RunName_Namikas03{8} = '13';

z_cm_14 = dataArray_Namikas03{:, 1};
q_kgperm2s_14 = dataArray_Namikas03{:, 10};
z_Namikas03{9} = z_cm_14(~isnan(q_kgperm2s_14));
q_Namikas03{9} = 1e3*q_kgperm2s_14(~isnan(q_kgperm2s_14));
RunName_Namikas03{9} = '14';

% Calculate saltation flux and height
zq_Namikas03 = zeros(N_Namikas03,1); %e-folding height (m)
sigma_zq_Namikas03 = zeros(N_Namikas03,1); %e-folding height (m)
Q_fit_Namikas03 = zeros(N_Namikas03,1); %total flux (g/m/s)

for i=1:N_Namikas03
    [q0,zq,~,sigma_q0,sigma_zq,~] = qz_profilefit(q_Namikas03{i}, z_Namikas03{i});
    Q = q0*zq; %calculate total flux (g/m/s)
    zq_Namikas03(i) = zq;
    sigma_zq_Namikas03(i) = sigma_zq;
    Q_fit_Namikas03(i) = Q;
end

%remove outlier values
N_Namikas03 = length(Namikas03_GoodInd);
q_Namikas03 = q_Namikas03(Namikas03_GoodInd);
Q_pub_Namikas03 = Q_pub_Namikas03(Namikas03_GoodInd);
Q_fit_Namikas03 =Q_fit_Namikas03(Namikas03_GoodInd);
RunName_Namikas03 = RunName_Namikas03(Namikas03_GoodInd);
ust_Namikas03 = ust_Namikas03(Namikas03_GoodInd);
z_Namikas03 = z_Namikas03(Namikas03_GoodInd);
zq_Namikas03 = zq_Namikas03(Namikas03_GoodInd);
sigma_zq_Namikas03 = sigma_zq_Namikas03(Namikas03_GoodInd);


%% Load and process Namikas 2003 data
Farrell12_GoodInd = setdiff(1:N_Farrell12,Farrell12_OutlierInd); %information about which values to include in analysis
fileID_Farrell12 = fopen(filename_Farrell12,'r');
dataArray_Farrell12 = textscan(fileID_Farrell12, formatSpec_Farrell12, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID_Farrell12);

% Initialize cell arrays with all values
z_Farrell12 = cell(N_Farrell12,1); %trap height (m)
q_Farrell12 = cell(N_Farrell12,1); %flux (g/m^2/s)
zq_Farrell12 = zeros(N_Farrell12,1); %e-folding height (m)
sigma_zq_Farrell12 = zeros(N_Farrell12,1); %e-folding height uncertainty (m)
Q_fit_Farrell12 = zeros(N_Farrell12,1); %total flux (g/m/s)
RunName_Farrell12 = cell(N_Farrell12,1); %run name

% Allocate imported array to column variable names and cell arrays
for i = 1:N_Farrell12
    ind_profile = find(dataArray_Farrell12{1}==i);
    
    %compute sediment flux and trap heights
    m_profile = dataArray_Farrell12{4}(ind_profile); %mass (g)
    T_profile = dataArray_Farrell12{3}(ind_profile); %duration (s)
    H_profile = dataArray_Farrell12{6}(ind_profile); %trap height (m)
    W_profile = dataArray_Farrell12{7}(ind_profile); %trap width (m)
    q_Farrell12{i} = m_profile./(T_profile.*H_profile.*W_profile);
    z_bottom_profile = dataArray_Farrell12{5}(ind_profile)-H_profile; %bottom height of profile traps (m)
    
    %iterative fit to profile to optimize trap heights for fitting
    %[z_profile,q0,zq,Q,~,sigma_zq] = BSNE_profilefit(q_Farrell12{i}, z_bottom_profile, H_profile);
    [z_profile,q0,zq,Q,~,sigma_zq,~,~,~,~,~,z_profile_geomean,q0_geomean,zq_geomean,Q_geomean] = BSNE_profilefit(q_Farrell12{i}, z_bottom_profile, H_profile);
%     %start with guess of trap midpoint heights as arithmetic mean of traps
%     z_profile = z_bottom_profile + H_profile/2;
%     
%     %height-integrated flux from exponential fit
%     [q0,zq,sigma_q0,sigma_zq,qz_pred,sigma_qz_pred,sigma_logqz_pred] = qz_profilefit(q_Farrell12{i}, z_profile);
%     
%     %now that we have zq, redo calculation of trap heights
%     z_profile_old = z_profile; %document previous z-profile to see difference
%     z_profile = z_profile_calc(z_bottom_profile,H_profile,zq); %calculate new trap midpoint heights based on zq
%     z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
%     
%     %iterate until the z_profile_difference is minutely small
%     while(z_profile_difference>1e-8)
%         [q0,zq,~,sigma_q0,sigma_zq,~] = qz_profilefit(q_Farrell12{i}, z_profile); %height-integrated flux from exponential fit
%         z_profile_old = z_profile; %document previous z-profile to see difference
%         z_profile = z_profile_calc(z_bottom_profile,H_profile,zq); %calculate new trap midpoint heights based on zq
%         z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
%     end
    
    %add final list of profile heights to array
    z_Farrell12{i} = z_profile;
    
    %add fluxes to array
    zq_Farrell12(i) = zq; %e-folding height (m)
    sigma_zq_Farrell12(i) = sigma_zq; %e-folding height uncertainty (m)
    Q_fit_Farrell12(i) = q0*zq; %total flux (g/m/s)
    
    %add run name to list
    RunName_Farrell12{i} = int2str(i);
end

%remove outlier values
N_Farrell12 = length(Farrell12_GoodInd);
q_Farrell12 = q_Farrell12(Farrell12_GoodInd);
Q_fit_Farrell12 =Q_fit_Farrell12(Farrell12_GoodInd);
RunName_Farrell12 = RunName_Farrell12(Farrell12_GoodInd);
ust_Farrell12 = ust_Farrell12(Farrell12_GoodInd);
z_Farrell12 = z_Farrell12(Farrell12_GoodInd);
zq_Farrell12 = zq_Farrell12(Farrell12_GoodInd);
sigma_zq_Farrell12 = sigma_zq_Farrell12(Farrell12_GoodInd);

%% save useful data
save(SaveData_Path,'*Greeley96','*Namikas03','*Farrell12');