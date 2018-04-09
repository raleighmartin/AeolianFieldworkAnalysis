%% SCRIPT TO IMPORT GREELEY (1996), NAMIKAS (2003), AND FARRELL (2012) FROM TEXT AND XLS FILES
clearvars; %initialize


%% assumptions
ust_relerr = 0.1; %10% relative error for u*
q_relerr = 0.1; %10% relative error for q
z_relerr = 0.1; %10% relative error for z

%% filepath information
folder_LiteratureData = '../../AnalysisData/Literature/'; %location where Namikas and Greeley data are stored
folder_SaveData = '../../AnalysisData/Literature/'; %folder for storing outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions
SaveData_Path = strcat(folder_SaveData,'LitData'); %path for saving data

%% processing variables
delimiter = '\t';
startRow = 2;


%% manually enter basic info - flux data

% Greeley et al. 1996
filename_Greeley96 = [folder_LiteratureData,'Greeley96.dat'];
formatSpec_Greeley96 = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
N_Greeley96 = 10;
ust_Greeley96 = [0.47; 0.51; 0.49; 0.31; 0.35; 0.31; 0.54; 0.49; 0.41; 0.42]; %u* values (m/s) used in Kok and Renno (2008)
sigma_ust_Greeley96 = ust_Greeley96*ust_relerr; %calculate u* uncertainty based on assumed relative error
Q_pub_Greeley96 = [NaN; NaN; NaN; NaN; 16.1; 2.7; 19; 19; 10.4; 10.4]; %values from Table 1 in paper
d50_Greeley96 = 0.230; %"modal" grain size (mm) value from paper (p. 10)
%Greeley96_OutlierInd = [4];
Greeley96_OutlierInd = []; %indices of outliers for Greeley data

% Namikas 2003
filename_Namikas03 = [folder_LiteratureData,'Namikas03verMF.dat'];
formatSpec_Namikas03 = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
N_Namikas03 = 9;
ust_Namikas03 = [0.32; 0.27; 0.32; 0.30; 0.38; 0.37; 0.38; 0.47; 0.63]; %u* values (m/s) from paper
sigma_ust_Namikas03 = ust_Namikas03*ust_relerr; %calculate u* uncertainty based on assumed relative error
Q_pub_Namikas03 = [0.35; 0.76; 2.73; 2.12; 5.61; 5.86; 5.57; 11.80; 29.11]; %values from Table 1 in paper
d50_Namikas03 = 0.250; %median grain size (mm) value from paper
%Namikas03_OutlierInd = [9];
Namikas03_OutlierInd = []; %indices of outliers for Namikas data

% Namikas 2006 (grain-size data associated with Namikas 2003 study)
filename_Namikas06_small = [folder_LiteratureData,'Namikas03_SizeResolvedVMF_SmallUfr.dat']; %data for small u*
filename_Namikas06_med = [folder_LiteratureData,'Namikas03_SizeResolvedVMF_MedUfr.dat']; %data for medium u*
filename_Namikas06_large = [folder_LiteratureData,'Namikas03_SizeResolvedVMF_LargeUfr.dat']; %data for large u*
formatSpec_Namikas06 = '%f%f%f%f%f%f%f%[^\n\r]';
ind_Namikas03_small = 1:3;
ind_Namikas03_med = 4:7;
ind_Namikas03_large = 8:9;
ust_Namikas06 = [0.30, 0.36, 0.55]; %shear velocity (m/s) - mean for each set of indices
z_bottom_Namikas06 = 1e-2*[0 1 2 4 7 13]'; %bottom of trap (m)
z_top_Namikas06 = 1e-2*[1 2 4 7 13 33]'; %top of trap (m)
H_Namikas06 = z_top_Namikas06-z_bottom_Namikas06; %heights of traps (m)
N_z_Namikas06 = length(z_bottom_Namikas06); %number of heights
d_lower_Namikas06 = [0.105 0.125 0.149 0.177 0.21 0.25 0.297 0.354 0.42 0.5]; %lower edge of grain size bins (mm)
d_upper_Namikas06 = [0.125 0.149 0.177 0.21 0.25 0.297 0.354 0.42 0.5 0.564]; %upper edge of grain size bins (mm)
% note - upper edge of last bin is assumed but probably there is no
% justification. to be safe, don't include this in any analyses

% Farrell 2012
filename_Farrell12 = [folder_LiteratureData,'FarrellEtAl2012.dat'];
formatSpec_Farrell12 = '%f%f%f%f%f%f%f%f%[^\n\r]';
N_Farrell12 = 14;
N_z_Farrell12 = 6; %number of heights for Farrell profiles
ust_Farrell12 = [0.54; 0.47; 0.53; 0.49; 0.50; 0.50; 0.47; 0.45; 0.51; 0.48; 0.41; 0.50; 0.50; 0.49]; %u* values (m/s) from paper
sigma_ust_Farrell12 = ust_Farrell12*ust_relerr; %calculate u* uncertainty based on assumed relative error
d50_Farrell12 = NaN; %median grain diameter not provided in paper
Farrell12_OutlierInd = [1,12:14];
%exclude first sample because of shorter time interval
%exclude last three because they are from coarser CSFC (as opposed to the rest of the samples from finer CSFF)

%% Load and process Greeley et al. 1996 data
GoodInd_Greeley96 = setdiff(1:N_Greeley96,Greeley96_OutlierInd); %information about which values to include in analysis
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

% Assign uncertainties
sigma_q_Greeley96 = cell(N_Greeley96,1);
sigma_z_Greeley96 = cell(N_Greeley96,1);
for i = 1:N_Greeley96
    sigma_q_Greeley96{i} = q_Greeley96{i}*q_relerr;
    sigma_z_Greeley96{i} = z_Greeley96{i}*z_relerr;
end


%% Load and process Namikas 2003 data
GoodInd_Namikas03 = setdiff(1:N_Namikas03,Namikas03_OutlierInd); %information about which values to include in analysis
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

% Assign uncertainties
sigma_q_Namikas03 = cell(N_Namikas03,1);
sigma_z_Namikas03 = cell(N_Namikas03,1);
for i = 1:N_Namikas03
    sigma_q_Namikas03{i} = q_Namikas03{i}*q_relerr;
    sigma_z_Namikas03{i} = z_Namikas03{i}*z_relerr;
end


%% Load and process Namikas 2006 grain-size data
fileID_Namikas06_small = fopen(filename_Namikas06_small,'r');
dataArray_Namikas06_small = textscan(fileID_Namikas06_small, formatSpec_Namikas06, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow, 'ReturnOnError', false);
fclose(fileID_Namikas06_small);
fileID_Namikas06_med = fopen(filename_Namikas06_med,'r');
dataArray_Namikas06_med = textscan(fileID_Namikas06_med, formatSpec_Namikas06, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow, 'ReturnOnError', false);
fclose(fileID_Namikas06_med);
fileID_Namikas06_large = fopen(filename_Namikas06_large,'r');
dataArray_Namikas06_large = textscan(fileID_Namikas06_large, formatSpec_Namikas06, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow, 'ReturnOnError', false);
fclose(fileID_Namikas06_large);

% get grain sizes
d_Namikas06 = 1e-3*dataArray_Namikas06_small{1}'; %size in mm
N_d_Namikas06 = length(d_Namikas06); %number of grain size bins

% get fluxes
qi_Namikas06_small = zeros(N_d_Namikas06,N_z_Namikas06);
qi_Namikas06_med = zeros(N_d_Namikas06,N_z_Namikas06);
qi_Namikas06_large = zeros(N_d_Namikas06,N_z_Namikas06);
for i = 1:N_z_Namikas06
    qi_Namikas06_small(:,i) = dataArray_Namikas06_small{i+1}; %height-specific size-specific flux - small u* (kg/m^2/s)
    qi_Namikas06_med(:,i) = dataArray_Namikas06_med{i+1}; %height-specific size-specific flux - med u* (kg/m^2/s)
    qi_Namikas06_large(:,i) = dataArray_Namikas06_large{i+1}; %height-specific size-specific flux - large u* (kg/m^2/s)
end

%combine into single cell array - transpose so that rows are grain sizes and columns are heights
qi_Namikas06 = {qi_Namikas06_small',...
    qi_Namikas06_med',...
    qi_Namikas06_large'};
N_Namikas06 = length(qi_Namikas06); %number of profiles

% Assign uncertainties
sigma_qi_Namikas06 = cell(N_Namikas06,1); %uncertainty in saltation flux
for i = 1:N_Namikas06
    sigma_qi_Namikas06{i} = qi_Namikas06{i}*q_relerr;
end
sigma_z_bottom_Namikas06 = z_bottom_Namikas06*z_relerr; %uncertainty in trap height


%% Load and process Farrell 2012 data
GoodInd_Farrell12 = setdiff(1:N_Farrell12,Farrell12_OutlierInd); %information about which values to include in analysis
fileID_Farrell12 = fopen(filename_Farrell12,'r');
dataArray_Farrell12 = textscan(fileID_Farrell12, formatSpec_Farrell12, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID_Farrell12);

% Initialize cell arrays with all values
z_bottom_Farrell12 = cell(N_Farrell12,1); %trap bottom height (m)
q_Farrell12 = cell(N_Farrell12,1); %flux (g/m^2/s)
H_Farrell12 = cell(N_Farrell12,1); %trap dimension height (m)
dbar_airborne_Farrell12 = cell(N_Farrell12,1); %mean height of grains in trap (mm)
%zq_Farrell12 = zeros(N_Farrell12,1); %e-folding height (m)
%sigma_zq_Farrell12 = zeros(N_Farrell12,1); %e-folding height uncertainty (m)
%Q_fit_Farrell12 = zeros(N_Farrell12,1); %total flux (g/m/s)
%sigma_Q_fit_Farrell12 = zeros(N_Farrell12,1); %total flux (g/m/s)
RunName_Farrell12 = cell(N_Farrell12,1); %run name


% Allocate imported array to column variable names and cell arrays
for i = 1:N_Farrell12
    ind_profile = find(dataArray_Farrell12{1}==i);
    
    %compute sediment flux and trap heights
    m_profile = dataArray_Farrell12{4}(ind_profile); %mass (g)
    T_profile = dataArray_Farrell12{3}(ind_profile); %duration (s)
    W_profile = dataArray_Farrell12{7}(ind_profile); %trap width (m)
    H_Farrell12{i} = dataArray_Farrell12{6}(ind_profile); %trap height (m)
    dbar_airborne_Farrell12{i} = dataArray_Farrell12{8}(ind_profile); %mean particle diameter (mm)
    z_bottom_Farrell12{i} = dataArray_Farrell12{5}(ind_profile)-H_Farrell12{i}; %bottom height of profile traps (m)
    q_Farrell12{i} = m_profile./(T_profile.*W_profile.*H_Farrell12{i});
    
    %add run name to list
    RunName_Farrell12{i} = int2str(i);
end

% Assign uncertainties
sigma_q_Farrell12 = cell(N_Farrell12,1);
sigma_z_bottom_Farrell12 = cell(N_Farrell12,1);
for i = 1:N_Farrell12
    sigma_q_Farrell12{i} = q_Farrell12{i}*q_relerr;
    sigma_z_bottom_Farrell12{i} = z_bottom_Farrell12{i}*z_relerr;
end


%% save useful data
save(SaveData_Path,'*Greeley96','*Namikas03','*Namikas06','*Farrell12');