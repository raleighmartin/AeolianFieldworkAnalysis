%% SCRIPT TO IMPORT GREELEY (1996), NAMIKAS (2003), AND FARRELL (2012) FLUX VALUES FROM TEXT FILE
clearvars; %initialize

%suppress figures
set(0,'DefaultFigureVisible', 'off');

%% constants
rho_a = 1.23; %air density kg/m^3

%% filepath information
folder_LiteratureData = '../../../Google Drive/Data/Literature/'; %location where Namikas and Greeley data are stored
filename_Greeley96 = [folder_LiteratureData,'Greeley96.dat'];
filename_Namikas03 = [folder_LiteratureData,'Namikas03verMF.dat'];
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_AnalysisData = '../AnalysisData/'; %folder for storing outputs of this analysis
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots
SaveData_Path = strcat(folder_AnalysisData,'GreeleyNamikasData'); %folder for saving data

%% manually enter basic info
N_Greeley96 = 10;
N_Namikas03 = 9;
ust_Greeley96 = [0.47; 0.51; 0.49; 0.31; 0.35; 0.31; 0.54; 0.49; 0.41; 0.42]; %u* values (m/s) used in Kok et al. (2008)
ust_Namikas03 = [0.32; 0.27; 0.32; 0.30; 0.38; 0.37; 0.38; 0.47; 0.63]; %u* values (m/s) from paper
ust_Farrell12 = [0.54; 0.47; 0.53; 0.49; 0.50; 0.50; 0.47; 0.45; 0.51; 0.48; 0.41; 0.50; 0.50; 0.49]; %u* values (m/s) from paper
Q_fit_Farrell12 = [30.16; 29.06; 19.37; 22.49; 19.48; 16.59; 15.10; 14.15; 16.44; 25.46; 19.94; 20.33; 27.51; 22.89]; %Q values (g/m/s) estimated from paper
zbar_Farrell12 = [0.113; 0.088; 0.088; 0.077; 0.087; 0.075; 0.069; 0.085; 0.081; 0.095; 0.083; 0.118; 0.123; 0.123]; %zQ values (m) estimated from paper 
Q_pub_Greeley96 = [NaN; NaN; NaN; NaN; 16.1; 2.7; 19; 19; 10.4; 10.4]; %values from Table 1 in paper
Q_pub_Namikas03 = [0.35; 0.76; 2.73; 2.12; 5.61; 5.86; 5.57; 11.80; 29.11]; %values from Table 1 in paper
d50_Greeley96 = 0.230; %median grain size (mm) value from paper
d50_Namikas03 = 0.250; %median grain size (mm) value from paper
d50_Farrell12 = 0.61; %use median grain size (mm) from Jeri deployment since they don't report it

%% give information about outliers to exclude
Greeley96_OutlierInd = 4;
Namikas03_OutlierInd = 9;
Greeley96_GoodInd = setdiff(1:N_Greeley96,Greeley96_OutlierInd);
Namikas03_GoodInd = setdiff(1:N_Namikas03,Namikas03_OutlierInd);

%% Initialize processing variables
delimiter = '\t';
startRow = 2;

%% Format string for each line of text:
formatSpec_Greeley96 = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
formatSpec_Namikas03 = '%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text files
fileID_Greeley96 = fopen(filename_Greeley96,'r');
fileID_Namikas03 = fopen(filename_Namikas03,'r');

%% Read columns of data according to format strings
dataArray_Greeley96 = textscan(fileID_Greeley96, formatSpec_Greeley96, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
dataArray_Namikas03 = textscan(fileID_Namikas03, formatSpec_Namikas03, 'Delimiter', delimiter, 'TreatAsEmpty' ,{'--'},'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text files
fclose(fileID_Greeley96);
fclose(fileID_Namikas03);

%% Process Greeley data

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
zbar_Greeley96 = zeros(N_Greeley96,1); %(m)
Q_fit_Greeley96 = zeros(N_Greeley96,1); %(g/m/s)

for i=1:N_Greeley96
    [q0,zbar] = qz_profilefit(q_Greeley96{i}, z_Greeley96{i}, 0, 0);
    Q = q0*zbar;
    zbar_Greeley96(i) = zbar;
    Q_fit_Greeley96(i) = Q;
    
    figure(i); clf; hold on;
    plot(q_Greeley96{i},z_Greeley96{i},'x')
    z_range = linspace(0, max(z_Greeley96{i}),50);
    q_pred = q0*exp(-z_range./zbar);
    plot(q_pred,z_range);
    set(gca,'xscale','log');
    set(gca,'FontSize',16);
    xlabel('q (g m^{-2} s^{-1})','FontSize',16);
    ylabel('z (m)','FontSize',16);
    title(['Run = ',RunName_Greeley96{i},', u_{*} = ',num2str(ust_Greeley96(i)),' m/s, Q = ',num2str(Q),' g/m/s, z_{q} = ',num2str(zbar),' m']);
    print([folder_Plots,'Greeley96_qzprofile_',RunName_Greeley96{i},'.png'],'-dpng');
end

%% Process Namikas data

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
zbar_Namikas03 = zeros(N_Namikas03,1); %e-folding height (m)
Q_fit_Namikas03 = zeros(N_Namikas03,1); %total flux (g/m/s)

for i=1:N_Namikas03
    [q0,zbar] = qz_profilefit(q_Namikas03{i}, z_Namikas03{i});
    Q = q0*zbar; %calculate total flux (g/m/s)
    zbar_Namikas03(i) = zbar;
    Q_fit_Namikas03(i) = Q;
    
    figure(i+N_Namikas03); clf; hold on;
    plot(q_Namikas03{i},z_Namikas03{i},'x')
    z_range = linspace(0, max(z_Namikas03{i}),50);
    q_pred = q0*exp(-z_range./zbar);
    plot(q_pred,z_range);
    set(gca,'xscale','log');
    set(gca,'FontSize',16);
    xlabel('q (g m^{-2} s^{-1})','FontSize',16);
    ylabel('z (m)','FontSize',16);
    title(['Run ',RunName_Namikas03{i},', u_{*} = ',num2str(ust_Namikas03(i)),' m/s, Q = ',num2str(Q),' g/m/s, z_{q} = ',num2str(zbar),' m']);
    print([folder_Plots,'Namikas03_qzprofile_',RunName_Namikas03{i},'.png'],'-dpng');
end

%plot zbar(e) versus u*
figure(20); clf;
plot(ust_Greeley96,zbar_Greeley96,'r^',...
    ust_Namikas03,zbar_Namikas03,'gd',...
    'MarkerSize',10);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('z_{q,e} (m)','FontSize',16);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_ust_zsalt_e.png'],'-dpng');

%plot z50 versus u*
figure(21); clf;
plot(ust_Greeley96,zbar_Greeley96*log(2),'r^',...
    ust_Namikas03,zbar_Namikas03*log(2),'gd',...
    'MarkerSize',10);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('z_{q,50} (m)','FontSize',16);
ylim([0.03 0.04]);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthEast');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_ust_zsalt_50.png'],'-dpng');

%plot Q versus u*
figure(22); clf;
plot(ust_Greeley96,Q_fit_Greeley96,'r^',...
    ust_Namikas03,Q_fit_Namikas03,'gd',...
    'MarkerSize',10);
xlabel('u_{*} (m/s)','FontSize',16);
ylabel('Q (g/m^2/s)','FontSize',16);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_ust_Q.png'],'-dpng');

%% Remove outlier value(s) from data
N_Greeley96 = length(Greeley96_GoodInd);
q_Greeley96 = q_Greeley96{Greeley96_GoodInd};
Q_pub_Greeley96 = Q_pub_Greeley96(Greeley96_GoodInd);
Q_fit_Greeley96 =Q_fit_Greeley96(Greeley96_GoodInd);
RunName_Greeley96 = RunName_Greeley96(Greeley96_GoodInd);
ust_Greeley96 = ust_Greeley96(Greeley96_GoodInd);
z_Greeley96 = z_Greeley96{Greeley96_GoodInd};
zbar_Greeley96 = zbar_Greeley96(Greeley96_GoodInd);

N_Namikas03 = length(Namikas03_GoodInd);
q_Namikas03 = q_Namikas03{Namikas03_GoodInd};
Q_pub_Namikas03 = Q_pub_Namikas03(Namikas03_GoodInd);
Q_fit_Namikas03 =Q_fit_Namikas03(Namikas03_GoodInd);
RunName_Namikas03 = RunName_Namikas03(Namikas03_GoodInd);
ust_Namikas03 = ust_Namikas03(Namikas03_GoodInd);
z_Namikas03 = z_Namikas03{Namikas03_GoodInd};
zbar_Namikas03 = zbar_Namikas03(Namikas03_GoodInd);

%% calculate shear stresses
tau_Greeley96 = rho_a*ust_Greeley96.^2;
tau_Namikas03 = rho_a*ust_Namikas03.^2;
tau_Farrell12 = rho_a*ust_Farrell12.^2;

%% fit to get threshold
[a, b, ~, ~, ~, ~] = linearfit(tau_Greeley96, Q_fit_Greeley96);
tauit_Greeley96 = -a/b;
tauex_Greeley96 = tau_Greeley96-tauit_Greeley96;
tauratio_Greeley96 = tauex_Greeley96./tau_Greeley96;
[a, b, ~, ~, ~, ~] = linearfit(tau_Namikas03, Q_fit_Namikas03);
tauit_Namikas03 = -a/b;
tauex_Namikas03 = tau_Namikas03-tauit_Namikas03;
tauratio_Namikas03 = tauex_Namikas03./tau_Namikas03;
%assume values from Jeri deployment for Namikas
tauit_Farrell12 = 0.1317;
tauex_Farrell12 = tau_Farrell12-tauit_Farrell12;
tauratio_Farrell12 = tauex_Farrell12./tau_Farrell12;

%% plot Q versus tau
figure(23); clf;
plot(tau_Greeley96,Q_fit_Greeley96,'r^',...
    tau_Namikas03,Q_fit_Namikas03,'gd',...
    'MarkerSize',10);
xlabel('\tau (Pa)','FontSize',16);
ylabel('Q (g/m^2/s)','FontSize',16);
set(gca,'FontSize',16);
h_legend = legend({'Greeley (1996)','Namikas (2003)'},'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'GreeleyNamikas_tau_Q.png'],'-dpng');

%% save useful data
save(SaveData_Path,'*Greeley96','*Namikas03','*Farrell12');

%re-allow figures
set(0,'DefaultFigureVisible', 'on');