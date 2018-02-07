%% initialize
clearvars;

%% basic parameters
rho_Site = [1.16; 1.22; 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/Thresholds/'; %folder for threshold data
folder_GrainSizeData = '../../AnalysisData/GrainSize/'; %folder for grain size data
folder_Functions = '../Functions/'; %folder with functions
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for threshold data

%% load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% load data
load(strcat(folder_AnalysisData,'ThresholdAnalysisData')); %load analysis windows
load(strcat(folder_GrainSizeData,'MeanGrainSize')); %load grain size data

%% plotting info
PlotFont = 12;
PlotMarkers_Site = {'s','d','o'};
PlotColors_Site = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};

%% information about wind tunnel experiments
rho_windtunnel = 1.2; %kg/m^3 (assumes T~20C)

d50_Bagnold1937 = [0.08,0.24,0.475,0.65,1]; %fig 12?
d50_raw_Chepil1945 = 1e-3*[23.389402888469,75.4247813642512,129.482509138756,162.452181441327,228.055591884665,344.191392057723,514.301508399084,713.494200125,1022.88310619382,1585.3099603592,2487.41229824242]; %geom. mean. values from table 1
rho_p_Chepil1945 = [2.58, 2.07, 2.09, 1.96, 1.94, 1.91, 1.91, 1.8, 1.78, 1.74, 1.65]; %particle diameters for correction, table 1
d50_Chepil1945 = d50_raw_Chepil1945.*rho_p_Chepil1945/2.65; %computed corrected grain diameters based on densities (Fig. 5, Kok et al., 2012)
d50_IversenRasmussen1994 = [0.320]; %fig 7
d50_LiMcKennaNeuman2012 = [0.55]; %quoted in paper
d50_McKennaNeumanBedard2016 = [0.29, 0.29, 0.31, 0.35, 0.37, 0.40, 0.73]; %table 1

d90_Bagnold1937 = [0.08,0.24,0.475,0.65,1]; %does not list, but assume the same as d50 because he calls it "uniform sand"
d90_raw_Chepil1945 = [0.05,0.1,0.15,0.18,0.25,0.42,0.59,0.83,1.19,2.0,3.0]; %max values from table 1
d90_Chepil1945 = d90_raw_Chepil1945.*rho_p_Chepil1945/2.65; %computed corrected grain diameters based on densities (Fig. 5, Kok et al., 2012)
d90_IversenRasmussen1994 = [0.350]; %table 1, sample 4
d90_LiMcKennaNeuman2012 = [0.550]; %use the same value as d50, for "well-sorted coarse sand"
d90_McKennaNeumanBedard2016 = [0.54, 0.54, 0.54, 0.54, 0.6, 0.85, 1]; %table 1
%"fine" - use upper limit screen size for fine fraction
%"1" - use upper limit screen size for fine fraction
%"2" - use lower limit screen size for coarse fraction
%"3" - use lower limit screen size for coarse fraction
%"4" - use modal value for coarse fraction
%"5" - use 75th percentile value for coarse fraction
%"coarse" - use second large coarse bin

ustit_Bagnold1937 = [0.147201,0.192,0.249072,0.301934,0.378681]; %fig 12?
ustit_Chepil1945 = [0.252616,0.165425,0.129336,0.139619,0.185891,0.197834,0.250146,0.286604,0.368084,0.435565,0.500015]; %"maximal" impact threshold velocity, fig 3?
ustit_IversenRasmussen1994 = [0.21]; %fig 3
ustit_LiMcKennaNeuman2012 = [0.28]; %quoted in paper
uit_McKennaNeumanBedard2016 = [5.875, 6.125, 6.125, 6.125, 6.625, 7.125, 12.875]; %midpoints of 0.25 m/s ranges
ustu_ratio_McKennaNeumanBedard2016 = [0.054, 0.054, 0.052, 0.048, 0.045, 0.043, 0.043]; %estimated ratios of u* and u
ustit_McKennaNeumanBedard2016 = ustu_ratio_McKennaNeumanBedard2016.*uit_McKennaNeumanBedard2016; %calculate impact threshold shear velocities

tauit_Bagnold1937 = rho_windtunnel*ustit_Bagnold1937.^2;
tauit_Chepil1945 = rho_windtunnel*ustit_Chepil1945.^2;
tauit_IversenRasmussen1994 = rho_windtunnel*ustit_IversenRasmussen1994.^2;
tauit_LiMcKennaNeuman2012 = rho_windtunnel*ustit_LiMcKennaNeuman2012.^2;
tauit_McKennaNeumanBedard2016 = rho_windtunnel*ustit_McKennaNeumanBedard2016.^2;

%% information about field studies
folder_FieldData = '../../AnalysisData/FluxLaw/'; %folder for loading lit data
load(strcat(folder_FieldData,'LitData')); %Literature data

rho_field = 1.22; %kg/m^3 (assumes T~15C)

%get threshold stress for Greeley fit to flux versus stress
d90_Greeley96 = 1.2; %mm, value chosen based on "coaser fraction (~1200 um) occurs as a lag deposit on large ripples
tau_Greeley96 = rho_field.*ust_Greeley96.^2;
sigma_tau_Greeley96 = 2*rho_field.*ust_Greeley96.*sigma_ust_Greeley96;
[intercept, slope, sigma_intercept, ~, Q_pred_Greeley96, ~] = ... %perform linear fit
    linearfit(tau_Greeley96, Q_fit_Greeley96, sigma_Q_fit_Greeley96);
tauit_Greeley96 = -intercept/slope;
sigma_tauit_Greeley96 = sigma_intercept/slope;

%get threshold stress for Namikas fit to flux versus stress
d90_Namikas03 = 0.4175; %mm, value chosen based on "sorting value of 0.37 phi", and, phi_50 = 2, thus phi_84 = 1.26, and apply phi formula
tau_Namikas03 = rho_field.*ust_Namikas03.^2;
sigma_tau_Namikas03 = 2*rho_field.*ust_Namikas03.*sigma_ust_Namikas03;
[intercept, slope, sigma_intercept, ~, Q_pred_Namikas03, ~] = ... %perform linear fit
    linearfit(tau_Namikas03(2:end-1), Q_fit_Namikas03(2:(length(Q_fit_Namikas03)-1)), sigma_Q_fit_Namikas03(2:(length(Q_fit_Namikas03)-1))); %exclude first and last values from fit
tauit_Namikas03 = -intercept/slope;
sigma_tauit_Namikas03 = sigma_intercept/slope;

%plot fits
figure(1); clf;
subplot(1,2,1); hold on;
errorbar(tau_Greeley96,Q_fit_Greeley96,sigma_Q_fit_Greeley96,'x');
plot(tau_Greeley96,Q_pred_Greeley96);
plot(tauit_Greeley96+[-1 1]*sigma_tauit_Greeley96,[0 0],'k','LineWidth',2);
xlabel('shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
ylabel('saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
title('Greeley (1996)');
set(gca,'Box','On');
set(gca,'FontSize',PlotFont);

subplot(1,2,2); hold on;
errorbar(tau_Namikas03,Q_fit_Namikas03,sigma_Q_fit_Namikas03,'x');
plot(tau_Namikas03(2:end-1),Q_pred_Namikas03);
plot(tauit_Namikas03+[-1 1]*sigma_tauit_Namikas03,[0 0],'k','LineWidth',2);
xlabel('shear stress, $$\tau$$ (Pa)','Interpreter','Latex');
ylabel('saltation flux, $$Q$$ (g/m/s)','Interpreter','Latex');
title('Namikas (2003)');
set(gca,'Box','On');
set(gca,'FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 6]);
print([folder_Plots,'tauit_GreeleyNamikas.png'],'-dpng');

%% information about COMSALT simulations
ustit_COMSALT_mono = [0.2987525;0.2927475;0.2526775];
sigma_ustit_COMSALT_mono = [0.0015571421472257;0.0016788364025916;2.78253481559513E-4];
ustit_COMSALT_poly = [0.3454;0.3205675;0.28254];
sigma_ustit_COMSALT_poly = [0.00255532124529708;9.19506208063144E-4;3.53741902145216E-4];
tauit_COMSALT_mono = rho_Site.*ustit_COMSALT_mono.^2;
sigma_tauit_COMSALT_mono = 2*rho_Site.*ustit_COMSALT_mono.*sigma_ustit_COMSALT_mono;
tauit_COMSALT_poly = rho_Site.*ustit_COMSALT_poly.^2;
sigma_tauit_COMSALT_poly = 2*rho_Site.*ustit_COMSALT_poly.*sigma_ustit_COMSALT_poly;

%% full COMSALT simulation of impact threshold from Kok et al. (2012) fig 21a and file "Figure21a_Earth_COMSALT_uit.dat"
d_COMSALT_mono_full = 1e-3*[50;57;65;80;100;125;150;200;250;325;400;500;650;800;1000;1150;1350;1650;2000];
ustit_COMSALT_mono_full = [0.1914325;0.1587225;0.139529375;0.13898875;0.14131125;0.1490075;0.156135;0.178005;0.19810625;0.22675;0.25335375;0.28567;0.326585625;0.36428125;0.40729875;0.436010625;0.47370875;0.521250625;0.5729775];
tauit_COMSALT_mono_full = rho_windtunnel.*ustit_COMSALT_mono_full.^2;

%get equivalent values for Lapotre
ustit_Lapotre_d_COMSALT = zeros(size(d_COMSALT_mono_full));
tauit_Lapotre_d_COMSALT = zeros(size(d_COMSALT_mono_full));
for i = 1:length(tauit_Lapotre_d_COMSALT)
    ustit_Lapotre_d_COMSALT(i) = ustit_Lapotreetal2016(d_COMSALT_mono_full(i)*1e-3);
    tauit_Lapotre_d_COMSALT(i) = rho_windtunnel.*ustit_Lapotre_d_COMSALT(i).^2;
end

%% Get info about sites
N_Sites = length(Sites);

%% estimate tauit for d50
d50_m = d50_surface_site/1000;
sigma_d50_m = sigma_d50_surface_site/1000;
d50_plus_m = d50_m+sigma_d50_m;
d50_minus_m = d50_m-sigma_d50_m;

tauit_d50 = zeros(N_Sites,1);
sigma_tauit_d50 = zeros(N_Sites,1);
for i = 1:N_Sites
    ustit_d50 = ustit_Lapotreetal2016(d50_m(i));
    ustit_d50_plus = ustit_Lapotreetal2016(d50_plus_m(i));
    ustit_d50_minus = ustit_Lapotreetal2016(d50_minus_m(i));
    tauit_d50(i) = rho_Site(i)*ustit_d50^2;
    tauit_d50_plus = rho_Site(i)*ustit_d50_plus^2;
    tauit_d50_minus = rho_Site(i)*ustit_d50_minus^2;
    sigma_tauit_d50(i) = (tauit_d50_plus-tauit_d50_minus)/2;
end

%% estimate tauit for d90
d90_m = d90_surface_site/1000;
sigma_d90_m = sigma_d90_surface_site/1000;
d90_plus_m = d90_m+sigma_d90_m;
d90_minus_m = d90_m-sigma_d90_m;

tauit_d90 = zeros(N_Sites,1);
sigma_tauit_d90 = zeros(N_Sites,1);
for i = 1:N_Sites
    ustit_d90 = ustit_Lapotreetal2016(d90_m(i));
    ustit_d90_plus = ustit_Lapotreetal2016(d90_plus_m(i));
    ustit_d90_minus = ustit_Lapotreetal2016(d90_minus_m(i));
    tauit_d90(i) = rho_Site(i)*ustit_d90^2;
    tauit_d90_plus = rho_Site(i)*ustit_d90_plus^2;
    tauit_d90_minus = rho_Site(i)*ustit_d90_minus^2;
    sigma_tauit_d90(i) = (tauit_d90_plus-tauit_d90_minus)/2;
end

%% plot Lapotre et al (2016) u*it curve
figure(2); clf; hold on;
plot(d_COMSALT_mono_full, ustit_Lapotre_d_COMSALT,'r');
plot(d_COMSALT_mono_full, ustit_COMSALT_mono_full,'b--');
plot(d50_Bagnold1937,ustit_Bagnold1937,'bd','MarkerSize',5); %plot values - literature
plot(d50_Chepil1945,ustit_Chepil1945,'bo','MarkerSize',5); %plot values - literature
plot(d50_IversenRasmussen1994,ustit_IversenRasmussen1994,'b^','MarkerSize',5); %plot values - literature
plot(d50_LiMcKennaNeuman2012,ustit_LiMcKennaNeuman2012,'bs','MarkerSize',5); %plot values - literature
xlim([0.05 2]);
ylim([0.02 4]);
legend('Lapotre et al. (2016), Eq. S1','COMSALT','Bagnold (1937)','Chepil (1945)','Iversen and Rasmussen (1994)','Li and McKenna Neuman (2012)','Location','NorthEast');
xlabel('Grain diameter, $$d$$ (mm)','Interpreter','Latex');
ylabel('Impact threshold shear velocity, $$u_{*,it}$$ (m/s)','Interpreter','Latex');
set(gca,'xscale','log','yscale','log','Box','On','FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'ustit_Fig21a.png'],'-dpng');

%% plot comparing various estimates of tauit - d50
figure(3); clf;
subplot(1,2,1); hold on;
 
%plot values from full COMSALT simulation for monodisperse grains - Kok et al., 2012
plot(d_COMSALT_mono_full,tauit_COMSALT_mono_full,'b--');

for i = 1:N_Sites
    %tauit error bars
    plot(d50_surface_site(i)*[1 1],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'-k','LineWidth',1); %plot TFEM value
    plot(d50_surface_site(i)*[1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],':k','LineWidth',1); %plot intercept SE dummy values for legend
    %plot(d50_surface_site(i)*[1 1],tauit_d50(i)+sigma_tauit_d50(i)*[-1 1],'--k','LineWidth',1); %plot it SE dummy values for legend
 
    %tauit values and error bars - COMSALT
    plot(d50_surface_site(i),tauit_COMSALT_mono(i),'kv','MarkerSize',10); %plot COMSALT values for homogeneous bed
    plot(d50_surface_site(i)*[1 1],tauit_COMSALT_mono(i)+sigma_tauit_COMSALT_mono(i)*[-1 1],'k-','LineWidth',1); %plot COMSALT values for homogeneous bed
    
    %tauit values
    plot(d50_surface_site(i)*[1 1],tauit_all(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot intercept average values
    plot(d50_surface_site(i)*[1 1],tauit_intercept(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot intercept average values
    %plot(d50_surface_site(i)*[1 1],tauit_d50(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot it average values
   
    %tauit error bars again - with color
    plot(d50_surface_site(i)*[1 1],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    plot(d50_surface_site(i)*[1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],':','Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    %plot(d50_surface_site(i)*[1 1],tauit_d50(i)+sigma_tauit_d50(i)*[-1 1],'--','Color',PlotColors_Site{i},'LineWidth',1); %plot it SE values
    
    %d50 error bars
    plot(d50_surface_site(i)+sigma_d50_surface_site(i)*[-1 1],tauit_all(i)*[1 1],'Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    plot(d50_surface_site(i)+sigma_d50_surface_site(i)*[-1 1],tauit_intercept(i)*[1 1],':','Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    %plot(d50_surface_site(i)+sigma_d50_surface_site(i)*[-1 1],tauit_d50(i)*[1 1],'--','Color',PlotColors_Site{i},'LineWidth',1); %plot it SE values
end

%plot values - literature
plot(d50_Bagnold1937,tauit_Bagnold1937,'bd','MarkerSize',5); 
plot(d50_Chepil1945,tauit_Chepil1945,'bo','MarkerSize',5); 
plot(d50_IversenRasmussen1994,tauit_IversenRasmussen1994,'b^','MarkerSize',5); 
plot(d50_LiMcKennaNeuman2012,tauit_LiMcKennaNeuman2012,'bs','MarkerSize',5); 
plot(d50_McKennaNeumanBedard2016,tauit_McKennaNeumanBedard2016,'bv','MarkerSize',5); 

%plot values - field
plot(d50_Greeley96,tauit_Greeley96,'r<','MarkerSize',5);
plot(d50_Namikas03,tauit_Namikas03,'r>','MarkerSize',5);

%organize plot
%legend('Kok et al. (2012)','TFEM fit','flux law fit','pred. by d_{50}, Lapotre et al. (2016)','COMSALT (monodisperse)','Location','SouthEast');
legend('Kok et al. (2012)','TFEM fit','flux law fit','COMSALT (monodisperse)','Location','SouthEast');
xlabel('median grain diameter, $$d_{50}$$ (mm)','interpreter','latex');
ylabel('impact threshold stress, $$\tau_{it}$$ (Pa)','interpreter','latex');
xlim([0.05 2]);
ylim([0.01 0.4]);
set(gca,'xscale','log','yscale','log');
set(gca,'Box','On');
set(gca,'FontSize',PlotFont);

%% plot comparing various estimates of tauit - d90
subplot(1,2,2); hold on;

%plot values - field
plot(d90_Greeley96,tauit_Greeley96,'r<','MarkerSize',5);
plot(d90_Namikas03,tauit_Namikas03,'r>','MarkerSize',5);

%plot values - literature
plot(d90_Bagnold1937,tauit_Bagnold1937,'bd','MarkerSize',5);
plot(d90_Chepil1945,tauit_Chepil1945,'bo','MarkerSize',5);
plot(d90_IversenRasmussen1994,tauit_IversenRasmussen1994,'b^','MarkerSize',5);
plot(d90_LiMcKennaNeuman2012,tauit_LiMcKennaNeuman2012,'bs','MarkerSize',5); 
plot(d90_McKennaNeumanBedard2016,tauit_McKennaNeumanBedard2016,'bv','MarkerSize',5); 

for i = 1:N_Sites
    %tauit error bars - d90 prediction
    %plot(d90_surface_site(i)*[1 1],tauit_d90(i)+sigma_tauit_d90(i)*[-1 1],'-.k','LineWidth',1); %plot it SE dummy values for legend
 
    %tauit values and error bars - COMSALT
    plot(d90_surface_site(i),tauit_COMSALT_poly(i),'k^','MarkerSize',10); %plot COMSALT values for heterogeneous bed
    plot(d90_surface_site(i)*[1 1],tauit_COMSALT_poly(i)+sigma_tauit_COMSALT_poly(i)*[-1 1],'k-','LineWidth',1); %plot COMSALT values for heterogeneous bed
  
    %tauit values
    plot(d90_surface_site(i)*[1 1],tauit_all(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot intercept average values
    plot(d90_surface_site(i)*[1 1],tauit_intercept(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot intercept average values
    %plot(d90_surface_site(i)*[1 1],tauit_d90(i),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot it average values
   
    %tauit error bars again - with color
    plot(d90_surface_site(i)*[1 1],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    plot(d90_surface_site(i)*[1 1],tauit_intercept(i)+sigma_tauit_intercept(i)*[-1 1],':','Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    %plot(d90_surface_site(i)*[1 1],tauit_d90(i)+sigma_tauit_d90(i)*[-1 1],'-.','Color',PlotColors_Site{i},'LineWidth',1); %plot it SE values
    
    %d90 error bars
    plot(d90_surface_site(i)+sigma_d90_surface_site(i)*[-1 1],tauit_all(i)*[1 1],'Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    plot(d90_surface_site(i)+sigma_d90_surface_site(i)*[-1 1],tauit_intercept(i)*[1 1],':','Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
    %plot(d90_surface_site(i)+sigma_d90_surface_site(i)*[-1 1],tauit_d90(i)*[1 1],'-.','Color',PlotColors_Site{i},'LineWidth',1); %plot it SE values
end

%plot values from full COMSALT simulation for monodisperse grains - Kok et al., 2012
plot(d_COMSALT_mono_full,tauit_COMSALT_mono_full,'b--');

%organize plot
%legend('Greeley et al. (1996)','Namikas (2003)','Bagnold (1937)','Chepil (1945)','Iversen & Rasmussen (1994)','Li & McKenna Neuman (2012)','McKenna Neuman & Bedard (2016)','pred. by d_{90}, Lapotre et al. (2016)','COMSALT (polydisperse)','Location','NorthWest');
legend('Greeley et al. (1996)','Namikas (2003)','Bagnold (1937)','Chepil (1945)','Iversen & Rasmussen (1994)','Li & McKenna Neuman (2012)','McKenna Neuman & Bedard (2016)','COMSALT (polydisperse)','Location','SouthEast');
xlabel('90th percentile grain diameter, $$d_{90}$$ (mm)','interpreter','latex');
ylabel('impact threshold stress, $$\tau_{it}$$ (Pa)','interpreter','latex');
xlim([0.05 2]);
ylim([0.01 0.4]);
set(gca,'xscale','log','yscale','log');
set(gca,'Box','On');
set(gca,'FontSize',PlotFont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 6]);
print([folder_Plots,'tauit_d50_d90.png'],'-dpng');

%% plot u*it comparison
figure(5); clf; hold on;

%plot literature values
plot(d50_Bagnold1937,ustit_Bagnold1937,'kd','MarkerSize',5); %plot values - literature
plot(d50_Chepil1945,ustit_Chepil1945,'ko','MarkerSize',5); %plot values - literature
plot(d50_IversenRasmussen1994,ustit_IversenRasmussen1994,'k^','MarkerSize',5); %plot values - literature
plot(d50_LiMcKennaNeuman2012,ustit_LiMcKennaNeuman2012,'ks','MarkerSize',5); %plot values - literature

%plot field values
for i = 1:N_Sites       
    %tauit versus d50 values
    plot(d50_surface_site(i),sqrt(tauit_all(i)/rho_Site(i)),PlotMarkers_Site{i},'Color',PlotColors_Site{i},'MarkerSize',10); %plot intercept average values
end

% %plot field error bars
% for i = 1:N_Sites
%     %tauit error bars
%     plot(d50_surface_site(i)*[1 1],tauit_all(i)+sigma_tauit_all(i)*[-1 1],'Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
%     
%     %d50 error bars
%     plot(d50_surface_site(i)+sigma_d50_surface_site(i)*[-1 1],tauit_all(i)*[1 1],'Color',PlotColors_Site{i},'LineWidth',1); %plot intercept SE values
% end

xlim([0.05 2]);
ylim([1e-1 1e0]);
legend('Bagnold (1937)','Chepil (1945)','Iversen and Rasmussen (1994)','Li and McKenna Neuman (2012)',...
    'Jericoacoara','Rancho Guadalupe','Oceano',...
    'Location','NorthWest');
xlabel('Median grain diameter, $$d_{50}$$ (mm)','Interpreter','Latex');
ylabel('Impact threshold shear velocity, $$u_{*,it}$$ (m/s)','Interpreter','Latex');
set(gca,'xscale','log','yscale','log','Box','On','FontSize',PlotFont);
set(gca,'xticklabel',{'0.1','1.0'})
set(gca,'yticklabel',{'0.1','1.0'})

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4.5]);
print([folder_Plots,'ustit_d50.png'],'-dpng');

