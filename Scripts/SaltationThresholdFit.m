%% SEPARATE AND PLOT RESULTS BY GRAIN SIZE BINS, ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;

%% INPUT INFORMATION

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_Plots = '../PlotOutput/SaltationFlux/'; %folder for plots

%load stress flux window data
load(strcat(folder_ProcessedData,'StressFluxWindows'));

% %load stress flux continuous window data
% load(strcat(folder_ProcessedData,'StressFluxContinuousWindows'));

%get info about sites, set markers for plotting
N_Sites = length(Sites);
Markers = {'bx','ro','gv'};

%set parameter values
rho_a = 1.23; %air density kg/m^3
rho_s = 2650; %particle density kg/m^3
kappa = 0.39; %von Karman parameter
g = 9.8; %gravity m/s^2

%set bin edges
d50_bins_min = [0.33, 0.39, 0.51];
d50_bins_max = [0.36, 0.44, 0.54];
N_bins = length(d50_bins_min);

%% CALCULATIONS

%initialize list of calculated values for bins
tauThr_bins = zeros(N_bins,1);
ustThr_bins = zeros(N_bins,1);
d50_threshold_bins = zeros(N_bins,1);

%create cell array of thresholds for each site
ustThr_Site = cell(N_Sites,1);
tauThr_Site = cell(N_Sites,1);

%go through each bin
for i = 1:N_bins

    %initialize lists for bins
    tauRe_bin = [];
    ustRe_bin = [];
    Q_bin = [];
    d50_bin = [];
    
    %go through each site
    for j = 1:N_Sites
                
        %get values for this bin for this site
        ind_bin_Site = find((d50_list{j}>=d50_bins_min(i))&(d50_list{j}<=d50_bins_max(i)));
        ind_bin_Site = intersect(ind_bin_Site,find(Q_list{j}>0)); %include only indices with positive Q
        
        %if bin contains values for this Site, get values and add to overall list
        if ~isempty(ind_bin_Site)
            tauRe_bin_Site = tauRe_cal_list{j}(ind_bin_Site);
            ustRe_bin_Site = ustRe_cal_list{j}(ind_bin_Site);
            Q_bin_Site = Q_list{j}(ind_bin_Site);
            d50_bin_Site = d50_list{j}(ind_bin_Site);
            
            %add to lists for bins
            tauRe_bin = [tauRe_bin; tauRe_bin_Site];
            ustRe_bin = [ustRe_bin; ustRe_bin_Site];
            Q_bin = [Q_bin; Q_bin_Site];
            d50_bin = [d50_bin; d50_bin_Site];
        end
    end
    
    %for each bin, fit to determine threshold
    P_linear = polyfit(tauRe_bin,Q_bin,1);
    tauThr = -P_linear(2)/P_linear(1);
    ustThr = sqrt(tauThr/rho_a);
    
    %add to lists for each bin
    tauThr_bins(i) = tauThr;
    ustThr_bins(i) = ustThr;
    d50_threshold_bins(i) = mean(d50_bin);
    d50_print = [int2str(round(1000*mean(d50_bin))),'um'];
    
    %get values for plotting fit
    tauRe_fit = linspace(tauThr,max(tauRe_bin),50);
    ustRe_fit = sqrt(tauRe_fit/rho_a);
    Q_fit = polyval(P_linear,tauRe_fit);
    
    %plot flux versus shear stress
    figure(i*2-1); clf;
    plot(tauRe_bin,Q_bin, 'o'); hold on;
    plot(tauRe_fit,Q_fit);
    title(['d = ',d50_print],'FontSize',16);
    xlabel('\tau_{Re} (Pa)','FontSize',16);
    ylabel('Q (g/m/s)','FontSize',16);
    ylim([0 35]);
    set(gca,'FontSize',16);
    print([folder_Plots,'FluxTau_Wenglor_d',d50_print,'.png'],'-dpng');

    %plot flux versus shear velocity
    figure(i*2); clf;
    plot(ustRe_bin,Q_bin, 'o'); hold on;
    plot(ustRe_fit,Q_fit);
    title(['d = ',d50_print],'FontSize',16);
    xlabel('u_{*,Re} (m/s)','FontSize',16);
    ylabel('Q (g/m/s)','FontSize',16);
    ylim([0 35]);
    set(gca,'FontSize',16);
    print([folder_Plots,'FluxUst_Wenglor_d',d50_print,'.png'],'-dpng');
end

%% get fit of u*thr versus shear velocity, save this value
d50_range = [min(d50_threshold_bins), max(d50_threshold_bins)];
C_ustThr_sqrtd50 = mean(ustThr_bins./sqrt(d50_threshold_bins));
C_ustThr_sqrtd50_units = 'm/s/sqrt(mm)';
savepath = strcat(folder_ProcessedData,'C_ustThr_sqrtd50');
save(savepath,'C_ustThr_sqrtd50','C_ustThr_sqrtd50_units');

%% plot u*thr versus shear velocity
figure(N_bins*2+1); clf; hold on;
plot(d50_threshold_bins,ustThr_bins,'o','MarkerSize',15);
plot(d50_range,C_ustThr_sqrtd50*sqrt(d50_range),'b-','LineWidth',3);
xlabel('d_{50} (mm)');
ylabel('u_{*,thr} (m/s)');
set(gca,'FontSize',16);
h_legend = legend({'data','sqrt fit'},'Location','NorthWest');
set(h_legend,'FontSize',16);
print([folder_Plots,'UstThr_D50_bins.png'],'-dpng');