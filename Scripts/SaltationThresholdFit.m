%% SEPARATE AND PLOT RESULTS BY GRAIN SIZE BINS, ANALYZE SALTATION THRESHOLD FROM FLUX-STRESS RELATIONSHIP

%initialize
clearvars;

%% INPUT INFORMATION

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for storing data output
folder_AnalysisData = '../AnalysisData/'; %folder for analysis data
folder_Plots = '../PlotOutput/SaltationThreshold/'; %folder for plots

%load stress flux window data
load(strcat(folder_AnalysisData,'StressFluxWindows_all'));

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

% %% GO THROUGH EACH SITE AND EACH DAY
% for i = 1:N_Sites
%     Dates = unique(date_all{i});
%     N_Dates = length(Dates);
%     for j = 1:N_Dates
%         %get all stress/flux values for date
%         ind_date = find(date_all{i}==Dates(j));
%         tau_date = tauRe_cal_all{i}(ind_date);
%         Q_date = Q_all{i}(ind_date);
%         
%         %get subset for this date with detected transport
%         ind_trans = find(Q_date>0);
%         tau_trans = tau_date(ind_trans);
%         Q_trans = Q_date(ind_trans);
%         
%         %fit line to this
%         P = polyfit(tau_trans,Q_trans,1);
%         tau_thr = -P(2)/P(1); %(Pa)
%         tau_fit = [tau_thr, max(tau_trans)];
%         Q_fit = [0, P(1)*(max(tau_trans)-tau_thr)];
%         
%         figure(1); clf;
%         plot(tau_date,Q_date,'o'); hold on;
%         plot(tau_fit,Q_fit,'k');
%         xlabel('\tau (Pa)');
%         ylabel('Q (g/m/s)');
%         legend('data',['\tau_{th} = ',num2str(tau_thr),' Pa'],'Location','NorthWest');
%         title([Sites{i},': ',datestr(Dates(j))]);
%         set(gca,'FontSize',16);
%         print([folder_Plots,'FluxTau_',Sites{i},'_',datestr(Dates(j)),'.png'],'-dpng');
%     end
% end

%% GO THROUGH EACH GROUP OF DAYS

datemin = {datetime(2014,11,13); %Jeri
    datetime(2015,3,23); %Rancho
    %[datetime(2015,5,15),datetime(2015,5,23),datetime(2015,6,1)]}; %Oceano
    [datetime(2015,5,15),datetime(2015,5,16),datetime(2015,5,23),datetime(2015,6,1)]}; %Oceano
datemax = {datetime(2014,11,20); %Jeri
    datetime(2015,3,24); %Rancho
    %[datetime(2015,5,19),datetime(2015,5,31),datetime(2015,6,4)]}; %Oceano
    [datetime(2015,5,15),datetime(2015,5,19),datetime(2015,5,31),datetime(2015,6,4)]}; %Oceano

for i = 1:N_Sites
   N_dategroups = length(datemin{i}); %number of date groups for site
   for j = 1:N_dategroups %go through each date group
        %get fluxes / stresses in date group
        date_ind = intersect(find(date_all{i}>=datemin{i}(j)),find(date_all{i}<=datemax{i}(j)));
        tau_group = tauRe_cal_all{i}(date_ind);
        Q_group = Q_all{i}(date_ind);
        
        %get subset for this date group with detected transport
        Q_thres = 1:10;
        tau_thres = zeros(size(Q_thres))*NaN;
        for k = 1:length(Q_thres)
            tau_Q = tau_group(round(Q_group)==Q_thres(k));
            tau_thres(k) = mean(tau_Q(~isnan(tau_Q)));
        end
        Q_thres = Q_thres(~isnan(tau_thres));
        tau_thres = tau_thres(~isnan(tau_thres));
        
%         %get subset for this date group with detected transport
%         %ind_thres = find(Q_group>1);
%         ind_thres = intersect(find(Q_group>=2),find(Q_group<=10));
%         tau_thres = tau_group(ind_thres);
%         Q_thres = Q_group(ind_thres);
        
        %fit line to this
        P = polyfit(tau_thres,Q_thres,1);
        tau_thr = -P(2)/P(1); %(Pa)
        tau_fit = [tau_thr, max(tau_thres)];
        Q_fit = [0, P(1)*(max(tau_thres)-tau_thr)];
        
        %plot this
        figure(1); clf; hold on;
        plot(tau_group,Q_group,'o'); 
        plot(tau_fit,Q_fit,'k');
        xlabel('\tau (Pa)');
        ylabel('Q (g/m/s)');
        legend('data',['\tau_{th} = ',num2str(tau_thr),' Pa'],'Location','NorthWest');
        title([Sites{i},': ',datestr(datemin{i}(j)),' - ',datestr(datemax{i}(j))]);
        set(gca,'FontSize',16);
        print([folder_Plots,'FluxTau_',Sites{i},'_',datestr(datemin{i}(j)),'_',datestr(datemax{i}(j)),'.png'],'-dpng');
    end
end

%PLOT ALL DATA TOGETHER VERSUS UST_EX AND TAU_EX
%Apply thresholds to compute ustex and tauex based on thresholds by date
ustex_all = cell(N_Sites,1);
tauex_all = cell(N_Sites,1);
for i=1:N_Sites
    ustex_all{i} = zeros(size(date_all{i}));
    tauex_all{i} = zeros(size(date_all{i}));
    %apply empirical threshold - Jericoacoara - all dates
    if i == 1;
        ustex_all{i} = ustRe_cal_all{i}-0.3474;
        tauex_all{i} = tauRe_cal_all{i}-0.14841;
    %apply empirical threshold - Rancho Guadalupe - all dates
    elseif i == 2;
       ustex_all{i} = ustRe_cal_all{i}-0.2839;
       tauex_all{i} = tauRe_cal_all{i}-0.099111;
    elseif i == 3;
        %apply empirical threshold - Oceano - May 15-19
        date_ind = intersect(find(date_all{i}>=datetime(2015,5,15)),find(date_all{i}<=datetime(2015,5,19)));
        ustex_all{i}(date_ind) = ustRe_cal_all{i}(date_ind)-0.2794;
        tauex_all{i}(date_ind) = tauRe_cal_all{i}(date_ind)-0.096036;
        %apply empirical threshold - Oceano - May 23-31
        date_ind = intersect(find(date_all{i}>=datetime(2015,5,23)),find(date_all{i}<=datetime(2015,5,31)));
        ustex_all{i}(date_ind) = ustRe_cal_all{i}(date_ind)-0.2903;
        tauex_all{i}(date_ind) = tauRe_cal_all{i}(date_ind)-0.10364;
        %apply empirical threshold - Oceano - June 1-4
        date_ind = intersect(find(date_all{i}>=datetime(2015,6,1)),find(date_all{i}<=datetime(2015,6,4)));
        ustex_all{i}(date_ind) = ustRe_cal_all{i}(date_ind)-0.2792;
        tauex_all{i}(date_ind) = tauRe_cal_all{i}(date_ind)-0.095865;
    end
end

% figure(2); clf; hold on;
% for i = 1:N_Sites
% %for i = 1:2
%     plot(ustex_all{i},Q_all{i},Markers{i});
%     xlabel('u_{*,ex} (m/s)');
%     ylabel('Q (g/m/s)');
% end
% 
% figure(3); clf; hold on;
% for i = 1:N_Sites
% %for i = 1:2
%     plot(tauex_all{i},Q_all{i},Markers{i});
%     xlabel('\tau_{ex} (Pa)');
%     ylabel('Q (g/m/s)');
% end
% 
% %% CALCULATIONS
% %set bin edges
% d50_bins_min = [0.33, 0.39, 0.51];
% d50_bins_max = [0.36, 0.44, 0.54];
% N_bins = length(d50_bins_min);
% 
% %initialize list of calculated values for bins
% tauThr_bins = zeros(N_bins,1);
% ustThr_bins = zeros(N_bins,1);
% d50_threshold_bins = zeros(N_bins,1);
% 
% %create cell array of thresholds for each site
% ustThr_Site = cell(N_Sites,1);
% tauThr_Site = cell(N_Sites,1);
% 
% %go through each bin
% for i = 1:N_bins
% 
%     %initialize lists for bins
%     tauRe_bin = [];
%     ustRe_bin = [];
%     Q_bin = [];
%     d50_bin = [];
%     
%     %go through each site
%     for j = 1:N_Sites
%                 
%         %get values for this bin for this site
%         ind_bin_Site = find((d50_all{j}>=d50_bins_min(i))&(d50_all{j}<=d50_bins_max(i)));
%         ind_bin_Site = intersect(ind_bin_Site,find(Q_all{j}>0)); %include only indices with positive Q
%         
%         %if bin contains values for this Site, get values and add to overall list
%         if ~isempty(ind_bin_Site)
%             tauRe_bin_Site = tauRe_cal_all{j}(ind_bin_Site);
%             ustRe_bin_Site = ustRe_cal_all{j}(ind_bin_Site);
%             Q_bin_Site = Q_all{j}(ind_bin_Site);
%             d50_bin_Site = d50_all{j}(ind_bin_Site);
%             
%             %add to lists for bins
%             tauRe_bin = [tauRe_bin; tauRe_bin_Site];
%             ustRe_bin = [ustRe_bin; ustRe_bin_Site];
%             Q_bin = [Q_bin; Q_bin_Site];
%             d50_bin = [d50_bin; d50_bin_Site];
%         end
%     end
%     
%     %for each bin, fit to determine threshold
%     P_linear = polyfit(tauRe_bin,Q_bin,1);
%     tauThr = -P_linear(2)/P_linear(1);
%     ustThr = sqrt(tauThr/rho_a);
%     
%     %add to lists for each bin
%     tauThr_bins(i) = tauThr;
%     ustThr_bins(i) = ustThr;
%     d50_threshold_bins(i) = mean(d50_bin);
%     d50_print = [int2str(round(1000*mean(d50_bin))),'um'];
%     
%     %get values for plotting fit
%     tauRe_fit = linspace(tauThr,max(tauRe_bin),50);
%     ustRe_fit = sqrt(tauRe_fit/rho_a);
%     Q_fit = polyval(P_linear,tauRe_fit);
%     
%     %plot flux versus shear stress
%     figure(i*2-1+1); clf;
%     plot(tauRe_bin,Q_bin, 'o'); hold on;
%     plot(tauRe_fit,Q_fit);
%     title(['d = ',d50_print],'FontSize',16);
%     xlabel('\tau_{Re} (Pa)','FontSize',16);
%     ylabel('Q (g/m/s)','FontSize',16);
%     ylim([0 35]);
%     set(gca,'FontSize',16);
%     print([folder_Plots,'FluxTau_Wenglor_d',d50_print,'.png'],'-dpng');
% 
%     %plot flux versus shear velocity
%     figure(i*2+1); clf;
%     plot(ustRe_bin,Q_bin, 'o'); hold on;
%     plot(ustRe_fit,Q_fit);
%     title(['d = ',d50_print],'FontSize',16);
%     xlabel('u_{*,Re} (m/s)','FontSize',16);
%     ylabel('Q (g/m/s)','FontSize',16);
%     ylim([0 35]);
%     set(gca,'FontSize',16);
%     print([folder_Plots,'FluxUst_Wenglor_d',d50_print,'.png'],'-dpng');
% end
% 
% %% get fit of u*thr versus shear velocity, save this value
% d50_range = [min(d50_threshold_bins), max(d50_threshold_bins)];
% C_ustThr_sqrtd50 = mean(ustThr_bins./sqrt(d50_threshold_bins));
% C_ustThr_sqrtd50_units = 'm/s/sqrt(mm)';
% savepath = strcat(folder_AnalysisData,'C_ustThr_sqrtd50');
% save(savepath,'C_ustThr_sqrtd50','C_ustThr_sqrtd50_units');
% 
% %% plot u*thr versus shear velocity
% figure(N_bins*2+1); clf; hold on;
% plot(d50_threshold_bins,ustThr_bins,'o','MarkerSize',15);
% plot(d50_range,C_ustThr_sqrtd50*sqrt(d50_range),'b-','LineWidth',3);
% xlabel('d_{50} (mm)');
% ylabel('u_{*,thr} (m/s)');
% set(gca,'FontSize',16);
% h_legend = legend({'data','sqrt fit'},'Location','NorthWest');
% set(h_legend,'FontSize',16);
% print([folder_Plots,'UstThr_D50_bins.png'],'-dpng');