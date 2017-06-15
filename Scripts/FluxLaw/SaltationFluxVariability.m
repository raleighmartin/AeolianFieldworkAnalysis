%% SCRIPT TO INTERPRET SOURCES OF SALTATION FLUX VARIABILITY RELATED TO TURBULENCE

%% initialize
clearvars;

%% Site-specific parameter values
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
z0_Site = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m)
kappa = 0.4; %von Karman parameter
tauit_intercept = [0.137, 0.106, 0.086];
sigma_tauit_intercept = [0.015, 0.020, 0.008];

%% info about times and time scales for sample plot
Site_sampleplot = 'RanchoGuadalupe';
starttime_sampleplot = datetime(2015,3,24,13,31,0);
%endtime_sampleplot = datetime(2015,3,24,13,41,0);
endtime_sampleplot = datetime(2015,3,24,13,36,0);
T_sampleplot = seconds(endtime_sampleplot-starttime_sampleplot); %duration of sample plot
Deltat_sampleplot = duration(0,1,0); %measurement interval for analysis
deltat_sampleplot = duration(0,0,1); %sampling interval for analysis

%% load data 
folder_LoadWindowData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
LoadWindowData_Path = strcat(folder_LoadWindowData,'DataWindows_30min'); %path for 30 minute data
load(LoadWindowData_Path);
folder_LoadThresholdData = '../../AnalysisData/Thresholds/'; %folder for retrieving processed data
LoadThresholdData1_Path = strcat(folder_LoadThresholdData,'ThresholdAnalysisWindows'); %path for loading threshold time window data
LoadThresholdData2_Path = strcat(folder_LoadThresholdData,'ThresholdAnalysisData'); %path for loading information from threshold analysis
load(LoadThresholdData1_Path);
load(LoadThresholdData2_Path);

%% info for plotting
folder_Plots = '../../PlotOutput/Variability/'; %folder for plots
plotname_etasample = ['eta_sample_',Site_sampleplot,'_',datestr(starttime_sampleplot,30),'_',int2str(seconds(Deltat_sampleplot)),'_',num2str(seconds(deltat_sampleplot)),'.png'];
plotpath_etasample = [folder_Plots,plotname_etasample];
plotname_windonly = ['wind_sample_',Site_sampleplot,'_',datestr(starttime_sampleplot,30),'_',int2str(seconds(Deltat_sampleplot)),'_',num2str(seconds(deltat_sampleplot)),'.png'];
plotpath_windonly = [folder_Plots,plotname_windonly];
plotname_windonly_singlepanel = ['wind_sample_singlepanel_',Site_sampleplot,'_',datestr(starttime_sampleplot,30),'_',int2str(seconds(Deltat_sampleplot)),'_',num2str(seconds(deltat_sampleplot)),'.png'];
plotpath_windonly_singlepanel = [folder_Plots,plotname_windonly_singlepanel];
plotfont = 14;

%%%%%%%%%
% PLOTS %
%%%%%%%%%

%% GET INFO FOR PLOTS
%get indices for values
ind_Site = find(strcmp(Site_sampleplot,Sites));
ind_window = intersect(find(StartTime_window{ind_Site}<=starttime_sampleplot),...
    find(EndTime_window{ind_Site}>=endtime_sampleplot));
ind_window_t = intersect(find(t_wind_window{ind_Site}{ind_window}>=starttime_sampleplot),...
    find(t_wind_window{ind_Site}{ind_window}<=endtime_sampleplot));
ind_Deltat = find(Deltat_all == Deltat_sampleplot);
ind_deltat = find(deltat_all == deltat_sampleplot);
ind_subwindow = intersect(find(starttime_all{ind_Site}{ind_Deltat,ind_deltat}>=starttime_sampleplot),...
    find(endtime_all{ind_Site}{ind_Deltat,ind_deltat}<=endtime_sampleplot));
N_ind_subwindow = length(ind_subwindow);

%get values
rho = rho_Site(ind_Site);
z0 = z0_Site(ind_Site);
zU = zU_window{ind_Site}(ind_window);
tauft = tauft_all(ind_Site);
tauit = tauit_all(ind_Site);

%compute fluid and impact threshold wind speed
uft = (sqrt(tauft/rho)/kappa)*log(zU/z0);
uit = (sqrt(tauit/rho)/kappa)*log(zU/z0);

%initialize stress partition values
tau_bar_total = zeros(N_ind_subwindow,1);
tau_bar_inactive = zeros(N_ind_subwindow,1);
tau0_bar_inactive = zeros(N_ind_subwindow,1);
eta_inactive = zeros(N_ind_subwindow,1);
tau_bar_active = zeros(N_ind_subwindow,1);
eta_active = zeros(N_ind_subwindow,1);
eta_deltat = cell(N_ind_subwindow,1);
eta_combined = zeros(N_ind_subwindow,1);

for i=1:N_ind_subwindow
    
    %get fQ and u for subwindow
    fQ = fQ_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)); 
    u = u_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)};

    %initialize eta_deltat    
    eta_deltat{i} = zeros(size(u));
    
    %get uth for subwindow
    if fQ == 0; uth = uft;
    elseif fQ == 1; uth = uit;
    else; uth = uth_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i));
    end
    
    %total computation
    u2_bar = mean(u.^2);
    tau_bar_total(i) = (rho*kappa^2*u2_bar)/log(zU/z0)^2;
    
    %computations for inactive indices
    ind_inactive = find(u<=uth);
    u2_bar_inactive = mean(u(ind_inactive).^2);
    tau_bar_inactive(i) = (rho*kappa^2*u2_bar_inactive)/log(zU/z0)^2;
    eta_deltat{i}(ind_inactive) = 0;
    
    %computations for active indices
    ind_active = find(u>uth);
    u2_bar_active = mean(u(ind_active).^2);
    tau_bar_active(i) = (rho*kappa^2*u2_bar_active)/log(zU/z0)^2;
    eta_deltat{i}(ind_active) = 1-(tauit/(rho*kappa^2))*(log(zU/z0)^2)./(u(ind_active).^2);
    eta_deltat{i}(eta_deltat{i}<0) = 0; %change negative points to 0
     
    %compute eta values
    if fQ == 0
        tau0_bar_inactive(i) = 0;
        eta_active(i) = 0;
    elseif fQ == 1
        tau0_bar_inactive(i) = tauit;
        eta_active(i) = 1-(tauit/tau_bar_active(i));
    else
        tau0_bar_inactive(i) = fQ*tauit+(1-fQ)*tau_bar_inactive(i);
        eta_active(i) = fQ*(1-(tauit/tau_bar_active(i)));
    end
    eta_inactive(i) = 1 - tau0_bar_inactive(i)/tau_bar_total(i);
    
    %computations for combined
    eta_combined(i) = mean(eta_deltat{i});
end

%% plot sample timeseries showing methods
figure(1); clf; hold on;

%wind
subplot(5,1,1); hold on;
plot(t_wind_window{ind_Site}{ind_window}(ind_window_t),u_window{ind_Site}{ind_window}(ind_window_t),'c'); %plot raw u
for i = 1:N_ind_subwindow
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        u_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'b'); %plot averaged u
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        uth_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))*[1 1],'r','LineWidth',3); %plot threshold u
    if i==1
        plot([starttime_sampleplot endtime_sampleplot],[uit, uit],'k--','LineWidth',2); %plot impact threshold
        plot([starttime_sampleplot endtime_sampleplot],[uft, uft],'k--','LineWidth',2); %plot fluid threshold
    end
end
xlim([starttime_sampleplot endtime_sampleplot])
ylim([floor(min(u_window{ind_Site}{ind_window}(ind_window_t))),ceil(max(u_window{ind_Site}{ind_window}(ind_window_t)))]);
legend('raw','{\delta}t avg','{\Delta}t u_{th}','it, ft','Location','EastOutside');
ylabel('wind speed, $$u$$ (m/s)','interpreter','latex');
set(gca,'FontSize',plotfont);

%particle counts
subplot(5,1,2); hold on;
for i = 1:length(ind_subwindow)
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        n_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'r'); %plot averaged counts
end
xlim([starttime_sampleplot endtime_sampleplot])
ylabel('particle rate, $$n$$ (1/s)','interpreter','latex');
legend('{\delta}t avg','Location','EastOutside');
set(gca,'FontSize',plotfont);

%saltation activity
subplot(5,1,3); hold on;
for i = 1:N_ind_subwindow
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        fQ_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))*[1 1],'r','LineWidth',3); %plot activity
end
xlim([starttime_sampleplot endtime_sampleplot])
ylim([0 1]);
ylabel('saltation activity, $$f_Q$$','interpreter','latex');
legend('{\Delta}t f_Q','Location','EastOutside');
set(gca,'FontSize',plotfont);

%stresses
subplot(5,1,4); hold on;
for i = 1:N_ind_subwindow
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        tau_bar_total(i)*[1 1],'k','LineWidth',3); %plot tau_total
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        tau_bar_inactive(i)*[1 1],'r','LineWidth',3); %plot tau_inactive
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        tau_bar_active(i)*[1 1],'b','LineWidth',3); %plot tau_inactive
    if i==1
        plot([starttime_sampleplot endtime_sampleplot],[tauit, tauit],'k--','LineWidth',1); %plot impact threshold
        plot([starttime_sampleplot endtime_sampleplot],[tauft, tauft],'k--','LineWidth',1); %plot fluid threshold
    end
end
xlim([starttime_sampleplot endtime_sampleplot])
ylabel('stress, $$\tau$$ (Pa)','interpreter','latex');
legend('\tau_{total}','\tau_{inactive}','\tau_{active}','it, ft','Location','EastOutside');
set(gca,'FontSize',plotfont);

%stress partition
subplot(5,1,5); hold on;
for i = 1:N_ind_subwindow
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        eta_deltat{i},'k'); %plot eta_deltat
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        eta_inactive(i)*[1 1],'r','LineWidth',2); %plot eta for inactive method
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        eta_active(i)*[1 1],'b','LineWidth',2); %plot eta for active method
    plot([starttime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i)),endtime_all{ind_Site}{ind_Deltat,ind_deltat}(ind_subwindow(i))],...
        eta_combined(i)*[1 1],'k','LineWidth',2); %plot combined eta
end
xlim([starttime_sampleplot endtime_sampleplot])
ylabel('stress partition, $$\eta$$','interpreter','latex');
legend('\eta_{{\delta}t}','\eta_{inactive}','\eta_{active}','<\eta_{{\delta}t}>','Location','EastOutside');
set(gca,'FontSize',plotfont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 12]);
print(plotpath_etasample,'-dpng');

%% plot sample timeseries with vertical wind
figure(2); clf; hold on;

%horizontal wind
subplot(3,1,1); hold on;
plot(t_wind_window{ind_Site}{ind_window}(ind_window_t),u_window{ind_Site}{ind_window}(ind_window_t),'c'); %plot raw u
for i = 1:N_ind_subwindow
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        u_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'b'); %plot averaged u
    if i==1
        plot([starttime_sampleplot endtime_sampleplot],[uit, uit],'k--','LineWidth',1); %plot impact threshold
        plot([starttime_sampleplot endtime_sampleplot],[uft, uft],'k--','LineWidth',1); %plot fluid threshold
    end
end
xlim([starttime_sampleplot endtime_sampleplot])
ylim([floor(min(u_window{ind_Site}{ind_window}(ind_window_t))),ceil(max(u_window{ind_Site}{ind_window}(ind_window_t)))]);
ylabel('horizontal wind speed, $$u$$ (m/s)','interpreter','latex');
grid on;
set(gca,'FontSize',plotfont);

%vertical wind
subplot(3,1,2); hold on;
plot(t_wind_window{ind_Site}{ind_window}(ind_window_t),w_window{ind_Site}{ind_window}(ind_window_t),'c'); %plot raw w
for i = 1:N_ind_subwindow
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        w_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'b'); %plot averaged w
end
xlim([starttime_sampleplot endtime_sampleplot])
ylim([floor(min(w_window{ind_Site}{ind_window}(ind_window_t))),ceil(max(w_window{ind_Site}{ind_window}(ind_window_t)))]);
ylim([-1 1])
ylabel('vertical wind speed, $$w$$ (m/s)','interpreter','latex');
grid on;
set(gca,'FontSize',plotfont);

%particle counts
subplot(3,1,3); hold on;
for i = 1:length(ind_subwindow)
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        n_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'r'); %plot averaged counts
end
xlim([starttime_sampleplot endtime_sampleplot])
ylabel('particle rate, $$n$$ (1/s)','interpreter','latex');
grid on;
set(gca,'FontSize',plotfont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 12]);
print(plotpath_windonly,'-dpng');

%% plot sample timeseries with vertical wind - single panel
figure(3); clf; hold on;

%horizontal wind
for i = 1:N_ind_subwindow
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        40*(u_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)}-4),'b'); %plot averaged u
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        200*w_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'r'); %plot averaged w
    plot(t_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},...
        n_all{ind_Site}{ind_Deltat,ind_deltat}{ind_subwindow(i)},'k','LineWidth',1); %plot averaged counts
end
xlim([starttime_sampleplot endtime_sampleplot])
ylabel('particle rate, $$n$$ (1/s)','interpreter','latex');
legend('u','w','n','Location','NorthWest');
grid on;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca,'FontSize',plotfont);

%print plot
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print(plotpath_windonly_singlepanel,'-dpng');