%% SCRIPT TO COMPARE ACTIVITY PREDICTION TO OBSERVATION

%% initialize
clearvars;
close all;

%% Site-specific parameter values
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
kappa = 0.4; %von Karman parameter

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,1); %sampling interval for analysis

%% plotting info
PlotFont = 14;
PlotMarkers = {'s','d','o'};
PlotColors = {'r','m','b'};

%% binning info
fQ_min = 0.05;
fQ_max = 0.95;
fQ_bin_minrange = 0.1;
fQ_bin_maxrange = 0.2;
fQ_bin_N_min = 3;

%% folders for loading and saving data
folder_LoadSubwindowData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_LoadSubwindowCalcs = '../../AnalysisData/Windowing/'; %folder for retrieving subwindow calcs
folder_LoadThresholdData = '../../AnalysisData/Thresholds/'; %folder for loading threshold values
folder_LoadRoughnessData = '../../AnalysisData/Roughness/'; %folder for loading roughness data
folder_Functions = '../Functions/'; %folder with functions
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% paths for loading and saving data
LoadSubwindowData_Path = strcat(folder_LoadSubwindowData,'DataSubwindows_30min_Unrestricted'); %path for loading time window data
LoadSubwindowCalcs_Path = strcat(folder_LoadSubwindowCalcs,'DataSubwindowCalcs_30min_Unrestricted'); %path for loading time window data
LoadThresholdData_Path = strcat(folder_LoadThresholdData,'ThresholdAnalysisData'); %path for threshold data
LoadRoughnessData_Path = strcat(folder_LoadRoughnessData,'RoughnessCalcs_30min_Restricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'ActivityPredictionData'); %path for saving output data

%% load data and rename variables as necessary
load(LoadSubwindowData_Path); %load primary data
load(LoadSubwindowCalcs_Path); %load 30-minute values
load(LoadThresholdData_Path); %load 30-minute values
load(LoadRoughnessData_Path); %load Roughness values
addpath(folder_Functions); %point MATLAB to location of functions

%% information about measurement and sampling interval indices for analysis
ind_Deltat = find(Deltat_all == Deltat_analysis); %index for measurement interval
ind_deltat = find(deltat_all == deltat_analysis); %index for sampling interval

%% get z0 for Site
z0_Site = z0Re_Q_fit;

%% initialize unbinned values
fplus_all = cell(N_Sites,1); %fraction of time above fluid threshold
fint_all = cell(N_Sites,1); %fraction of time between impact and fluid threshold
fint_down_all = cell(N_Sites,1); %fraction of time between impact and fluid threshold from downward crossing
ind_fplus_all = cell(N_Sites,1); %indices for times above fluid threshold
ind_fint_all = cell(N_Sites,1); %indices for times between impact and fluid threshold
ind_fint_down_all = cell(N_Sites,1); %indices for times between impact and fluid threshold from downward crossing
fQpred_ft_all = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold
fQpred_it_all = cell(N_Sites,1); %predicted transport fraction based on time above impact threshold - average
fQpred_hyst_all = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold and hysteresis - average

%% initialize binned values
fQ_bin_avg_all = cell(N_Sites,1); %observed transport activity - average
fQ_bin_SE_all = cell(N_Sites,1); %observed transport activity - standard error
fQpred_ft_bin_avg_all = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold - average
fQpred_ft_bin_SE_all = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold - standard error
fQpred_it_bin_avg_all = cell(N_Sites,1); %predicted transport fraction based on time above impact threshold - average
fQpred_it_bin_SE_all = cell(N_Sites,1); %predicted transport fraction based on time above impact threshold - standard error
fQpred_hyst_bin_avg_all = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold and hysteresis - average
fQpred_hyst_bin_SE_all = cell(N_Sites,1); %predicted transport fraction based on time above fluid threshold and hysteresis - standard error

%% initialize normalized chi2 values
chi2nu_ft_all = zeros(N_Sites,1); %chi2nu for fluid threshold control
chi2nu_it_all = zeros(N_Sites,1); %chi2nu for impact threshold control
chi2nu_hyst_all = zeros(N_Sites,1); %chi2nu for hysteresis control

%% go through sites
for i = 1:N_Sites
%for i = 3    
    %get values for analysis
    StartTime = StartTime_subwindow{i}{ind_Deltat,ind_deltat};
    EndTime = EndTime_subwindow{i}{ind_Deltat,ind_deltat};
    fQ = fQ_subwindow_all{i}{ind_Deltat,ind_deltat};
    u = u_subwindow{i}{ind_Deltat,ind_deltat};
    t = t_subwindow{i}{ind_Deltat,ind_deltat};
    zU = zU_subwindow{i}{ind_Deltat,ind_deltat};
    ustit = sqrt(tauit_all(i)/rho_Site(i));
    uit = (ustit/kappa).*log(zU./z0_Site(i));
    ustft = sqrt(tauft_all(i)/rho_Site(i));
    uft = (ustft/kappa).*log(zU./z0_Site(i));
    N_subwindows = length(StartTime);
    
    %initialize unbinned values
    fplus_all{i} = zeros(N_subwindows,1);
    fint_all{i} = zeros(N_subwindows,1);
    fint_down_all{i} = zeros(N_subwindows,1);
    ind_fplus_all{i} = cell(N_subwindows,1); %indices for times above fluid threshold
    ind_fint_all{i} = cell(N_subwindows,1); %indices for times between impact and fluid threshold
    ind_fint_down_all{i} = cell(N_subwindows,1); %indices for times between impact and fluid threshold from downward crossing
    fQpred_ft_all{i} = zeros(N_subwindows,1); %predicted transport fraction based on time above fluid threshold
    fQpred_it_all{i} = zeros(N_subwindows,1); %predicted transport fraction based on time above impact threshold
    fQpred_hyst_all{i} = zeros(N_subwindows,1); %predicted transport fraction based on time above fluid threshold and hysteresis
    
    %initialize list for final state of wind in each subwindow
    final_state_all = zeros(N_subwindows,1); 
    % -1 if most recent time was predicted non-transport
    % +1 if most recent time was 
    % 0 if most recent time unknown
    
    for j = 1:N_subwindows
        %get final state from previous window if it lines up with current window, use this to set init state
        ind_last = find(EndTime==StartTime(j));
        if ~isempty(ind_last)
            init_state = final_state_all(ind_last);
        else
            init_state = 0;
        end
                
        %determine fractions of u in different regions - window-averaged u
        [fplus,~,fint,~,fint_down,ind_fplus,~,ind_fint,~,ind_fint_down,final_state] = ...
            CalculateWindHysteresis(u{j},uft(j),uit(j),init_state);

        %record last state
        final_state_all(j) = final_state;
        
        %add these fractions to list of values from window-averaged winds
        fplus_all{i}(j) = fplus; %fraction of time with u_avg above fluid threshold
        fint_all{i}(j) = fint; %fraction of time with u_avg in intermediate zone
        fint_down_all{i}(j) = fint_down; %fraction of time with u_avg in intermediate zone with downward crossing
        
        %add indices to list of values
        ind_fplus_all{i}{j} = ind_fplus; %indices for times above fluid threshold
        ind_fint_all{i}{j} = ind_fint; %indices for times between impact and fluid threshold
        ind_fint_down_all{i}{j} = ind_fint_down; %indices for times between impact and fluid threshold from downward crossing
        
        %compute predicted values
        fQpred_ft_all{i}(j) = fplus; %predicted transport fraction based on time above fluid threshold
        fQpred_it_all{i}(j) = fplus + fint; %predicted transport fraction based on time above impact threshold
        fQpred_hyst_all{i}(j) = fplus + fint_down; %predicted transport fraction based on time above fluid threshold and hysteresis      
    end
    
    %% perform binning
    
    %get indices for binning
    ind_binning = intersect(find(fQ>=fQ_min),find(fQ<=fQ_max)); 

    %keep only values for binning
    fQ_binning = fQ(ind_binning);
    fQpred_ft_binning = fQpred_ft_all{i}(ind_binning);
    fQpred_it_binning = fQpred_it_all{i}(ind_binning);
    fQpred_hyst_binning = fQpred_hyst_all{i}(ind_binning);
    
    %get binned values for fQ
    [~, ~, fQ_bin_min, fQ_bin_max, fQ_bin_avg, fQ_bin_SE] = Binning(fQ_binning, fQ_bin_minrange, fQ_bin_maxrange, fQ_bin_N_min);
    N_fQ_bins = length(fQ_bin_min);
    fQ_bin_avg_all{i} = fQ_bin_avg;
    fQ_bin_SE_all{i} = fQ_bin_SE;
    
    % initialize binned values
    fQpred_ft_bin_avg_all{i} = zeros(N_fQ_bins,1); %predicted transport fraction based on time above fluid threshold - average
    fQpred_ft_bin_SE_all{i} = zeros(N_fQ_bins,1); %predicted transport fraction based on time above fluid threshold - standard error
    fQpred_it_bin_avg_all{i} = zeros(N_fQ_bins,1); %predicted transport fraction based on time above impact threshold - average
    fQpred_it_bin_SE_all{i} = zeros(N_fQ_bins,1); %predicted transport fraction based on time above impact threshold - standard error
    fQpred_hyst_bin_avg_all{i} = zeros(N_fQ_bins,1); %predicted transport fraction based on time above fluid threshold and hysteresis - average
    fQpred_hyst_bin_SE_all{i} = zeros(N_fQ_bins,1); %predicted transport fraction based on time above fluid threshold and hysteresis - standard error
    
    %get binned values for fQpred
    for k=1:N_fQ_bins   
        fQ_bin_ind = find(fQ_binning>=fQ_bin_min(k)&fQ_binning<=fQ_bin_max(k));
        fQpred_ft_bin_avg_all{i}(k) = mean(fQpred_ft_binning(fQ_bin_ind));
        fQpred_ft_bin_SE_all{i}(k) = std(fQpred_ft_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
        fQpred_it_bin_avg_all{i}(k) = mean(fQpred_it_binning(fQ_bin_ind));
        fQpred_it_bin_SE_all{i}(k) = std(fQpred_it_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
        fQpred_hyst_bin_avg_all{i}(k) = mean(fQpred_hyst_binning(fQ_bin_ind));
        fQpred_hyst_bin_SE_all{i}(k) = std(fQpred_hyst_binning(fQ_bin_ind))/sqrt(length(fQ_bin_ind));
    end
    
    %compute chi2
    chi2nu_ft_all(i) = sum((fQ_bin_avg_all{i}-fQpred_ft_bin_avg_all{i}).^2./fQpred_ft_bin_SE_all{i}.^2)./length(fQ_bin_avg_all{i}); %chi2nu for fluid threshold control
    chi2nu_it_all(i) = sum((fQ_bin_avg_all{i}-fQpred_it_bin_avg_all{i}).^2./fQpred_it_bin_SE_all{i}.^2)./length(fQ_bin_avg_all{i}); %chi2nu for impact threshold control
    chi2nu_hyst_all(i) = sum((fQ_bin_avg_all{i}-fQpred_hyst_bin_avg_all{i}).^2./fQpred_hyst_bin_SE_all{i}.^2)./length(fQ_bin_avg_all{i}); %chi2nu for hysteresis control
end

%display chi2 values
round(chi2nu_ft_all,1)
round(chi2nu_it_all,1)
round(chi2nu_hyst_all,1)

%%%%%%%%%%%%%
% SAVE DATA %
%%%%%%%%%%%%%
clear *binning
save(SaveData_Path,'fint_all','fint_down_all','fplus_all',...
    'fQ_bin_avg_all','fQ_bin_SE_all','fQpred*',...
    'ind_fint_all','ind_fint_down_all','ind_fplus_all',...
    'rho_Site','z0_Site','kappa');

%%%%%%%%%
% PLOTS %
%%%%%%%%%

%% plot fQ_pred versus fQ for sites

%initialize plot
figure(1); clf; 
for i = 1:N_Sites
    subplot(1,N_Sites,i); clf;
end

%plot average values
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(fQ_bin_avg_all{i},fQpred_it_bin_avg_all{i},PlotMarkers{1},'Color',PlotColors{1});
    plot(fQ_bin_avg_all{i},fQpred_hyst_bin_avg_all{i},PlotMarkers{2},'Color',PlotColors{2});
    plot(fQ_bin_avg_all{i},fQpred_ft_bin_avg_all{i},PlotMarkers{3},'Color',PlotColors{3});
end

%plot SE values
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],fQpred_it_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{1});
        plot(fQ_bin_avg_all{i}(k)*[1 1],fQpred_it_bin_avg_all{i}(k)+fQpred_it_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{1});
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],fQpred_hyst_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{2});
        plot(fQ_bin_avg_all{i}(k)*[1 1],fQpred_hyst_bin_avg_all{i}(k)+fQpred_hyst_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{2});
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],fQpred_ft_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{3});
        plot(fQ_bin_avg_all{i}(k)*[1 1],fQpred_ft_bin_avg_all{i}(k)+fQpred_ft_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{3});
    end
end

%plot 1-1 lines
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot([0,1],[0,1],'k--');
end

%annotate plot
paneltext = {'(a)';'(b)';'(c)'};
for i = 1:N_Sites
    subplot(1,N_Sites,i);
    text(0.05,0.95,paneltext{i},'FontSize',PlotFont)
    xlabel('observed activity, $$f_Q$$','interpreter','latex');
    if i==1
        ylabel('predicted activities, $$f_{Q,ft}, f_{Q,dual}, f_{Q,it}$$','interpreter','latex');
    end
    title(SiteNames{i});
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
end

%add legend just to last subplot
subplot(1,N_Sites,N_Sites);
legend({'impact thres.','dual thres.','fluid thres.'},'Location','SouthEast');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 5]);
print([folder_Plots,'fQpred_fQ.png'],'-dpng');