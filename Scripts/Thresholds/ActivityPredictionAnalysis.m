%% SCRIPT TO COMPARE ACTIVITY PREDICTION TO OBSERVATION

%% initialize
clearvars;

%% Site-specific parameter values
rho_Site = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
z0_Site = [1e-4, 1e-4, 1e-4]; %aerodynamic roughness length (m)
kappa = 0.4; %von Karman parameter

%% info about time scales for core analysis
Deltat_analysis = duration(0,1,0); %measurement interval for analysis
deltat_analysis = duration(0,0,1); %sampling interval for analysis

%% plotting info
PlotFont = 14;
PlotMarkers = {'s','d','o','<','>'};
PlotColors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};

%% binning info
fQ_min = 0.05;
fQ_max = 0.95;
fQ_bin_minrange = 0.1;
fQ_bin_maxrange = 0.2;
fQ_bin_N_min = 3;

%% load data 
folder_LoadData_1 = '../../AnalysisData/Thresholds/'; %folder for retrieving thresholds data
LoadData_1_Path = strcat(folder_LoadData_1,'ThresholdAnalysisData.mat'); %path for loading flux law windows
load(LoadData_1_Path);
folder_LoadData_2 = '../../AnalysisData/Thresholds/'; %folder for retrieving activity and wind data
LoadData_2_Path = strcat(folder_LoadData_2,'ThresholdAnalysisWindows'); %path for loading activity and wind data
load(LoadData_2_Path);

%% load functions
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% specify where to save data and plots
folder_SaveData = '../../AnalysisData/Thresholds/'; %folder for outputs of this analysis
SaveData_Path = strcat(folder_SaveData,'ActivityPredictionData'); %path for saving output data
folder_Plots = '../../PlotOutput/Thresholds/'; %folder for plots

%% information about measurement and sampling interval indices for analysis
ind_Deltat = find(Deltat_all == Deltat_analysis); %index for measurement interval
ind_deltat = find(deltat_all == deltat_analysis); %index for sampling interval

%% initialize unbinned values
fplus_all = cell(N_Sites,1); %fraction of time above fluid threshold
fint_all = cell(N_Sites,1); %fraction of time between impact and fluid threshold
fint_down_all = cell(N_Sites,1); %fraction of time between impact and fluid threshold from downward crossing
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

%% go through sites
for i = 1:N_Sites
    
    %get values for analysis
    starttime = starttime_all{i}{ind_Deltat,ind_deltat};
    endtime = endtime_all{i}{ind_Deltat,ind_deltat};
    fQ = fQ_all{i}{ind_Deltat,ind_deltat};
    u = u_all{i}{ind_Deltat,ind_deltat};
    zU = zU_all{i}{ind_Deltat,ind_deltat};
    ustit = sqrt(tauit_all(i)/rho_Site(i));
    uit = (ustit/kappa).*log(zU./z0_Site(i));
    ustft = sqrt(tauft_all(i)/rho_Site(i));
    uft = (ustft/kappa).*log(zU./z0_Site(i));
    N_subwindows = length(starttime);
    
    %initialize unbinned values
    fplus_all{i} = zeros(N_subwindows,1);
    fint_all{i} = zeros(N_subwindows,1);
    fint_down_all{i} = zeros(N_subwindows,1);
    fQpred_ft_all{i} = zeros(N_subwindows,1); %predicted transport fraction based on time above fluid threshold
    fQpred_it_all{i} = zeros(N_subwindows,1); %predicted transport fraction based on time above impact threshold
    fQpred_hyst_all{i} = zeros(N_subwindows,1); %predicted transport fraction based on time above fluid threshold and hysteresis
    
    %initialize list for last state of wind
    laststate_all = zeros(N_subwindows,1); 
    % -1 if most recent time was below it
    % +1 if most recent time was above ft
    % 0 if most recent time unknown
    
    for j = 1:N_subwindows
        
        %get last state from previous window if it exists, use this to set init state
        ind_last = find(starttime(j)==endtime);
        if ~isempty(ind_last)
            init_state = laststate_all(ind_last);
        else
            init_state = 0;
        end
                
        %determine fractions of u in different regions - window-averaged u
        [fplus,fminus,fint,fint_up,fint_down] = CalculateWindHysteresis(u{j},uft(j),uit(j),init_state);

        %add these fractions to list of values from window-averaged winds
        fplus_all{i}(j) = fplus; %fraction of time with u_avg above fluid threshold
        fint_all{i}(j) = fint; %fraction of time with u_avg in intermediate zone
        fint_down_all{i}(j) = fint_down; %fraction of time with u_avg in intermediate zone with downward crossing
        
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
end

%% plot fQ_pred versus fQ for sites

%initialize plot
figure(3); clf; 
for i = 1:N_Sites
    subplot(1,N_Sites,i); clf;
end

%plot average values
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(fQ_bin_avg_all{i},fQpred_ft_bin_avg_all{i},PlotMarkers{1},'Color',PlotColors{1});
    plot(fQ_bin_avg_all{i},fQpred_it_bin_avg_all{i},PlotMarkers{2},'Color',PlotColors{2});
    plot(fQ_bin_avg_all{i},fQpred_hyst_bin_avg_all{i},PlotMarkers{3},'Color',PlotColors{3});
end

%plot SE values
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    N_fQ_bins = length(fQ_bin_avg_all{i});
    for k = 1:N_fQ_bins
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],fQpred_ft_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{1});
        plot(fQ_bin_avg_all{i}(k)*[1 1],fQpred_ft_bin_avg_all{i}(k)+fQpred_ft_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{1});
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],fQpred_it_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{2});
        plot(fQ_bin_avg_all{i}(k)*[1 1],fQpred_it_bin_avg_all{i}(k)+fQpred_it_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{2});
        plot(fQ_bin_avg_all{i}(k)+fQ_bin_SE_all{i}(k)*[-1 1],fQpred_hyst_bin_avg_all{i}(k)*[1 1],'Color',PlotColors{3});
        plot(fQ_bin_avg_all{i}(k)*[1 1],fQpred_hyst_bin_avg_all{i}(k)+fQpred_hyst_bin_SE_all{i}(k)*[-1 1],'Color',PlotColors{3});
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
    xlabel('transport activity, $$f_Q$$','interpreter','latex');
    ylabel('predicted activity, $$f_{Q,pred}$$','interpreter','latex');
    set(gca,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
end

%add legend just to last subplot
subplot(1,N_Sites,N_Sites);
legend({'fluid thres.','impact thres.','hysteresis'},'Location','SouthEast');

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 5]);
print([folder_Plots,'fQpred_fQ.png'],'-dpng');