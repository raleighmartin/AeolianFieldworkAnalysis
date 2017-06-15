%% INITIALIZATION
%initialize
clearvars;
close all;

%% PARAMETERS AND INPUTS
%% flux profile interval times
T_subwindow = [...
%    duration(0,0,0.04),...    
%    duration(0,0,0.08),...    
%    duration(0,0,0.12),...    
%    duration(0,0,0.2),...    
%    duration(0,0,0.4),...    
%    duration(0,0,0.6),...    
    duration(0,0,1),...    
    duration(0,0,2),...    
    duration(0,0,3),...
    duration(0,0,5),...
    duration(0,0,10),...
    duration(0,0,15),...
    duration(0,0,30),...
    duration(0,1,0),...
    duration(0,2,0),...
    duration(0,3,0),...
    duration(0,5,0),...
    duration(0,10,0),...
    duration(0,15,0),...
    duration(0,30,0)
    ];
N_T_subwindow = length(T_subwindow); %number of subwindow durations
T_subwindow_s = seconds(T_subwindow); %duration of subwindows

%% other parameter values
T_window = duration(0,30,0);
kappa = 0.4; %von Karman parameter
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
g = 9.8; %gravity m/s^2
fD_min = 0.005; %minimum detection rate for zero flux
Q_min = 0.05; %detection limit for Q, set to zero if below this value
zq_Q_min = 0.10; %assumed saltation height for detection limit for exponential profile for detection limit for individual Wenglor
zW_limit = 3; %limit on number of unique Wenglor heights in profile

%set info for plotting
Markers_Field = {'s','d','o'};
Colors_Field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
PlotFont = 14;

%% LOAD DATA AND FUNCTIONS
%folders for loading data, saving data, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_DataBSNE = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder containing BSNE data
folder_SaveData = '../../AnalysisData/Methods/'; %folder for data output
folder_Plots = '../../PlotOutput/Methods/'; %folder containing plot output
folder_Functions = '../Functions/'; %folder with functions

%paths for loading and saving data
LoadData_Path = strcat(folder_LoadData,'DataWindows_30min_Restricted'); %path for 30 minute data
SaveData_Path = strcat(folder_SaveData,'FluxSubwindows_30min_Restricted'); %path for 30 minute data

%load data
load(LoadData_Path); %load window data

%load functions
addpath(folder_Functions); %point MATLAB to location of functions

%% initialize variable arrays
timeofday_subwindow = cell(N_Sites,1); %time of day array
StartTime_subwindow = cell(N_Sites,1); %start time array
EndTime_subwindow = cell(N_Sites,1); %end time array
ind_window_subwindow = cell(N_Sites,1); %index of corresponding window for subwindow

%% initialize wind arrays
zU_subwindow = cell(N_Sites,1); %anemometer height for subwindow
ubar_subwindow = cell(N_Sites,1); %horizontal winds for subwindow
vbar_subwindow = cell(N_Sites,1); %lateral winds for subwindow
wbar_subwindow = cell(N_Sites,1); %vertical winds for subwindow
theta_subwindow = cell(N_Sites,1); %wind angle for subwindow

%% initialize flux arrays
t_flux_subwindow = cell(N_Sites,1); %times for subwindow with flux and wind
qbar_subwindow = cell(N_Sites,1); %partial flux for subwindow
sigma_qbar_subwindow = cell(N_Sites,1); %partial flux uncertainty for subwindow
Qfit_subwindow = cell(N_Sites,1); %total flux by fitting for subwindow
sigma_Qfit_subwindow = cell(N_Sites,1); %total flux by fitting uncertainty for subwindow
zq_subwindow = cell(N_Sites,1); %flux height by fitting for subwindow
sigma_zq_subwindow = cell(N_Sites,1); %uncertainty in flux height by fitting for subwindow
Chi2_Qfit_subwindow = cell(N_Sites,1); %goodness of Q fit
df_Qfit_subwindow = cell(N_Sites,1); %degrees of freedom for fitting
Qsum_subwindow = cell(N_Sites,1); %total flux by summing for subwindow
sigma_Qsum_subwindow = cell(N_Sites,1); %total flux by summing uncertainty for subwindow

%% GO THROUGH SITES
for i = 1:N_Sites

    %% load BSNE data
    load(strcat(folder_DataBSNE,'FluxBSNE_',Sites{i})); 
    
    %% get information about windows and subwindows
    ind_windows = find(hasfluxdata_window{i}==1); %get indices of subwindows with flux
    N_windows = length(ind_windows); %get number of windows
    
    %% initialize variable arrays - times
    timeofday_subwindow{i} = cell(N_T_subwindow,1); %time of day array
    StartTime_subwindow{i} = cell(N_T_subwindow,1); %start time array
    EndTime_subwindow{i} = cell(N_T_subwindow,1); %end time array
    ind_window_subwindow{i} = cell(N_T_subwindow,1); %index of corresponding window for subwindow
    
    %% initialize variable arrays - wind
    zU_subwindow{i} = cell(N_T_subwindow,1); %anemometer height for subwindow
    ubar_subwindow{i} = cell(N_T_subwindow,1); %wind timeseries array
    vbar_subwindow{i} = cell(N_T_subwindow,1); %lateral winds for subwindow
    wbar_subwindow{i} = cell(N_T_subwindow,1); %vertical wind timeseries array
    theta_subwindow{i} = cell(N_T_subwindow,1); %wind angle for subwindow

    %% initialize variable arrays - flux
    qbar_subwindow{i} = cell(N_T_subwindow,1); %partial flux for subwindow
    sigma_qbar_subwindow{i} = cell(N_T_subwindow,1); %partial flux uncertainty for subwindow
    Qfit_subwindow{i} = cell(N_T_subwindow,1); %total flux by fitting for subwindow
    sigma_Qfit_subwindow{i} = cell(N_T_subwindow,1); %total flux by fitting uncertainty for subwindow
    Chi2_Qfit_subwindow{i} = cell(N_T_subwindow,1); %goodness of Q fit
    df_Qfit_subwindow{i} = cell(N_T_subwindow,1); %degrees of freedom for fitting
    zq_subwindow{i} = cell(N_T_subwindow,1); %flux height by fitting for subwindow
    sigma_zq_subwindow{i} = cell(N_T_subwindow,1); %uncertainty in flux height by fitting for subwindow
    Qsum_subwindow{i} = cell(N_T_subwindow,1); %total flux by summing for subwindow
    sigma_Qsum_subwindow{i} = cell(N_T_subwindow,1); %total flux by summing uncertainty for subwindow
    
    %% GO THROUGH MEASUREMENT INTERVALS
    for m = 1:N_T_subwindow

        %% get information about subwindows
        N_subwindows_per_window = T_window/T_subwindow(m);
        N_subwindows = N_windows*N_subwindows_per_window; %get number of subwindows
        
        %% display processing status
        processing_status = [Sites{i},', ',int2str(m),' of ',int2str(N_T_subwindow),', ',datestr(now)]

        %% initialize subwindow values - time info
        timeofday_subwindow{i}{m} = zeros(N_subwindows,1); %time of day array
        StartTime_subwindow{i}{m} = datetime(zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1)); %start time array
        EndTime_subwindow{i}{m} = datetime(zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1),zeros(N_subwindows,1)); %end time array            
        ind_window_subwindow{i}{m} = zeros(N_subwindows,1); %array of associated window numbers
        
        %% initialize wind values
        zU_subwindow{i}{m} = zeros(N_subwindows,1); %anemometer height
        ubar_subwindow{i}{m} = zeros(N_subwindows,1); %wind timeseries array
        vbar_subwindow{i}{m} = zeros(N_subwindows,1); %transverse wind timeseries array
        wbar_subwindow{i}{m} = zeros(N_subwindows,1); %vertical wind timeseries array
        theta_subwindow{i}{m} = zeros(N_subwindows,1); %wind angle for subwindow

        %% initialize flux values
        qbar_subwindow{i}{m} = cell(N_subwindows,1); %partial flux for subwindow
        sigma_qbar_subwindow{i}{m} = cell(N_subwindows,1); %partial flux uncertainty for subwindow
        Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by fitting for subwindow
        sigma_Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by fitting uncertainty for subwindow
        Chi2_Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %goodness of Q fit
        df_Qfit_subwindow{i}{m} = zeros(N_subwindows,1); %degrees of freedom for fitting
        zq_subwindow{i}{m} = zeros(N_subwindows,1); %flux height by fitting for subwindow
        sigma_zq_subwindow{i}{m} = zeros(N_subwindows,1); %uncertainty in flux height by fitting for subwindow
        Qsum_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by summing for subwindow
        sigma_Qsum_subwindow{i}{m} = zeros(N_subwindows,1); %total flux by summing uncertainty for subwindow
    
        %% GO THROUGH 30-MINUTE WINDOWS
        for j = 1:N_windows
                        
            %% get information about BSNE
            ind_BSNE = find([FluxBSNE.StartTime] <= StartTime_window{i}(ind_windows(j)) & [FluxBSNE.EndTime] >= EndTime_window{i}(ind_windows(j)));
            zq_BSNE = FluxBSNE(ind_BSNE).z.zq; %get best fit saltation layer height for BSNE profile
            
            %% get wind values for window
            zU = zU_base_window{i}(ind_windows(j)); %anemometer height
            t_wind = t_wind_int_window{i}{ind_windows(j)}; %times for interpolated wind
            u = u_int_window{i}{ind_windows(j)}; %rotated interpolated streamwise velocities
            v = v_int_window{i}{ind_windows(j)}; %rotated interpolated lateral winds
            w = w_int_window{i}{ind_windows(j)}; %rotated interpolated vertical velocities

            %% get flux values for window
            t_flux = t_flux_int_window{i}{ind_windows(j)}; %times for interpolated flux
            zW = zW_window{i}{ind_windows(j)}; %get Wenglor heights
            sigma_zW = sigma_zW_window{i}{ind_windows(j)}; %get uncertainty in Wenglor heights
            q = q_int_window{i}{ind_windows(j)}; %get partial fluxes
            Cqn = Cqn_int_window{i}{ind_windows(j)}; %get calibration coefficients
            sigma_Cqn = sigma_Cqn_int_window{i}{ind_windows(j)}; %get uncertainties in calibration coefficients
            n = n_int_window{i}{ind_windows(j)}; %get counts rates
            
            %% get subwindow times within window
            window_StartTime_subwindow = StartTime_window{i}(ind_windows(j))+((1:N_subwindows_per_window)-1)*T_subwindow(m);
            window_EndTime_subwindow = EndTime_window{i}(ind_windows(j))+(1:N_subwindows_per_window)*T_subwindow(m);

            %% go through subwindows
            for k = 1:N_subwindows_per_window
                
                %% get index for adding values to array
                subwindow_array_ind = (j-1)*N_subwindows_per_window+k;

                %% get indices of all subwindow points within window
                subwindow_wind_ind = find(t_wind>=window_StartTime_subwindow(k)&t_wind<window_EndTime_subwindow(k));
                subwindow_flux_ind = find(t_flux>=window_StartTime_subwindow(k)&t_flux<window_EndTime_subwindow(k));

                %% get time of day, start time, end time, and other basic info about subwindow
                StartTime_subwindow{i}{m}(subwindow_array_ind) = window_StartTime_subwindow(k); %start time
                EndTime_subwindow{i}{m}(subwindow_array_ind) = window_EndTime_subwindow(k); %start time
                ind_window_subwindow{i}{m}(subwindow_array_ind) = ind_windows(j); %window corresponding to subwindow
                zU_subwindow{i}{m}(subwindow_array_ind) = zU; %anemometer height - from 30-minute window

                %% get wind data in subwindow
                ubar_subwindow{i}{m}(subwindow_array_ind) = mean(u(subwindow_wind_ind)); %winds for subwindow
                vbar_subwindow{i}{m}(subwindow_array_ind) = mean(v(subwindow_wind_ind)); %lateral winds for subwindow
                wbar_subwindow{i}{m}(subwindow_array_ind) = mean(w(subwindow_wind_ind)); %vertical winds for subwindow
                theta_subwindow{i}{m}(subwindow_array_ind) = 180/pi*atan(mean(v(subwindow_wind_ind))./mean(u((subwindow_wind_ind))));
                timeofday_subwindow{i}{m}(subwindow_array_ind) = hour(window_StartTime_subwindow(k))+minute(window_StartTime_subwindow(k))/60+second(window_StartTime_subwindow(k))/3600;
               
                %% get flux data in subwindow
                q_subwindow = q(subwindow_flux_ind,:); %get partial fluxes
                n_subwindow = n(subwindow_flux_ind,:); %get counts rates
                
                %get mean q profile information for fitting
                qbar = mean(q_subwindow); %compute mean qz profile
                nbar = mean(n_subwindow)/dt_flux_window(i); %compute mean particle counts per second for each qz
                sigma_nbar = sqrt(nbar/T_subwindow_s(m)); %compute uncertainty in mean particle counts
                Cqnbar = mean(Cqn); %compute mean calibration coefficient for each qz
                sigma_Cqnbar = mean(sigma_Cqn); %compute mean calibration coefficient uncertainty for each qz
                sigma_qbar = sqrt((sigma_Cqnbar.*nbar).^2+(sigma_nbar.*Cqnbar).^2); %compute uncertainty in qz profile, from counts and calibration coefficient

                %deal with repeated values
                zW_unique = unique(zW);
                N_zW_unique = length(zW_unique);
                qbar_unique = zeros(size(zW_unique));
                sigma_qbar_unique = zeros(size(zW_unique));
                sigma_zW_unique = zeros(size(zW_unique));
                for s = 1:N_zW_unique
                    ind_zW = find(zW==zW_unique(s));
                    if length(ind_zW)>1 %compute mean and uncertainty for repeated heights
                        q_min = (Q_min/zW_unique(s))*exp(-zW_unique(s)/zq_Q_min); %get min q detection limit for this height
                        qbar_z = mean(qbar(ind_zW)); %get mean q for this height
                        if qbar_z < q_min %if mean q at height is below detection limit
                            qbar_z = 0; sigma_qbar_z = 0; %then just set values to zero
                        else %otherwise, use script "MeanUncertainty.m"
                            [qbar_z, sigma_qbar_z] = MeanUncertainty(qbar(ind_zW), sigma_qbar(ind_zW)); %compute flux mean and uncertainy for repeated heights
                        end
                        sigma_zW_unique(s) = mean(sigma_zW(ind_zW)); %no reduction in uncertainty for height, because all uncertainties are correlated, so just take mean of values
                        qbar_unique(s) = qbar_z;
                        sigma_qbar_unique(s) = sigma_qbar_z;
                    else %otherwise, if only Wenglor at given height, just use existing values
                        qbar_unique(s) = qbar(ind_zW);
                        sigma_qbar_unique(s) = sigma_qbar(ind_zW);
                        sigma_zW_unique(s) = sigma_zW(ind_zW);
                    end
                end

                %Remove 0 values for fitting
                ind_fit = find(qbar_unique>0);
                N_fit = length(ind_fit);
                qbar_fit = qbar_unique(ind_fit);
                zW_fit = zW_unique(ind_fit);
                sigma_qbar_fit = sigma_qbar_unique(ind_fit);
                sigma_zW_fit = zeros(1,N_fit); %neglect uncertainty in Wenglor height, which is already accounted for by calibration

                %Perform profile fit to get q0, zq, and Q if sufficient points for fitting
                if N_fit>=zW_limit
                    [q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,q_pred,sigma_q_pred] = qz_profilefit(qbar_fit,zW_fit,sigma_qbar_fit,sigma_zW_fit);
                    q_residuals = q_pred - qbar_fit; %residuals between observed and predicted q
                    Chi2_Qfit = sum((q_residuals./sigma_q_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
                    df_Qfit = N_fit-2; %compute degrees of freedom for Qfit
                else %otherwise, set to NaN
                    Q=NaN;
                    sigma_Q = NaN;
                    zq=NaN;
                    sigma_q0=NaN;
                    sigma_zq=NaN;
                    Chi2_Qfit=NaN;
                    df_Qfit=NaN;
                end
                
                %Perform summation to get Q
                %compute delta_z weights for additive flux calculations
                z1_Qsum = [0 sqrt(zW_unique(1:end-1).*zW_unique(2:end))]; %bottom of summation bins
                z2_Qsum = [sqrt(zW_unique(1:end-1).*zW_unique(2:end)) Inf];
                q0_Qsum = qbar_unique./(exp(-zW_unique/zq_BSNE));
                sigma_q0_Qsum = sigma_qbar_unique./(exp(-zW_unique/zq_BSNE));
                deltaQ = q0_Qsum.*zq_BSNE.*(exp(-z1_Qsum/zq_BSNE)-exp(-z2_Qsum/zq_BSNE));
                sigma_deltaQ = sigma_q0_Qsum.*zq_BSNE.*(exp(-z1_Qsum/zq_BSNE)-exp(-z2_Qsum/zq_BSNE));
                deltaz = deltaQ./qbar_unique; %equivalent deltaz for summation
                Qsum = sum(deltaQ);
                sigma_Qsum = sqrt(sum(sigma_deltaQ.^2));
                
                qbar_subwindow{i}{m}{subwindow_array_ind} = qbar; %partial flux for subwindow
                sigma_qbar_subwindow{i}{m}{subwindow_array_ind} = sigma_qbar; %partial flux uncertainty for subwindow
                Qfit_subwindow{i}{m}(subwindow_array_ind) = Q; %total flux by fitting for subwindow
                sigma_Qfit_subwindow{i}{m}(subwindow_array_ind) = sigma_Q; %total flux by fitting uncertainty for subwindow
                Chi2_Qfit_subwindow{i}{m}(subwindow_array_ind) = Chi2_Qfit; %total flux by fitting uncertainty for subwindow
                df_Qfit_subwindow{i}{m}(subwindow_array_ind) = df_Qfit; %degrees of freedom for fitting
                zq_subwindow{i}{m}(subwindow_array_ind) = zq; %flux height by fitting for subwindow
                sigma_zq_subwindow{i}{m}(subwindow_array_ind) = sigma_zq; %uncertainty in flux height by fitting for subwindow
                Qsum_subwindow{i}{m}(subwindow_array_ind) = Qsum; %total flux by summing for subwindow
                sigma_Qsum_subwindow{i}{m}(subwindow_array_ind) = sigma_Qsum; %total flux by summing uncertainty for subwindow
            end           
        end
    end
end

%% save analysis data
save(SaveData_Path,'SiteNames','Sites','N_Sites','*subwindow'); 

%% load analysis data
load(SaveData_Path); %load saved data
addpath(folder_Functions); %point MATLAB to location of functions


%% plot comparison of flux estimates
for m = 1:N_T_subwindow
    
    %% initialize figure for calculation timescale
    figure(m); clf; hold on;
    
    %% plot comparison of flux estimates
    for i = 1:N_Sites
        plot(Qsum_subwindow{i}{m},Qfit_subwindow{i}{m},Markers_Field{i},'Color',Colors_Field{i})
    end
    
    %% determine maximum Q
    Qmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
    for i = 1:N_Sites
        Qmax(i) = ceil(max(max(Qsum_subwindow{i}{m}),max(Qfit_subwindow{i}{m}))/10)*10;
    end
    Qmax = max(Qmax);
    
    %% plot 1-1 line
    plot([0 Qmax],[0 Qmax],'k--');
    
    %% format plot
    xlim([0 Qmax]);
    ylim([0 Qmax]);
    xlabel('$$Q_{sum}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
    ylabel('$$Q_{fit}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
    title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
    legend(SiteNames,'Location','NorthWest');
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    set(gca,'FontSize',PlotFont);
    
    %% print plot
    set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
    print([folder_Plots,'Q_comparison_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
end

%% compute relative difference of Q estimates
Qrel_subwindow = cell(N_Sites,1);
median_Qrel_subwindow = cell(N_Sites,1);
for i = 1:N_Sites
    Qrel_subwindow{i} = cell(N_T_subwindow,1);
    median_Qrel_subwindow{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        Qrel_subwindow{i}{m} = abs(Qfit_subwindow{i}{m}-Qsum_subwindow{i}{m})./Qsum_subwindow{i}{m};
        median_Qrel_subwindow{i}(m) = median(Qrel_subwindow{i}{m}(~isnan(Qrel_subwindow{i}{m})));
    end
end

%% plot comparison of relative difference in flux estimates
for m = 1:N_T_subwindow
    
    %% initialize figure for calculation timescale
    figure(m); clf; hold on;
    
    %% plot comparison of flux estimates
    for i = 1:N_Sites
        ind_nonzero = find(Qsum_subwindow{i}{m}>=Q_min);
        plot(Qsum_subwindow{i}{m}(ind_nonzero),100*Qrel_subwindow{i}{m}(ind_nonzero),Markers_Field{i},'Color',Colors_Field{i})
    end
    
    %% determine maximum Q
    Qmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
    for i = 1:N_Sites
        Qmax(i) = ceil(max(Qsum_subwindow{i}{m}/10))*10;
    end
    Qmax = max(Qmax); 
        
    %% format plot
    xlim([0 Qmax]);
    xlabel('$$Q_{sum}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
    ylabel('$$|Q_{fit}-Q_{sum}|/Q_{sum}$$ ($$\%$$)','Interpreter','LaTeX');
    title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
    legend(SiteNames,'Location','NorthEast');
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    set(gca,'FontSize',PlotFont);
    
    %% print plot
    set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
    print([folder_Plots,'Q_comparison_relative_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
end

%% plot median relative difference
figure(N_T_subwindow+1); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, 100*median_Qrel_subwindow{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$|Q_{fit}-Q_{sum}|/Q_{sum}$$ ($$\%$$)','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'Q_diff_T.png'],'-dpng'); %for draft

%% recompute df
for i = 1:N_Sites
    for m = 1:N_T_subwindow
        for j = 1:length(qbar_subwindow{i}{m})
            df_Qfit_subwindow{i}{m}(j) = length(qbar_subwindow{i}{m}{j});
        end
    end
end

%% compute median chi2nu for fit
Chi2nu_subwindow = cell(N_Sites,1);
median_Chi2nu_subwindow = cell(N_Sites,1);
for i = 1:N_Sites
    Chi2nu_subwindow{i} = cell(N_T_subwindow,1);
    median_Chi2nu_subwindow{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        Chi2nu_subwindow{i}{m} = Chi2_Qfit_subwindow{i}{m}./df_Qfit_subwindow{i}{m};
        median_Chi2nu_subwindow{i}(m) = median(Chi2nu_subwindow{i}{m}(~isnan(Chi2nu_subwindow{i}{m})));
    end
end

%% compute median zq
median_zq_subwindow = cell(N_Sites,1);
for i = 1:N_Sites
    median_zq_subwindow{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        zq_all = zq_subwindow{i}{m}(~isnan(qz_subwindow{i}{m}));
        median_zq_subwindow{i}{m} = median(zq_all);
    end
end

%% compute median sigma_q / q
sigmaq_rel_subwindow = cell(N_Sites,1);
median_sigmaq_rel_subwindow = cell(N_Sites,1);
for i = 1:N_Sites
    sigmaq_rel_subwindow{i} = cell(N_T_subwindow,1);
    median_sigmaq_rel_subwindow{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        N_subwindow = length(sigma_qbar_subwindow{i}{m});
        sigmaq_rel_subwindow{i}{m} = zeros(N_subwindow,1);
        for j = 1:N_subwindow
            sigmaq_rel_subwindow{i}{m}(j) = mean(sigma_qbar_subwindow{i}{m}{j}./qbar_subwindow{i}{m}{j});
        end
        median_sigmaq_rel_subwindow{i}(m) = median(sigmaq_rel_subwindow{i}{m}(~isnan(sigmaq_rel_subwindow{i}{m})));
    end
end




%% compute fraction NaN for fit
fNaN_subwindow = cell(N_Sites,1);
for i = 1:N_Sites
    fNaN_subwindow{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        fNaN_subwindow{i}(m) = length(find(isnan(Chi2nu_subwindow{i}{m})))/length(Chi2nu_subwindow{i}{m});
    end
end

%% plot median Chi2nu
figure(N_T_subwindow+2); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, median_Chi2nu_subwindow{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$\chi^{2}_{\nu}$$','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'Chi2nu_T.png'],'-dpng'); %for draft

%% plot median_sigmaq_rel
figure(N_T_subwindow+3); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, median_sigmaq_rel_subwindow{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('median $$\sigma_{q}/{q}$$','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'sigmaq_rel_T.png'],'-dpng'); %for draft


%% plot percent NaN for subwindow
figure(N_T_subwindow+4); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, 100*fNaN_subwindow{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('percent NaN for profile fit','Interpreter','LaTeX');
legend(SiteNames,'Location','NorthOutside');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
ylims = ylim;
ylim([0 ylims(2)]);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
print([folder_Plots,'PercentNaN_T.png'],'-dpng'); %for draft

%% fit to Q versus u2
Chi2nu_Q_u2 = cell(N_Sites,1);
uth_Q_u2 = cell(N_Sites,1);
C_Q_u2 = cell(N_Sites,1);
for i = 1:N_Sites
    Chi2nu_Q_u2{i} = zeros(N_T_subwindow,1);
    uth_Q_u2{i} = zeros(N_T_subwindow,1);
    C_Q_u2{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        ind_fit = find(Qsum_subwindow{i}{m}>1);
        u2_fit = ubar_subwindow{i}{m}(ind_fit).^2;
        Q_fit = Qsum_subwindow{i}{m}(ind_fit);
        [a, b, sigma_a, sigma_b, Q_pred, sigma_Q_pred] = linearfit(u2_fit, Q_fit);
        uth = sqrt(-a/b); C = b;
        Q_residuals = Q_pred - Q_fit; %residuals between observed and predicted q
        Chi2 = sum((Q_residuals./sigma_Q_pred).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
        df = length(ind_fit)-2; %compute degrees of freedom for Qfit
        Chi2nu_Q_u2{i}(m) = Chi2/df;
        uth_Q_u2{i}(m) = uth;
        C_Q_u2{i}(m) = C;
    end
end

%% plot Q versus u for different timescales
for m = 1:N_T_subwindow
    
    %% initialize figure for calculation timescale
    figure(m); clf; hold on;
  
    %% determine maximum u
    umax = zeros(N_Sites,1); %compute for each site, then take max of all sites
    for i = 1:N_Sites
        umax(i) = ceil(max(ubar_subwindow{i}{m}/10))*10;
    end
    umax = max(umax);
    
    %% plot comparison of flux estimates
    for i = 1:N_Sites
        plot(ubar_subwindow{i}{m},Qsum_subwindow{i}{m},Markers_Field{i},'Color',Colors_Field{i})
    end

    %% plot fit
    for i = 1:N_Sites
        uth = uth_Q_u2{i}(m);
        C = C_Q_u2{i}(m);
        u_fit = linspace(uth,umax,100);
        Q_fit = C*(u_fit.^2-uth.^2);
        plot(u_fit,Q_fit,'Color',Colors_Field{i},'LineWidth',2)
    end
        
    %% determine maximum Q
    Qmax = zeros(N_Sites,1); %compute for each site, then take max of all sites
    for i = 1:N_Sites
        Qmax(i) = ceil(max(Qsum_subwindow{i}{m}/10))*10;
    end
    Qmax = max(Qmax);  
       
    %% format plot
    ylim([0 Qmax]);
    xlabel('$$u$$ (m s$$^{-1}$$)','Interpreter','LaTeX');
    ylabel('$$Q_{sum}$$ (g m$$^{-1}$$ s$$^{-1}$$)','Interpreter','LaTeX');
    title(['T_{HF} = ',num2str(T_subwindow_s(m)),' s']);
    legend(SiteNames,'Location','NorthWest');
    set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
    set(gca,'FontSize',PlotFont);
    
    %% print plot
    set(gcf,'PaperUnits','inches','PaperSize',[7.3 6.5],'PaperPosition',[0 0 7.3 6.5],'PaperPositionMode','Manual');
    print([folder_Plots,'Q_u_',num2str(T_subwindow_s(m)),'.png'],'-dpng'); %for draft
end

%% plot parameters for Q versus u^2 fit
figure(N_T_subwindow+4); clf;

%% plot Chi2nu
subplot('Position',[0.15 0.71 0.8 0.25]); hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, Chi2nu_Q_u2{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
%xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('$$\chi^{2}_{\nu}$$ for $$Q$$ vs $$u^2$$ fit','Interpreter','LaTeX');
legend(SiteNames,'Location','SouthWest');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log','YScale','Log');
set(gca,'FontSize',PlotFont);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% plot uth
subplot('Position',[0.15 0.39 0.8 0.25]); hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, uth_Q_u2{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
%xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('$$u_{th}$$ (m s$$^{-1}$$) for $$Q$$ vs $$u^2$$ fit','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% plot C
subplot('Position',[0.15 0.07 0.8 0.25]); hold on;
for i = 1:N_Sites
    plot(T_subwindow_s, C_Q_u2{i},Markers_Field{i},'Color',Colors_Field{i});
end

%% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('$$C$$ (g m$$^{-2}$$) for $$Q$$ vs $$u^u$$ fit','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On','XScale','Log');
set(gca,'FontSize',PlotFont);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% print plot
set(gcf,'PaperUnits','inches','PaperSize',[4 9],'PaperPosition',[0 0 5 9],'PaperPositionMode','Manual');
print([folder_Plots,'Q_u_fit.png'],'-dpng'); %for draft

%% compute standard deviation of zq versus averaging time
std_zq_subwindow = cell(N_Sites,1);
for i = 1:N_Sites
    std_zq_subwindow{i} = zeros(N_T_subwindow,1);
    for m = 1:N_T_subwindow
        std_zq_subwindow{i}(m) = std(zq_subwindow{i}{m}(~isnan(zq_subwindow{i}{m})));
    end
end

%% plot standard deviation of zq versus averaging time
figure(N_T_subwindow+5); clf; hold on;
for i = 1:N_Sites
    plot(T_subwindow_s,std_zq_subwindow{i},Markers_Field{i},'Color',Colors_Field{i})
end

%% format plot
xlabel('$$T_{HF}$$ (s)','Interpreter','LaTeX');
ylabel('std. dev. in $$z_q$$ from fit (m)','Interpreter','LaTeX');
set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
set(gca,'yscale','log','xscale','log');
set(gca,'FontSize',PlotFont);
legend(Sites);
xlim([min(T_subwindow_s) max(T_subwindow_s)]);

%% print plot
set(gcf,'PaperUnits','inches','PaperSize',[7 4],'PaperPosition',[0 0 7 4],'PaperPositionMode','Manual');
print([folder_Plots,'std_zq_T.png'],'-dpng'); %for draft