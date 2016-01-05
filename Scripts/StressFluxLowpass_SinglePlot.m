%% PLOT LOW-PASS TIMESERIES FOR DATA PAIRS WITH SIMILAR STRESS

%initialize
% clearvars;
% close all;

%information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
folder_Plots = '../PlotOutput/LowPassTimeseries/'; %folder for plots

%load flux stress window data
load(strcat(folder_ProcessedData,'StressFluxWindows_all'));

%get info about number of sites
N_Sites = length(Sites);

%indices for plot
indices_plot = {[];...
    [];...
    [180]};
    %[8, 16, 22, 27, 49, 107, 139, 152, 155, 164, 169, 174, 180, 196, 206, 218, 221, 238, 281, 288, 388]};
    
for i = 1:N_Sites
    N_Plots = length(indices_plot{i});
    for j = 1:N_Plots
        %get info for plotting
        u = u_LP_all{i}{indices_plot{i}(j)};
        q = q_LP_all{i}{indices_plot{i}(j)};

        %get info for labeling plot
        zW = zW_all{i}{indices_plot{i}(j)};
        zU = zU_all{i}(indices_plot{i}(j));
        starttime = datestr(StartTime_all{i}(indices_plot{i}(j)));
        ust = ustRe_all{i}(indices_plot{i}(j));
        ubar = u_bar_all{i}(indices_plot{i}(j));
        ustd = u_std_all{i}(indices_plot{i}(j));
        eta = eta_active_fQ_1s_TFEMthr_all{i}(indices_plot{i}(j));
        zL = zL_all{i}(indices_plot{i}(j));
        theta = theta_all{i}(indices_plot{i}(j));
        Q = Q_all{i}(indices_plot{i}(j));
        fQ = fQ_all{i}(indices_plot{i}(j));
        
        %get info for plot limits
        %ulims = [floor(min(u)) ceil(max(u))];
        %uticks = floor(min(u)):ceil(max(u));
        ulims = [0 ceil(max(u))];
        uticks = 0:ceil(max(u));
        qmax = max(q(:,1));
        qstep = 10^(floor(log10(qmax)));
        qmax = ceil(qmax/qstep)*qstep;
        qlims = [0 qmax];
        qticks = 0:qstep:qmax;
        
        %plot comparison
        figure(j); clf; hold on;
        [ax1,p11,p12] = plotyy(t_u_LP_all{i},u,t_q_LP_all{i},q(:,1));
        p11.LineWidth = 1;
        p12.LineWidth = 1;
        set(ax1(1),'ylim',ulims);
        set(ax1(1),'ytick',uticks);
        set(ax1(2),'ylim',qlims);
        set(ax1(2),'ytick',qticks);
        xlabel('t (s)');
        ylabel(ax1(1),'$$\hat{u}$$ (m/s)','Interpreter','Latex','FontSize',16);
        ylabel(ax1(2),'$$\hat{q}$$ (g/m$$^2$$/s)','Interpreter','Latex','FontSize',16);
        title([Sites{i}, ', ', starttime,...
            ', u_{*} = ',num2str(round(ust,3)),' m/s']);
        set(gca,'FontSize',16);
        set(gcf, 'PaperPosition',[0 0 6 6]);
        print([folder_Plots,'StressFluxLowpassSample_',Sites{i},'_',int2str(j),'.pdf'],'-dpdf');
    end
end