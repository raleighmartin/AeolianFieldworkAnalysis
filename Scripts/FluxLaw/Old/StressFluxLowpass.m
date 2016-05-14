% %% PLOT LOW-PASS TIMESERIES FOR DATA PAIRS WITH SIMILAR STRESS
% 
% %initialize
% clearvars;
% close all;
% 
% %information about where to load data and save plots
% folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for general data files
% folder_Plots = '../PlotOutput/LowPassTimeseries/'; %folder for plots
% 
% %load flux stress window data
% load(strcat(folder_ProcessedData,'StressFluxWindows_all'));
% 
% %get info about number of sites
% N_Sites = length(Sites);

%indices for comparison
indices_comparison = {...
    [];...
    [];...
    [288 388;...
    90 829;...
    169 451;...
    226 840;...
    107 254;...
    221 238;...
    164 281;...
    218 736;...
    174 820;...
    8 206;...
    16 180;...
    49 196;...
    22 746;...
    27 155;...
    152 780;...
    139 772;...
    760 766]...
    };
    
for i = 1:N_Sites
    N_Comparisons = size(indices_comparison{i},1);
    for j = 1:N_Comparisons
        %get info for plotting
        u1 = u_LP_all{i}{indices_comparison{i}(j,1)};
        q1 = q_LP_all{i}{indices_comparison{i}(j,1)};
        u2 = u_LP_all{i}{indices_comparison{i}(j,2)};
        q2 = q_LP_all{i}{indices_comparison{i}(j,2)};

        %get info for labeling plot
        zW1 = zW_all{i}{indices_comparison{i}(j,1)};
        datetime1 = datestr(StartTime_all{i}(indices_comparison{i}(j,1)));
        ust1 = ustRe_all{i}(indices_comparison{i}(j,1));
        ubar1 = u_bar_all{i}(indices_comparison{i}(j,1));
        ustd1 = u_std_all{i}(indices_comparison{i}(j,1));
        eta1 = eta_zs_all{i}(indices_comparison{i}(j,1));
        zL1 = zL_all{i}(indices_comparison{i}(j,1));
        theta1 = theta_all{i}(indices_comparison{i}(j,1));
        Q1 = Q_all{i}(indices_comparison{i}(j,1));
        fQ1 = fQ_all{i}(indices_comparison{i}(j,1));
        
        zW2 = zW_all{i}{indices_comparison{i}(j,2)};
        datetime2 = datestr(StartTime_all{i}(indices_comparison{i}(j,2)));
        ust2 = ustRe_all{i}(indices_comparison{i}(j,2));
        ubar2 = u_bar_all{i}(indices_comparison{i}(j,2));
        ustd2 = u_std_all{i}(indices_comparison{i}(j,2));
        eta2 = eta_zs_all{i}(indices_comparison{i}(j,2));
        zL2 = zL_all{i}(indices_comparison{i}(j,2));
        theta2 = theta_all{i}(indices_comparison{i}(j,2));
        Q2 = Q_all{i}(indices_comparison{i}(j,2));
        fQ2 = fQ_all{i}(indices_comparison{i}(j,2));

        %get info for plot limits
        ulims = [floor(min([u1; u2])) ceil(max([u1;u2]))];
        uticks = floor(min([u1; u2])):ceil(max([u1;u2]));
        qmax = max([q1(:,1); q2(:,1)]);
        qstep = 10^(floor(log10(qmax)));
        qmax = ceil(qmax/qstep)*qstep;
        qlims = [0 qmax];
        qticks = 0:qstep:qmax;
        
        %plot first comparison
        figure(j); clf;
        if Q1>Q2;
            subplot(2,1,1); %on top panel if larger
        elseif Q2>Q1;
            subplot(2,1,2); %on bottom panel is smaller
        end
        [ax1,p11,p12] = plotyy(t_u_LP_all{i},u1,t_q_LP_all{i},q1(:,1));
        p12.LineWidth = 1;
        set(ax1(1),'ylim',ulims);
        set(ax1(1),'ytick',uticks);
        set(ax1(2),'ylim',qlims);
        set(ax1(2),'ytick',qticks);
        xlabel('t (s)');
        ylabel(ax1(1),'$$\hat{u}$$ (m/s)','Interpreter','Latex');
        ylabel(ax1(2),'$$\hat{q}$$ (g/m$$^2$$/s)','Interpreter','Latex');
        title([datetime1,...
            ', u_{*} = ',num2str(round(ust1,3)),...
            ' m/s, \mu_{u} = ',num2str(round(ubar1,2)),...
            ' m/s, \sigma_{u} = ',num2str(round(ustd1,2)),...
            ' m/s, \eta = ',num2str(round(eta1,2)),...
            ' , z/L = ',num2str(round(zL1,3)),...
            ' , \theta = ',num2str(round(theta1)),...
            '^{\circ}, Q = ',num2str(round(Q1,2)),...
            ' g/m/s, f_{Q} = ',num2str(round(fQ1,2)),...
            ' , z_W = ',num2str(round(zW1(1)*100,1)),' cm']);

        %plot second comparison
        if Q1>Q2;
            subplot(2,1,2); %on top panel if larger
        elseif Q2>Q1;
            subplot(2,1,1); %on bottom panel if smaller
        end
        [ax2,p21,p22] = plotyy(t_u_LP_all{i},u2,t_q_LP_all{i},q2(:,1));
        p22.LineWidth = 1;
        set(ax2(1),'ylim',ulims);
        set(ax2(1),'ytick',uticks);
        set(ax2(2),'ylim',qlims);
        set(ax2(2),'ytick',qticks);
        xlabel('t (s)');
        ylabel(ax2(1),'$$\hat{u}$$ (m/s)','Interpreter','Latex');
        ylabel(ax2(2),'$$\hat{q}$$ (g/m$$^2$$/s)','Interpreter','Latex');
        title([datetime2,...
            ', u_{*} = ',num2str(round(ust2,3)),...
            ' m/s, \mu_{u} = ',num2str(round(ubar2,2)),...
            ' m/s, \sigma_{u} = ',num2str(round(ustd2,2)),...
            ' m/s, \eta = ',num2str(round(eta2,2)),...
            ' , z/L = ',num2str(round(zL2,3)),...
            ' , \theta = ',num2str(round(theta2)),...
            '^{\circ}, Q = ',num2str(round(Q2,2)),...
            ' g/m/s, f_{Q} = ',num2str(round(fQ2,2)),...
            ' , z_W = ',num2str(round(zW2(1)*100,1)),' cm']);
        
        set(gcf, 'PaperPosition',[0 0 11 6]);
        print([folder_Plots,'StressFluxLowpass_',Sites{i},'_',int2str(j),'.png'],'-dpng');
    end
end
        

