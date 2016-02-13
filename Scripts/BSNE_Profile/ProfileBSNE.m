%% initialize
clearvars;

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/StressFlux/'; %folder for outputs
folder_Plots = '../../PlotOutput/BSNE_Profile/'; %folder for plots

%% load BSNE data
load(strcat(folder_AnalysisData,'FluxBSNE'));

%% Plot all profiles
for ind_Site = 1:length(FluxBSNE)
    for ind_Interval = 1:length(FluxBSNE{ind_Site})

        %% Get data for profile
        FluxInterval = FluxBSNE{ind_Site}(ind_Interval);

        z = FluxInterval.z.z;
        q = FluxInterval.qz.qz;
        sigma_q = FluxInterval.qz.sigma;
        q_pred = FluxInterval.qz.qz_pred;
        sigma_q_pred = FluxInterval.qz.sigma_qz_pred;
        q_pred_minus = q_pred - sigma_q_pred;
        q_pred_plus = q_pred + sigma_q_pred;
        [z_sort, ind_sort] = sort(z);
        q_pred_sort = q_pred(ind_sort);
        q_pred_minus_sort = q_pred_minus(ind_sort);
        q_pred_plus_sort = q_pred_plus(ind_sort);

        %make plot
        figure(1); clf;
        for j=1:2;
            subplot(2,1,j);
            errorbar(z,q,sigma_q,'bx','MarkerSize',8); hold on;
            plot(z_sort,q_pred_sort,'k','LineWidth',1);
            plot(z_sort,q_pred_minus_sort,'r--',z_sort,q_pred_plus_sort,'r--');
            xlims = xlim;
            xlim([0 xlims(2)]);
            set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
            xlabel('\textbf{BSNE trap height, $$z_{B}$$ (m)}','Interpreter','Latex');
            ylabel('\textbf{BSNE flux, $$q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
            set(gca,'FontSize',16);
            if j == 1
                ylims = ylim;
                ylim([0 ylims(2)]);
            elseif j == 2
                set(gca,'yscale','log');
            end
        end
        
        %print plot        
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 8]);
        print([folder_Plots,'BSNE_sample_',FluxInterval.Site,'_',datestr(FluxInterval.Date),'_',int2str(ind_Interval),'.png'],'-dpng');
    end
end