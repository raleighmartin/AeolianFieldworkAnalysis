%% initialize
clearvars;

%% information about where to load data and save plots
folder_DataBSNE = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder containing BSNE data
folder_Plots = '../../PlotOutput/BSNE_Profile/'; %folder for plots

%% information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% information for plotting
Markers_field = {'s','d','o'};
Colors_field = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.2116 0.1898 0.5777],[0.6473 0.7456 0.4188]};
MarkerSize_field = 8;
LineWidth_field = 1;
PlotFont = 14;

%% initialize plot
close all;
h = figure;
set(h,'visible','off');
set(gca,'FontSize',14,'XMinorTick','On','YMinorTick','On','Box','On');
xlabel('$$z$$ (m)','Interpreter','Latex');
ylabel('$$q$$ (g/m$$^2$$/s)','Interpreter','Latex');
set(gcf, 'PaperPosition',[0 0 8 6]);

%% initialize list of Chi2
Q_all = cell(N_Sites,1); %flux list
Chi2_Qfit_all = cell(N_Sites,1); %normalized Chi2 for flux profile fit
df_Qfit_all = cell(N_Sites,1); %degrees of freedom for flux profile fit

%% Go through each site
for i = 1:N_Sites

    %% load BSNE data
    load(strcat(folder_DataBSNE,'FluxBSNE_',Sites{i})); 
    
    %% get number of profiles
    N_profile = length(FluxBSNE);
    
    %% initialize list of Chi2
    Q_all{i} = zeros(N_profile,1); %flux list
    Chi2_Qfit_all{i} = zeros(N_profile,1); %normalized Chi2 for flux profile fit
    df_Qfit_all{i} = zeros(N_profile,1); %degrees of freedom for flux profile fit
    
    %% go through each interval
    for j = 1:N_profile

        %% Get data for profile
        FluxInterval = FluxBSNE(j);
        z = FluxInterval.z.z;
        q = FluxInterval.qz.qz;
        sigma_q = FluxInterval.qz.sigma;
        q_pred = FluxInterval.qz.qz_pred;
        %sigma_q_pred = FluxInterval.qz.sigma_qz_pred;
        %q_pred_minus = q_pred - sigma_q_pred;
        %q_pred_plus = q_pred + sigma_q_pred;
        [z_sort, ind_sort] = sort(z);
        q_pred_sort = q_pred(ind_sort);
        %q_pred_minus_sort = q_pred_minus(ind_sort);
        %q_pred_plus_sort = q_pred_plus(ind_sort);

        %get and record list of Chi2
        q_residuals = q_pred - q; %residuals between observed and predicted q
        Chi2_Qfit = sum((q_residuals./sigma_q).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
        df_Qfit = length(q)-2; %compute degrees of freedom for fitting
        Q_all{i}(j) = FluxInterval.Q.Q; %flux for time interval
        Chi2_Qfit_all{i}(j) = Chi2_Qfit; %normalized Chi2 for flux profile fit
        df_Qfit_all{i}(j) = df_Qfit; %degrees of freedom for flux profile fit
                
        %make plot 
        cla; hold on;
        errorbar(z,q,sigma_q,'b+','MarkerSize',10);
        plot(z_sort,q_pred_sort,'k');
        legend('data',['fit, \chi^2_{\nu} = ',num2str((Chi2_Qfit/df_Qfit),'%.2f')],'Location','NorthEast');
        set(gca,'FontSize',PlotFont);
        
        %print plot
        set(gca,'LooseInset',get(gca,'TightInset'));
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
        print([folder_Plots,'BSNE_FluxProfile_',Sites{i},'_',int2str(j),'.png'],'-dpng');
    end
end

%% PLOT CHI2 FOR BSNE PROFILE VERSUS Q
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    plot(Q_all{i},Chi2_Qfit_all{i}./df_Qfit_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    xlim([0 55]);
    ylim([1e-1 2e2]);
    set(gca,'Yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Quality of fit for BSNE profile, $$\chi^2_{\nu}$$}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot for draft
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'Chi2_BSNE_profilefit.png'],'-dpng');