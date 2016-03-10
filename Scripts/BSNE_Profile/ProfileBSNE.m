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
sigma_Q_all = cell(N_Sites,1); %flux uncertainty list
zq_all = cell(N_Sites,1); %flux height list
sigma_zq_all = cell(N_Sites,1); %flux height uncertainty list
Chi2_Qfit_all = cell(N_Sites,1); %normalized Chi2 for flux profile fit
df_Qfit_all = cell(N_Sites,1); %degrees of freedom for flux profile fit

%% Go through each site
for i = 1:N_Sites

    %% load BSNE data
    load(strcat(folder_DataBSNE,'FluxBSNE_',Sites{i})); 
    
    %% get number of profiles
    N_profile = length(FluxBSNE);
    
    %% initialize list of fluxes
    zq_all{i} = zeros(N_profile,1); %flux height list
    sigma_zq_all{i} = zeros(N_profile,1); %uncertainty in flux height list
    Q_all{i} = zeros(N_profile,1); %flux list
    sigma_Q_all{i} = zeros(N_profile,1); %uncertainty in flux list
    
    %% initialize list of Chi2
    if strcmp(Sites{i},'Oceano') %%for Oceano, two profiles
        Q_Oceano = cell(2,1);
        Q_1_all = zeros(N_profile,1); %flux list;
        Q_2_all = zeros(N_profile,1); %flux list;
        Chi2_Qfit_all{i} = cell(2,1);
        Chi2_Qfit_all{i}{1} = zeros(N_profile,1); %Chi2 for flux profile fit
        Chi2_Qfit_all{i}{2} = zeros(N_profile,1); %Chi2 for flux profile fit
        df_Qfit_all{i} = cell(2,1);
        df_Qfit_all{i}{1} = zeros(N_profile,1); %degrees of freedom for flux profile fit
        df_Qfit_all{i}{2} = zeros(N_profile,1); %degrees of freedom for flux profile fit
    else
        Chi2_Qfit_all{i} = zeros(N_profile,1); %Chi2 for flux profile fit
        df_Qfit_all{i} = zeros(N_profile,1); %degrees of freedom for flux profile fit
    end
        
    %% go through each interval
    for j = 1:N_profile

        %extract general profile info
        FluxInterval = FluxBSNE(j);
        zq_all{i}(j) = FluxInterval.z.zq; %flux height for time interval
        sigma_zq_all{i}(j) = FluxInterval.z.sigma_zq; %flux height uncertainty
        Q_all{i}(j) = FluxInterval.Q.Q; %flux for time interval
        sigma_Q_all{i}(j) = FluxInterval.Q.sigma_Q; %flux uncertainty
        
        %%for Oceano, two profiles
        if strcmp(Sites{i},'Oceano')
            % Get data for profile
            z_1 = FluxInterval.z.z_1;
            z_2 = FluxInterval.z.z_2;
            q_1 = FluxInterval.qz.qz_1;
            q_2 = FluxInterval.qz.qz_2;
            sigma_q_1 = FluxInterval.qz.sigma_1;
            sigma_q_2 = FluxInterval.qz.sigma_2;
            q_pred_1 = FluxInterval.qz.qz_pred_1;
            q_pred_2 = FluxInterval.qz.qz_pred_2;
            [z_sort_1, ind_sort_1] = sort(z_1);
            [z_sort_2, ind_sort_2] = sort(z_2);
            q_pred_sort_1 = q_pred_1(ind_sort_1);
            q_pred_sort_2 = q_pred_2(ind_sort_2);

            %get and record list of Chi2
            q_residuals_1 = q_pred_1 - q_1; %residuals between observed and predicted q           
            q_residuals_2 = q_pred_2 - q_2; %residuals between observed and predicted q           
            Chi2_Qfit_1 = sum((q_residuals_1./sigma_q_1).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
            Chi2_Qfit_2 = sum((q_residuals_2./sigma_q_2).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
            df_Qfit_1 = length(q_1)-2; %compute degrees of freedom for fitting
            df_Qfit_2 = length(q_2)-2; %compute degrees of freedom for fitting
            Q_1_all(j) = FluxInterval.Q.Q_1; %flux for time interval
            Q_2_all(j) = FluxInterval.Q.Q_2; %flux for time interval
            Chi2_Qfit_all{i}{1}(j) = Chi2_Qfit_1; %normalized Chi2 for flux profile fit
            Chi2_Qfit_all{i}{2}(j) = Chi2_Qfit_2; %normalized Chi2 for flux profile fit
            df_Qfit_all{i}{1}(j) = df_Qfit_1; %degrees of freedom for flux profile fit
            df_Qfit_all{i}{2}(j) = df_Qfit_2; %degrees of freedom for flux profile fit

            %make plot
            cla; hold on;
            errorbar(z_1,q_1,sigma_q_1,'b+','MarkerSize',10);
            errorbar(z_2,q_2,sigma_q_2,'ro','MarkerSize',10);
            plot(z_sort_1,q_pred_sort_1,'b',z_sort_2,q_pred_sort_2,'r');
            legend('profile 1','profile 2',...
                ['fit 1, \chi^2_{\nu} = ',num2str((Chi2_Qfit_1/df_Qfit_1),'%.2f')],...
                ['fit 2, \chi^2_{\nu} = ',num2str((Chi2_Qfit_2/df_Qfit_2),'%.2f')],...
                'Location','NorthEast');
            set(gca,'FontSize',PlotFont,'yscale','log');

            %print plot
            set(gca,'LooseInset',get(gca,'TightInset'));
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
            print([folder_Plots,'BSNE_FluxProfile_',Sites{i},'_',int2str(j),'.png'],'-dpng');
        
        %% for other sites, one profile
        else
            % Get data for profile
            z = FluxInterval.z.z;
            q = FluxInterval.qz.qz;
            sigma_q = FluxInterval.qz.sigma;
            q_pred = FluxInterval.qz.qz_pred;
            [z_sort, ind_sort] = sort(z);
            q_pred_sort = q_pred(ind_sort);

            %get and record list of Chi2
            q_residuals = q_pred - q; %residuals between observed and predicted q
            Chi2_Qfit = sum((q_residuals./sigma_q).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
            df_Qfit = length(q)-2; %compute degrees of freedom for fitting
            Chi2_Qfit_all{i}(j) = Chi2_Qfit; %normalized Chi2 for flux profile fit
            df_Qfit_all{i}(j) = df_Qfit; %degrees of freedom for flux profile fit

            %make plot 
            cla; hold on;
            errorbar(z,q,sigma_q,'b+','MarkerSize',10);
            plot(z_sort,q_pred_sort,'k');
            legend('data',['fit, \chi^2_{\nu} = ',num2str((Chi2_Qfit/df_Qfit),'%.2f')],'Location','NorthEast');
            set(gca,'FontSize',PlotFont,'yscale','log');

            %print plot
            set(gca,'LooseInset',get(gca,'TightInset'));
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
            print([folder_Plots,'BSNE_FluxProfile_',Sites{i},'_',int2str(j),'.png'],'-dpng');
        end
    end
end

%% PLOT RELATIVE UNCERTAINTY IN Q VERSUS TOTAL UNCERTAINTY
figure; hold on;
for i = 1:N_Sites
    plot(Q_all{i},sigma_Q_all{i}./Q_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field);
end
xlim([0 55]);
ylim([0 0.3]);
xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
ylabel('\textbf{BSNE relative uncertainty, $$\sigma_Q/Q$$}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'sigma_Q_BSNE_profilefit.png'],'-dpng');

%% PLOT RELATIVE UNCERTAINTY IN zq VERSUS TOTAL UNCERTAINTY
figure; hold on;
for i = 1:N_Sites
    plot(Q_all{i},sigma_zq_all{i}./zq_all{i},Markers_field{i},'Color',Colors_field{i},'MarkerSize',MarkerSize_field);
end
xlim([0 55]);
ylim([0 0.15]);
xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
ylabel('\textbf{BSNE relative uncertainty, $$\sigma_{zq}/{z_q}$$}','Interpreter','Latex');
set(gca,'FontSize',PlotFont);

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
print([folder_Plots,'sigma_zq_BSNE_profilefit.png'],'-dpng');

%% PLOT CHI2 FOR BSNE PROFILE VERSUS Q
figure;
for i = 1:N_Sites
    subplot(1,N_Sites,i); hold on;
    if strcmp(Sites{i},'Oceano')
        plot(Q_1_all,Chi2_Qfit_all{i}{1}./df_Qfit_all{i}{1},'b+','MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
        plot(Q_2_all,Chi2_Qfit_all{i}{2}./df_Qfit_all{i}{2},'ro','MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
        legend('profile 1','profile 2','Location','NorthEast');
    else
        plot(Q_all{i},Chi2_Qfit_all{i}./df_Qfit_all{i},'b+','MarkerSize',MarkerSize_field/2,'LineWidth',LineWidth_field);
    end
    xlim([0 55]);
    ylim([0.7e-1 2e2]);
    %ylim([1e-1 2e2]);
    set(gca,'Yscale','log','XMinorTick','On','YMinorTick','On','Box','On');
    xlabel('\textbf{Saltation flux, $$Q$$ (g m$$^{-2}$$ s$$^{-1}$$)}','Interpreter','Latex');
    if i==1
        ylabel('\textbf{Quality of fit for BSNE profile, $$\chi^2_{\nu}$$}','Interpreter','Latex');
    end
    title(Sites{i});
    set(gca,'FontSize',PlotFont);
end

%print plot
set(gca, 'LooseInset', get(gca,'TightInset'));
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 5]);
print([folder_Plots,'Chi2_BSNE_profilefit.png'],'-dpng');