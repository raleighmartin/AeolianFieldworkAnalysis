%% initialize
clearvars;
close all;

%%
%%%%%%%%%%%%%%%%
% DATA SOURCES %
%%%%%%%%%%%%%%%%

%% information about where to load data and save plots
folder_AnalysisData = '../../AnalysisData/SizeSelective/'; %folder for loading/saving analysis data
folder_LitData = '../../AnalysisData/Literature/'; %folder for loading/saving literature data
AnalysisData_Path = strcat(folder_AnalysisData,'SizeSelectiveAnalysis');
LitData_Path = strcat(folder_LitData,'LitAnalysis');
folder_Plots = '../../PlotOutput/SizeSelective/'; %folder for plots
folder_Functions = '../Functions/'; %folder with functions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA LOADING AND EXTRACTION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load data and functions
load(AnalysisData_Path); %load analysis data
load(LitData_Path); %load literature data
addpath(folder_Functions); %point MATLAB to location of functions

%%
%%%%%%%%%%%%%%%%%
% PLOTTING INFO %
%%%%%%%%%%%%%%%%%

%% plotting information
PlotFont = 10; %font for labels
LineWidth_Plot = 1; %width of lines
Marker_bin = {'o','s','d','^','v','>','<'}; %markers for bins
Color_bin = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]}; %colors for bins

%%
%%%%%%%%%
% PLOTS %
%%%%%%%%%

%% plot tau/tau_it-conditioned d_ref,air/d_50,sfc versus z/z_q
% Fig 1: d_ref = d_50
% Fig 2: d_ref = d_10
% Fig 3: d_ref = d_90

for n = 1:3
    if n == 1
        drefnorm_profile_airborne_Namikas06 = d50norm_profile_airborne_Namikas06;
        drefnorm_profile_airborne_taunorm_Farrell12 = d50norm_profile_airborne_taunorm_Farrell12;
        drefnorm_profile_airborne_taunorm_Cluster = d50norm_profile_airborne_taunorm_Cluster;
    elseif n == 2
        drefnorm_profile_airborne_Namikas06 = d10norm_profile_airborne_Namikas06;
        drefnorm_profile_airborne_taunorm_Farrell12 = d10norm_profile_airborne_taunorm_Farrell12;
        drefnorm_profile_airborne_taunorm_Cluster = d10norm_profile_airborne_taunorm_Cluster;
    elseif n == 3
        drefnorm_profile_airborne_Namikas06 = d90norm_profile_airborne_Namikas06;
        drefnorm_profile_airborne_taunorm_Farrell12 = d90norm_profile_airborne_taunorm_Farrell12;
        drefnorm_profile_airborne_taunorm_Cluster = d90norm_profile_airborne_taunorm_Cluster;
    end
    
    %initialize figure and subplots
    figure(n); clf;
    h_subplot = gobjects(8,1);

    % Namikas
    h_subplot(1) = subplot('position',[0.06 0.73 0.275 0.23]); hold on;

    %plot airborne profiles
    for i = 1:N_taunorm_bins
        %get index of Namakis06 associated with taunorm_bin 
        ind_Namikas06 = intersect(find(taunorm_Namikas06>taunorm_min_bins(i)),...
            find(taunorm_Namikas06<taunorm_max_bins(i)));
        if ~isempty(ind_Namikas06) %plot only if there is a profile
            plot(drefnorm_profile_airborne_Namikas06{ind_Namikas06},znorm_Namikas06,...     
                [Marker_bin{i},'-']);
        else %for those profiles with no data, create a dummy plot so that legend renders properly
            plot([0 0],[0 0],[Marker_bin{i},'-']);
        end
    end

    %format plot
    htitle = title('Namikas, 2006');
    set(htitle,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','YScale','log','Box','On');
    set(gca,'YTickLabel',{'0.1','1','10'})
    ylabel('Non-dim. ht., $$z/<z_q>$$','Interpreter','Latex')
    ylim([0.1 10]);
    if n == 1
        xlim([0.4 1.0]);   
        text(0.42,7.5,'(a)'); %label panel
    elseif n == 2
        xlim([0.2 0.8]);
        text(0.22,7.5,'(a)'); %label panel          
    elseif n == 3
        xlim([0.6 1.2]);
        text(0.62,7.5,'(a)'); %label panel
    end
    set(gca,'FontSize',PlotFont);

    %% Farrell
    h_subplot(2) = subplot('position',[0.385 0.73 0.275 0.23]); hold on;

    %plot airborne profiles
    for i = 1:N_taunorm_bins
        ind_plot = find(~isnan(drefnorm_profile_airborne_taunorm_Farrell12{i}));
        if ~isempty(ind_plot) %plot only profiles with data
            plot(drefnorm_profile_airborne_taunorm_Farrell12{i},...
                znorm_taunorm_Farrell12{i},[Marker_bin{i},'-']);
        else %for those profiles with no data, create a dummy plot so that legend renders properly
            plot([0 0],[0 0],[Marker_bin{i},'-']); 
        end
    end

    %format plot
    htitle = title('Farrell et al., 2012');
    set(htitle,'FontSize',PlotFont);
    set(gca,'XMinorTick','On','YMinorTick','On','YScale','log','Box','On');
    set(gca,'YTickLabel',{'0.1','1','10'})
    ylim([0.1 10]);
    if n == 1
        xlim([0.4 1.0]);   
        text(0.42,7.5,'(b)'); %label panel
    elseif n == 2
        xlim([0.2 0.8]);
        text(0.22,7.5,'(b)'); %label panel          
    elseif n == 3
        xlim([0.6 1.2]);
        text(0.62,7.5,'(b)'); %label panel
    end
    set(gca,'FontSize',PlotFont);

    %% our data
    panel_labels = {'(c)','(d)','(e)','(f)','(g)','(h)'};
    for i = 1:N_Cluster
        if i == 1
            h_subplot(3) = subplot('position',[0.06 0.41 0.275 0.23]); hold on;
        elseif i == 2
            h_subplot(4) = subplot('position',[0.385 0.41 0.275 0.23]); hold on;
        elseif i == 3 %make this subplot bigger to accommodate legend
            h_subplot(5) = subplot('position',[0.71 0.41 0.275 0.455]); hold on;
        elseif i == 4
            h_subplot(6) = subplot('position',[0.06 0.09 0.275 0.23]); hold on;
        elseif i == 5
            h_subplot(7) = subplot('position',[0.385 0.09 0.275 0.23]); hold on;
        else
            h_subplot(8) = subplot('position',[0.71 0.09 0.275 0.23]); hold on;
        end

        %plot airborne profiles - d50air/d50sfc
        for j = 1:N_taunorm_bins
            ind_plot = find(~isnan(drefnorm_profile_airborne_taunorm_Cluster{i}(j,:)));
            if ~isempty(ind_plot) %plot only profiles with data
                h_plot = plot(drefnorm_profile_airborne_taunorm_Cluster{i}(j,ind_plot),...
                    znorm_profile_airborne_taunorm_Cluster{i}(j,ind_plot),...
                    [Marker_bin{j},'-']);
                get(h_plot,'Color')
            else %for those profiles with no data, create a dummy plot so that legend renders properly
                plot([0 0],[0 0],[Marker_bin{j},'-']);
            end
        end

        %format plot
        htitle = title(ClusterNames{i});
        set(htitle,'FontSize',PlotFont);
        if mod(i,3) == 1
            ylabel('Non-dim. ht., $$z/<z_q>$$','Interpreter','Latex')
        end
        if i>=2*round(N_Cluster/2)-2
            if n == 1
                xlabel('Non-dim. airborne g.s., $${d}_{50,air}/{d}_{50,bed}$$','Interpreter','Latex')
            elseif n == 2
                xlabel('Non-dim. airborne g.s., $${d}_{10,air}/{d}_{50,bed}$$','Interpreter','Latex')
            elseif n == 3
                xlabel('Non-dim. airborne g.s., $${d}_{90,air}/{d}_{50,bed}$$','Interpreter','Latex')
            end
        end
        ylim([0.1 10]);
        if n == 1
            xlim([0.4 1.0]);   
            text(0.42,7.5,panel_labels{i}); %label panel
        elseif n == 2
            xlim([0.2 0.8]);
            text(0.22,7.5,panel_labels{i}); %label panel          
        elseif n == 3
            xlim([0.6 1.2]);
            text(0.62,7.5,panel_labels{i}); %label panel
        end
        set(gca,'XMinorTick','On','YMinorTick','On','YScale','log','Box','On');
        set(gca,'YTickLabel',{'0.1','1','10'})
        set(gca,'FontSize',PlotFont);

        %create legend
        legend_items = cell(N_taunorm_bins,1);
        for j = 1:N_taunorm_bins
            if j == 1
                legend_items{j} = ['\tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];       
            else
                legend_items{j} = [num2str(taunorm_min_bins(j),'%10.1f'),' < \tau/\tau_{it} \leq ',num2str(taunorm_max_bins(j),'%10.1f')];
            end
        end
        if i == 3
            h_legend = legend(legend_items,'Location','NorthOutside');
            set(h_legend,'FontSize',8);
        end
    end

    %print plot
    set(gcf,'PaperUnits','inches','PaperSize',[8 6],'PaperPosition',[0 0 8 6],'PaperPositionMode','Manual');
    if n == 1
        print([folder_Plots,'d50norm_profile_airborne_comparison.png'],'-dpng');
    elseif n == 2
        print([folder_Plots,'d10norm_profile_airborne_comparison.png'],'-dpng');
    elseif n == 3
        print([folder_Plots,'d90norm_profile_airborne_comparison.png'],'-dpng');
    end
end

%% plot q_norm vs z_norm and exponential profile for each site
% znorm = z / zq for Cluster
% qnorm = q / q0 for profile

ind_d_all = [6,13,19]; %indices of d's to plot

%create three separate plots
for n = 1:3
    if n == 1
        ind_Cluster_all = [1,2];
    elseif n == 2
        ind_Cluster_all = [3,4];
    elseif n == 3
        ind_Cluster_all = [5,6];
    end

    %initialize plot
    figure(n+3); clf;

    %go through first two sites
    for i = 1:2

        %go through three different diameters
        for j = 1:3

            %initialize subplot
            if i == 1
                if j == 1
                    h_subplot(1) = subplot('position',[0.06 0.57 0.275 0.38]);
                elseif j == 2
                    h_subplot(3) = subplot('position',[0.385 0.57 0.275 0.38]);
                elseif j == 3
                    h_subplot(3) = subplot('position',[0.71 0.57 0.275 0.38]);
                end
            elseif i == 2
                if j == 1
                    h_subplot(4) = subplot('position',[0.06 0.09 0.275 0.38]);
                elseif j == 2
                    h_subplot(5) = subplot('position',[0.385 0.09 0.275 0.38]);
                elseif j == 3
                    h_subplot(6) = subplot('position',[0.71 0.09 0.275 0.38]);
                end
            end
            hold on;

            %information on profiles
            ind_Cluster = ind_Cluster_all(i); %which cluster to plot?
            ind_d = ind_d_all(j); %which grain-size bin to plot?

            %get Non-dimensionalized profiles
            N_profile = length(qi_Cluster{ind_Cluster}); %get number of profiles
            znorm_plot = cell(N_profile,1);
            qnorm_plot = cell(N_profile,1);
            L_profile = zeros(N_profile,1);
            for k = 1:N_profile
                [~,ind_sort] = sort(z_profile_Cluster{ind_Cluster}{k}); %sort by ascending z
                znorm_plot{k} = z_profile_Cluster{ind_Cluster}{k}(ind_sort)./zqi_bar_Cluster{ind_Cluster}(ind_d); %get sorted znorms for plotting
                qnorm_plot{k} = (qi_Cluster{ind_Cluster}{k}(ind_sort,ind_d)./q0i_Cluster{ind_Cluster}(k,ind_d)); %get sorted qnorms for plotting
                L_profile(k) = length(znorm_plot{k});
            end
                
            %get reference profile
            znorm_min = min(vertcat(znorm_plot{:})); %get min znorm
            znorm_max = max(vertcat(znorm_plot{:})); %get max znorm

            %plot profiles
            plot(exp(-[znorm_min znorm_max]),[znorm_min znorm_max],'k','LineWidth',2); %plot reference profile
            for k = 1:N_profile %plot individual profiles
                if L_profile(k)==max(L_profile) %plot only full profiles
                    plot(qnorm_plot{k},znorm_plot{k});
                end
            end

            %format profile
            if i == 2
                xlabel('Non-dim. size-selective flux, $$q_{i}(z)/q_{i,0}$$','Interpreter','Latex');
            end
            if j == 1
                ylabel('Non-dimensionalized height, $$z/\langle z_{q,i} \rangle$$','Interpreter','Latex');
            end
            htitle = title([ClusterNames{ind_Cluster},', d = ',num2str(d_f_mid_Cluster{ind_Cluster}(ind_d),3),' mm']);
            set(htitle,'FontSize',PlotFont-2);
            set(gca,'YMinorTick','On','XMinorTick','On','XScale','log','Box','On');
            set(gca,'FontSize',PlotFont);
            
            %get and update plotting limits
            xlims_init = xlim;
            xlims_new = [10^(floor(log10(xlims_init(1)))) 1];
            xlim(xlims_new);
            ylims_init = ylim;
            ylims_new = [0 ceil(ylims_init(2))];
            ylim(ylims_new);
            
            %correct some labels
            if xlims_new(1) <= 1e-5
                set(gca,'XMinorTick','Off');
                set(gca,'XTick',[10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0)]);
                set(gca,'XTickLabel',{'10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
            elseif xlims_new(1) <= 1e-4
                set(gca,'XMinorTick','Off');
                set(gca,'XTick',[10^(-4),10^(-3),10^(-2),10^(-1),10^(0)]);
                set(gca,'XTickLabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
            end
                
            %panel labels
            x_text = 10^(0.9*(log10(xlims_new(2))-log10(xlims_new(1)))+log10(xlims_new(1)));
            y_text = 0.9*ylims_new(2);
            if i == 1
                if j == 1
                    text(x_text,y_text,'(a)')
                elseif j == 2
                    text(x_text,y_text,'(b)')
                elseif j == 3
                    text(x_text,y_text,'(c)')
                end
            elseif i == 2
                if j == 1
                    text(x_text,y_text,'(d)')
                elseif j == 2
                    text(x_text,y_text,'(e)')
                elseif j == 3
                    text(x_text,y_text,'(f)')
                end
            end
        end
    end

    %print plot
    set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
    print([folder_Plots,'profile_demo_',int2str(n),'.png'],'-dpng');
end

%% plot full profiles at each site

%initialize plot
figure(7); clf;

%go through each site
for i = 1:6
    if i == 1
        h_subplot(1) = subplot('position',[0.06 0.57 0.275 0.38]);
    elseif i == 2
        h_subplot(3) = subplot('position',[0.385 0.57 0.275 0.38]);
    elseif i == 3
        h_subplot(3) = subplot('position',[0.71 0.57 0.275 0.38]);
    elseif i == 4
        h_subplot(4) = subplot('position',[0.06 0.09 0.275 0.38]);
    elseif i == 5
        h_subplot(5) = subplot('position',[0.385 0.09 0.275 0.38]);
    elseif i == 6
        h_subplot(6) = subplot('position',[0.71 0.09 0.275 0.38]);
    end
    hold on;

    %get Non-dimensionalized profiles
    N_profile = length(q_profile_Cluster{i}); %get number of profiles
    znorm_plot = cell(N_profile,1);
    qnorm_plot = cell(N_profile,1);
    L_profile = zeros(N_profile,1);
    for k = 1:N_profile
        [~,ind_sort] = sort(z_profile_Cluster{i}{k}); %sort by ascending z
        znorm_plot{k} = z_profile_Cluster{i}{k}(ind_sort)./zq_bar_Cluster(i); %get sorted znorms for plotting
        qnorm_plot{k} = ((q_profile_Cluster{i}{k}(ind_sort).*zq_bar_Cluster(i))/Q_profile_Cluster{i}(k)); %get sorted qnorms for plotting
        L_profile(k) = length(znorm_plot{k});
    end

    %get reference profile
    znorm_min = min(vertcat(znorm_plot{:})); %get min znorm
    znorm_max = max(vertcat(znorm_plot{:})); %get max znorm

    %plot profiles
    plot(exp(-[znorm_min znorm_max]),[znorm_min znorm_max],'k','LineWidth',2); %plot reference profile
    for k = 1:N_profile %plot individual profiles
        if L_profile(k)==max(L_profile) %plot only full profiles
            plot(qnorm_plot{k},znorm_plot{k});
        end
    end
    
    %format profile
    if i >= 4 
        xlabel('Non-dim. saltation flux, $$q(z)/q_{0}$$','Interpreter','Latex');
    end
    if i == 1 || i == 4
        ylabel('Non-dim. height, $$z/\langle z_{q} \rangle$$','Interpreter','Latex');
    end
    htitle = title(ClusterNames{i});
    set(htitle,'FontSize',PlotFont);
    set(gca,'YMinorTick','On','XMinorTick','On','XScale','log','Box','On');
    set(gca,'FontSize',PlotFont);

    %get and update plotting limits
    xlims_init = xlim;
    xlims_new = [10^(floor(log10(xlims_init(1)))) 1];
    xlim(xlims_new);
    ylims_init = ylim;
    ylims_new = [0 ceil(ylims_init(2))];
    ylim(ylims_new);

    %correct some labels
    if xlims_new(1) <= 1e-5
        set(gca,'XMinorTick','Off');
        set(gca,'XTick',[10^(-5),10^(-4),10^(-3),10^(-2),10^(-1),10^(0)]);
        set(gca,'XTickLabel',{'10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
    elseif xlims_new(1) <= 1e-4
        set(gca,'XMinorTick','Off');
        set(gca,'XTick',[10^(-4),10^(-3),10^(-2),10^(-1),10^(0)]);
        set(gca,'XTickLabel',{'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
    end

    %panel labels
    x_text = 10^(0.9*(log10(xlims_new(2))-log10(xlims_new(1)))+log10(xlims_new(1)));
    y_text = 0.9*ylims_new(2);
    %x_text = 0.5012;
    %y_text = 5.4000;
    if i == 1
        text(x_text,y_text,'(a)')
    elseif i == 2
        text(x_text,y_text,'(b)')
    elseif i == 3
        text(x_text,y_text,'(c)')
    elseif i == 4
        text(x_text,y_text,'(d)')
    elseif i == 5
        text(x_text,y_text,'(e)')
    elseif i == 6
        text(x_text,y_text,'(f)')
    end

    %print plot
    set(gcf,'PaperUnits','inches','PaperSize',[8 5],'PaperPosition',[0 0 8 5],'PaperPositionMode','Manual');
    print([folder_Plots,'Cluster_profiles.png'],'-dpng');
end