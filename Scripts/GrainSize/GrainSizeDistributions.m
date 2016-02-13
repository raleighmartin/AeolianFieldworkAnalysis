%% initialize
clearvars;

%% information about where to load data and save plots
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Plots = '../../PlotOutput/GrainSize/'; %folder for plots

%% Information about sites
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);
zq_Sites = [0.0960, 0.1079, 0.0550]; %Saltation height by site.  Use this to identify BSNE height for grain size distribution

%% initialize combined surface grain size distribution
GrainSizeSurface_Combined = cell(N_Sites,1);
GrainSizeBSNE_Combined = cell(N_Sites,1);

%% plotting information
PlotFont = 12;
LineWidth_Surface = 1;

%% go through each site
for i = 1:N_Sites
%for i = 1;
    %load grain size data
    GrainSizeData_Path = strcat(folder_ProcessedData,'GrainSize_',Sites{i});
    load(GrainSizeData_Path);
    
    %load processed data
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path);
        
    %get unique BSNE names
    DateBSNE = [GrainSize_BSNE.Date]; %BSNE dates
    StartTimeBSNE = [GrainSize_BSNE.StartTime]; %BSNE start times
    EndTimeBSNE = [GrainSize_BSNE.EndTime]; %BSNE end times
    TimeBSNE = mean([StartTimeBSNE;EndTimeBSNE]); %BSNE midpoint times
    NameBSNE_day = {GrainSize_BSNE.NameBSNE}; %BSNE names
    NameBSNE = unique(NameBSNE_day); %unique BSNE names
    
    N_Samples = length(GrainSize_Surface);
    for j = 1:N_Samples
        %get surface size distribution
        d_surface = [GrainSize_Surface(j).gsd(2:end-1).Sizeclass_mid_mm];
        dlogd_surface = log([GrainSize_Surface(j).gsd(2:end-1).Sizeclass_upper_mm]) - log([GrainSize_Surface(j).gsd(2:end-1).Sizeclass_lower_mm]);
        dV_surface = [GrainSize_Surface(j).gsd(2:end-1).retained];
        dV_dlogd_surface = dV_surface./dlogd_surface;
        
        %get BSNE size distribution - average over all of day's samples for each BSNE
        indBSNE_day = find(DateBSNE==GrainSize_Surface(j).Date); %find BSNE samples from same day
        NameBSNE_day = {GrainSize_BSNE(indBSNE_day).NameBSNE}; %BSNE names
        NameBSNE = unique(NameBSNE_day); %unique BSNE names
        N_BSNE = length(NameBSNE);
        
        %get height from StressFluxWindows for nearest BSNE height
        
        %initialize grain sizes for BSNE
        d_BSNE = cell(N_BSNE,1);
        dV_dlogd_BSNE = cell(N_BSNE,1);
        
        for k = 1:N_BSNE
            indBSNE = indBSNE_day(strcmp(NameBSNE{k},NameBSNE_day));
            N_ind = length(indBSNE); %how many samples?
            d = cell(N_ind,1); %initialize grain size
            dlogd = cell(N_ind,1); %initialize dlogd
            dV = cell(N_ind,1); %initialize pdf
            for m = 1:N_ind
                d{m} = [GrainSize_BSNE(indBSNE(m)).gsd(2:end-1).Sizeclass_mid_mm];
                dlogd{m} = log([GrainSize_BSNE(indBSNE(m)).gsd(2:end-1).Sizeclass_upper_mm]) - log([GrainSize_BSNE(indBSNE(m)).gsd(2:end-1).Sizeclass_lower_mm]);
                dV{m} = [GrainSize_BSNE(indBSNE(m)).gsd(2:end-1).retained];
            end
            d_BSNE{k} = d{1};
            if N_ind==1
                dV_dlogd_BSNE{k} = dV{1}./dlogd{1};
            else
                dV_dlogd_BSNE{k} = mean(cell2mat(dV))./dlogd{1};
            end
        end
                    
        %plot
        figure(1); clf; hold on;
        plot(d_surface,dV_dlogd_surface,'LineWidth',LineWidth_Surface);
        for k=1:N_BSNE
            plot(d_BSNE{k},dV_dlogd_BSNE{k});
        end
        set(gca,'xscale','log','yscale','log');
        legend('Surface',NameBSNE{:},'Location','EastOutside');
        ylim([5 2e2]);
        xlim([5e-2 2]);
        ax = gca;
        ax.XTick = [0.05:0.01:0.1, 0.2:0.1:1, 2];
        ax.XTickLabel = {'0.05','','','','','0.1','0.2','0.3','','0.5','','','','','1','2'};
        ax.YTick = [5:10, 20:10:100, 200];
        ax.YTickLabel = {'5','','','','','10','20','30','','50','','','','','100','200'};
        
        %label plot
        xlabel('Particle size, d (mm)');
        ylabel('Norm. size distr., dV/dln(d)');
        title([Sites{i},', ',GrainSize_Surface(j).Location,', ',datestr(GrainSize_Surface(j).CollectionTime)]);
        set(gca,'XMinorTick','On','YMinorTick','On','Box','On');
        set(gca,'FontSize',PlotFont);
        
        %print plot
        PrintName = [Sites{i},'_',GrainSize_Surface(j).Location,'_',datestr(GrainSize_Surface(j).CollectionTime,'yyyymmdd_HHMM')];
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 3]);
        print([folder_Plots,PrintName,'.png'],'-dpng');
%         set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.4 3]);
%         print([folder_Plots,PrintName,'.eps'],'-depsc');
    end
end