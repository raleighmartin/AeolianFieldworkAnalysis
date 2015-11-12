% %% SCRIPT TO GENERATE WINDOWS OF SALTATION FLUX AND STRESS VALUES FOR ANALYSIS
% % SCRIPT DEPENDENCIES:
% % CreateTimeBlocks.m
% % ExtractVariableTimeInterval.m
% % IntersectingTimeIntervals.m
% % qz_profilefit
% % window_average
% 
% %% initialize
% clearvars;
% 
%% information about sites for analysis
Sites = {'Jericoacoara';'RanchoGuadalupe';'Oceano'};
N_Sites = length(Sites);

%% set time interval for computing velocity profiles
ProfileTimeInterval = duration(0,30,0); %duration of window for computations
RunningTimeInterval = duration(0,5,0); %offset to use for running averages, to enrich dataset
N_RunningPerProfile = floor(ProfileTimeInterval/RunningTimeInterval)-1; %number of offsets for starting times

%% set durations of window average for flux frequencies
dt_min_s = 0.04; %mininum window average time
dt_max_s = 600; %maximum window average time
N_dt = 30; %number of durations of window average
dt_s = unique(round(logspace(0,log10(dt_max_s/dt_min_s),N_dt))*dt_min_s); %create window average dt's
N_dt = length(dt_s); %recalculate number of window average dt's after removing repeats
dt_windowaverage = duration(0,0,dt_s);

%% information about where to load data and save plots
folder_ProcessedData = '../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_AnalysisData = '../AnalysisData/'; %folder for storing outputs of this analysis
folder_Plots = '../PlotOutput/IntermittencyFactor/'; %folder for plots

SaveData_Path = strcat(folder_AnalysisData,'FluxFrequencyWindows_all');

% %% load processed and metadata for each site, add to structured arrays of all data and metadata
% Data = cell(N_Sites,1);
% for i = 1:N_Sites
%     ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
%     Data{i} = load(ProcessedData_Path); %load processed data
% end

%% initialize flux variable lists
Q_all = cell(N_Sites,1); %total flux
fD_all = cell(N_Sites,1); %Wenglor total transport detection frequency
fQ_all = cell(N_Sites,1); %Estimated actual transport frequency
fQalt_all = cell(N_Sites,1); %Estimated actual transport frequency (alternative calc)

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
    
    %% get flux data from overall processed data file
    FluxData = Data{i}.ProcessedData.FluxWenglor; %Wenglor flux data
    
    %% get start times and end times for flux observations
    FluxStartTimes = [FluxData.StartTime]';
    FluxEndTimes = [FluxData.EndTime]';
    
    %create time blocks based on flux start and end times
    [BlockStartTimes, BlockEndTimes] = ...
            CreateTimeBlocks(FluxStartTimes, FluxEndTimes, ProfileTimeInterval);
    
    %Add in additional blocks for running average offsets
    for j = 1:N_RunningPerProfile
        [BlockStartTimes_j, BlockEndTimes_j] = ...
            CreateTimeBlocks(FluxStartTimes+(j*RunningTimeInterval), FluxEndTimes, ProfileTimeInterval);
        BlockStartTimes = [BlockStartTimes; BlockStartTimes_j];
        BlockEndTimes = [BlockEndTimes; BlockEndTimes_j];
    end
    BlockStartTimes = sort(BlockStartTimes);
    BlockEndTimes = sort(BlockEndTimes);
    N_Blocks = length(BlockStartTimes);
    
    %% initialize lists of flux values
    Q_all{i} = zeros(N_Blocks,1)*NaN; %total flux
    fD_all{i} = zeros(N_Blocks,N_dt)*NaN; %Wenglor detection frequency
    fQ_all{i} = zeros(N_Blocks,N_dt)*NaN; %Wenglor estimated transport frequency
    fQalt_all{i} = zeros(N_Blocks,N_dt)*NaN; %Wenglor estimated transport frequency - alt calc
    
    %% go through time blocks
    for j = 1:N_Blocks

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Blocks),', ',datestr(now)]
        
        %get specific start and end time
        StartTime = BlockStartTimes(j);
        EndTime = BlockEndTimes(j);
        
        %% FLUX CALCULATIONS FOR INTERVAL
        
        %extract time interval - get interval number and indices within interval for analysis
        [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(FluxData,StartTime,EndTime,'t','t','t');
        
        %determine if values span entire analysis interval
        if ~isempty(IntervalN)
            %use only first interval
            IntervalN = IntervalN(1);
            IntervalInd = IntervalInd{1};
            StartTime_Extraction = FluxData(IntervalN).t.t(IntervalInd(1));
            EndTime_Extraction = FluxData(IntervalN).t.t(IntervalInd(end));
            t_Interval = FluxData(IntervalN).t.t(IntervalInd);
        else %generate arbitrary values if there are no data
            StartTime_Extraction = datetime(0,0,0);
            EndTime_Extraction = datetime(0,0,0);
        end
        
        %restrict to time intervals with data spanning entire interval
        if (StartTime_Extraction==StartTime)&&(EndTime_Extraction==EndTime)
                    
            %get raw q and qcal values
            q = FluxData(IntervalN).qz.qz(IntervalInd,:);
            qcal = FluxData(IntervalN).qz.qzPerCount(IntervalInd,:);

            %reverse calculate particle arrival numbers per Wenglor per time interval
            n = q./qcal;

            %get total number of detected particles during whole time window
            N = sum(sum(n));

            %compute mean q profile
            q_profile = mean(q);
            z_profile = FluxData(IntervalN).z.z;
            N_z = length(z_profile);

            %Perform profile fit to get q0, ze, and Q
            [q0,ze] = qz_profilefit(q_profile,z_profile,0,0); %For now, ignore uncertainty in profile values
            Q = q0*ze; %get total flux [g/m/s]

            %convert to 0 if NaN and expected Q<0.05 g/m/s based on lowest Wenglor and ze = 10 cm
            if isnan(Q)&&((0.1*q_profile(1))/exp(-z_profile(1)/0.1)<0.05);
                Q=0;
            end

            %add to list
            Q_all{i}(j) = Q; %total flux

            %Determine flux frequencies for different window averaging times
            if Q==0 %if no flux, set frequencies to zero
                fD_all{i}(j,:) = 0;
                fQ_all{i}(j,:) = 0;
                fQalt_all{i}(j,:) = 0;
            else %otherwise, do calculations
                qsum = sum(q')'; %create vector qsum, which is sum of all q's for each time increment
                for k = 1:N_dt; %go through all window averaging times
                    qsum_avg = window_average(qsum, t_Interval, dt_windowaverage(k)); %compute window average
                    T = length(qsum_avg); %number of timesteps

                    %get frequency of detection
                    fD = (sum(qsum_avg>0))/T;

                    %estimate fQ iteratively 
                    fQ_prev = 0; %previous estimate of fQ (arbitrarily set to 0 for while loop functioning) 
                    fQ = fD; %start by assuming fQ=fD
                    while abs(fQ-fQ_prev)>0.001 %compare current and previous estimations of fQ for while loop
                        fQ_prev = fQ; %update fQ_prev with last fQ
                        fQ = fD/(1-exp(-N/(fQ_prev*T))); %calculate new fQ
                    end

                    %alternate estimate of fQ
                    fQalt = fD/(1-exp(-N/T)); %calculate new fQ

                    %remove points with fQ>1
                    if fQ>1
                        fQ = NaN;
                    end
                    if fQalt>1
                        fQalt = NaN;
                    end 
                    
                    %add to lists of all fD and fQ
                    fD_all{i}(j,k) = fD;
                    fQ_all{i}(j,k) = fQ;
                    fQalt_all{i}(j,k) = fQalt;
                end
            end
        end
         
        %KEEP ONLY INTERVALS WHERE FLUX IS WELL DEFINED
        ind_good = find(~isnan(Q_all{i}));
        Q_all{i} = Q_all{i}(ind_good); %total flux
        fD_all{i} = fD_all{i}(ind_good,:); %detection frequency
        fQ_all{i} = fQ_all{i}(ind_good,:); %flux frequency
        fQalt_all{i} = fQalt_all{i}(ind_good,:); %alternate calculation of flux frequency
    end
end

%plot example of timescale dependence of fQ
figure(1); clf; hold on;
i = 3;
j = 79;
%j = 1105;
%j = 403;
plot(dt_s,fD_all{i}(j,:),'k','LineWidth',2)
plot(dt_s,fQalt_all{i}(j,:),'b--','LineWidth',2)
plot(dt_s,fQ_all{i}(j,:),'g-.','LineWidth',2)
set(gca,'xscale','log');
xlabel('dt (s)');
ylabel('f');
legend('f_{D}','f_{Q} (method 1)','f_{Q} (method 2)','Location','SouthEast');
title(['Q = ',num2str(round(Q_all{i}(j),3)),' g/m/s']);
set(gca,'FontSize',16);
print([folder_Plots,'fQ_methods_',Sites{i},'_',int2str(j),'.png'],'-dpng');

%get min dt for which fQ>90% and plot this
dt_90 = zeros(size(Q_all{3}))*NaN;
for j = 1:length(fQ_all{3})
    dt_min = dt_s(find(fQ_all{3}(j,:)>0.9,1));
    if ~isempty(dt_min)
        dt_90(j)=dt_min;
    end
end
figure(2); clf; hold on;
plot(Q_all{3},dt_90,'o');
ylim([min(dt_s) max(dt_s)]);
xlabel('Q (g/m/s)');
ylabel('dt_{90%} (s)');
set(gca,'FontSize',16);

%get binned mean value and plot this
Q_bin = 1:2:29;
dt_90_Q_avg = zeros(size(Q_bin));
for j = 1:length(Q_bin)
    ind_Q = find(round((Q_all{3}+1)/2)==(Q_bin(j)+1)/2);
    dt_90_Q = dt_90(ind_Q);
    dt_90_Q_avg(j) = mean(dt_90_Q(~isnan(dt_90_Q)));
end
plot(Q_bin,dt_90_Q_avg,'LineWidth',2);
set(gca,'yscale','log');
print([folder_Plots,'dt90_Q_',Sites{i},'.png'],'-dpng');

%% SAVE DATA
save(SaveData_Path,'Sites','dt_windowaverage','*all');