%% function to take info about BSNE trap weights ('WeightBSNE') and grain sizes ('GrainSize_BSNE')
% and aggregate this into mass and grain-size profiles
% Dependencies: IntervalsStartsWithin, qzCalc, qz_profilefit,
% z_profile_calc
% Used by: DataExtraction_Jericoacoara2014.m

function FluxBSNE = ProcessBSNEs(WeightBSNE,GrainSize_BSNE)

%% assumed error values
sigma_FluxSampling = 0.1; %flux sampling fractional error (Goossens et al., 2000)
sigma_HeightWidthBSNE_cm = 0.1; %error in height and width of BSNE dimensions
sigma_DurationBSNE = duration(0,0,30); %error in duration of BSNE sample

%% get information about flux calculation time intervals, based on A1
indices_A1 = find(strcmp(WeightBSNE.NameBSNE,'A1')); %indices of A1 BSNE records
N_Intervals = length(indices_A1); %count intervals based on A1, which is in all profiles
StartTimes_A1 = WeightBSNE.StartTime(indices_A1); %get start times associated with A1
EndTimes_A1 = WeightBSNE.EndTime(indices_A1); %get end times associated with A1
Site_A1 = WeightBSNE.Site(indices_A1); %get sites from A1
Date_A1 = WeightBSNE.Date(indices_A1); %get dates from A1

%% determine which weights are not in error and can be included in flux profiles
indices_nonerror = find(WeightBSNE.ErrorCode==0); % get nonerror indices

%% Initialize 'FluxBSNE' structured array with info from A1 BSNE
FluxBSNE = struct;

%% Compute BSNE flux and uncertainty for each interval
for i = 1:N_Intervals
    
    %look for all time intervals that start within A1 interval
    [StartTimesContained, EndTimesContained, indices_Intervals] = IntervalsStartsWithin(StartTimes_A1(i),EndTimes_A1(i),WeightBSNE.StartTime,WeightBSNE.EndTime);

    %add site and date info to structured array
    FluxBSNE(i).Site = Site_A1{i};
    FluxBSNE(i).Date = Date_A1(i);
    
    %generate minimal time intervals (i.e., include only portions that are mutually overlapping for all traps)
    FluxBSNE(i).StartTime = max(StartTimesContained);
    FluxBSNE(i).EndTime = min(EndTimesContained);
    
    %duration of each individual BSNE
    DurationBSNEs = EndTimesContained-StartTimesContained;
    
    %reduce intervals based on error code
    indices_Intervals = intersect(indices_Intervals,indices_nonerror); %indices to BSNE non-error intervals to be included in profile
    N_ProfileBSNEs = length(indices_Intervals); %how many BSNEs are there for profile?
    
    %create lists for flux and grain-size profiles
    name = WeightBSNE.NameBSNE(indices_Intervals); %get names of BSNEs for interval
    FluxBSNE(i).name = name; %add names to structured array
    d_StartTime_profile = WeightBSNE.StartTime_GrainSize(indices_Intervals);
    d_EndTime_profile = WeightBSNE.EndTime_GrainSize(indices_Intervals);
    
    %get associated grain-size (d10, d50, d90) values for profile and info about combined samples
    d_10_profile = zeros(N_ProfileBSNEs,1)*NaN;
    d_50_profile = zeros(N_ProfileBSNEs,1)*NaN;
    d_90_profile = zeros(N_ProfileBSNEs,1)*NaN;
    d_IsCombined_profile = zeros(N_ProfileBSNEs,1)*NaN;

    for j = 1:N_ProfileBSNEs
        %find index in array associated with grain size sample for BSNE
        ind_GrainSize = intersect(...
            intersect(find([GrainSize_BSNE(:).StartTime]==d_StartTime_profile(j)),...
            find([GrainSize_BSNE(:).EndTime]==d_EndTime_profile(j))),...
            find(strcmp({GrainSize_BSNE(:).NameBSNE},name(j))));
        
        %get d10, d50, d90, and sample information for this (could add other values later if desired)
        if ~isempty(ind_GrainSize) %add only if a value is found, otherwise default is NaN
            d_10_profile(j) = GrainSize_BSNE(ind_GrainSize).d_10_mm;
            d_50_profile(j) = GrainSize_BSNE(ind_GrainSize).d_50_mm;
            d_90_profile(j) = GrainSize_BSNE(ind_GrainSize).d_90_mm;
            d_IsCombined_profile(j) = strcmp(GrainSize_BSNE(ind_GrainSize).Notes,'combined');
        end
    end
    
    %Get BSNE bottom heights, z_bottom (m), uncertainties, sigma_z_bottom,
    %and trap heights, z_trapheight (m). Add these to structured array
    z_bottom = WeightBSNE.BottomHeight_cm(indices_Intervals)/100;
    sigma_z_bottom = WeightBSNE.BottomHeightErr_cm(indices_Intervals)/100;
    z_trapheight = WeightBSNE.HeightBSNE_cm(indices_Intervals)/100;
    
    %put all height info together into structured array, add to FluxBSNE
    z = struct('bottom', z_bottom, 'sigma_bottom', sigma_z_bottom, 'trapheight', z_trapheight, 'units',{'m'});
    FluxBSNE(i).z = z;
    
    %Compute fluxes, qz (g/m^2/s), and associated errors, sigma_qz (g/m^2/s) for each BSNE
    qz = zeros(N_ProfileBSNEs,1); %initialize list of qz's for profile
    sigma_qz = zeros(N_ProfileBSNEs,1); %initialize list of qz uncertainties for profile
    for j = 1:N_ProfileBSNEs %go through each BSNE in profile
        %compute mass uncertainty
        sigma_MassBSNE_g = ...
            sqrt((WeightBSNE.Weight_g(indices_Intervals(j))*sigma_FluxSampling).^2+...
            (WeightBSNE.WeightErr_g(indices_Intervals(j))).^2);
        
        %perform flux computation
        [qz_sample, sigma_qz_sample] = qzCalc(WeightBSNE.Weight_g(indices_Intervals(j)),...
            WeightBSNE.HeightBSNE_cm(indices_Intervals(j)),...
            WeightBSNE.WidthBSNE_cm(indices_Intervals(j)),...
            DurationBSNEs(j),...
            sigma_MassBSNE_g,...
            sigma_HeightWidthBSNE_cm,...
            sigma_HeightWidthBSNE_cm,...
            sigma_DurationBSNE);
        
        %add values to lists
        qz(j) = qz_sample;
        sigma_qz(j) = sigma_qz_sample;
    end
    
    %put all flux info together into structured array, add to FluxBSNE
    qz = struct('qz', qz, 'sigma', sigma_qz,'units',{'g/m^2/s'});
    FluxBSNE(i).qz = qz;
    
    %put all surface and profile grain size values together into structured array
    %d_Surface = struct('d_10',
    d = struct('d_10',d_10_profile,'d_50',d_50_profile,'d_90',d_90_profile,'IsCombined',d_IsCombined_profile,'StartTime',d_StartTime_profile','EndTime',d_EndTime_profile','units',{'mm'});
    FluxBSNE(i).d = d;
end

%% Iteratively compute BSNE heights, flux profiles, and uncertainties for each interval
for i = 1:N_Intervals
    %get flux and height profiles for fitting    
    qz_profile = FluxBSNE(i).qz.qz;
    z_bottom_profile = FluxBSNE(i).z.bottom; %bottom heights of BSNEs
    z_trapheight_profile = FluxBSNE(i).z.trapheight; %trap heights of BSNEs
    sigma_qz_profile = FluxBSNE(i).qz.sigma;
    sigma_z_profile = FluxBSNE(i).z.sigma_bottom;

    %start with guess of BSNE midpoint heights as arithmetic mean of traps
    z_profile = z_bottom_profile+z_trapheight_profile/2;
    
    %height-integrated flux from exponential fit
    [q0,zq,sigma_q0,sigma_zq,qz_pred,sigma_qz_pred,sigma_logqz_pred] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile);
    
    %now that we have zq, redo calculation of BSNE heights
    z_profile_old = z_profile; %document previous z-profile to see difference
    z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
    z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
    
    %iterate until the z_profile is minutely small
    while(z_profile_difference>1e-8)
        [q0,zq,sigma_q0,sigma_zq,qz_pred,sigma_qz_pred,sigma_logqz_pred] = qz_profilefit(qz_profile, z_profile, sigma_qz_profile, sigma_z_profile); %height-integrated flux from exponential fit
        z_profile_old = z_profile; %document previous z-profile to see difference
        z_profile = z_profile_calc(z_bottom_profile,z_trapheight_profile,zq); %calculate new BSNE midpoint heights based on zq
        z_profile_difference = mean(abs((z_profile-z_profile_old)./z_profile)); %get mean relative difference between profile heights
    end
   
    %% add to structured array
    FluxBSNE(i).qz.q0 = q0;
    FluxBSNE(i).qz.qz_pred = qz_pred;
    FluxBSNE(i).qz.sigma_q0 = sigma_q0;
    FluxBSNE(i).qz.sigma_qz_pred = sigma_qz_pred;
    FluxBSNE(i).qz.sigma_logqz_pred = sigma_logqz_pred;
    FluxBSNE(i).z.zq = zq;
    FluxBSNE(i).z.sigma_zq = sigma_zq;
    FluxBSNE(i).z.z = z_profile;
    FluxBSNE(i).z.sigma_z = sigma_z_profile;
    
    %% calculate Q and sigma_Q, add to structured array
    Q = q0*zq; %get total flux [g/m/s]
    sigma_Q = sqrt((sigma_q0*zq)^2+(sigma_zq*q0)^2); %estimate uncertainty in total flux
    FluxBSNE(i).Q = struct('Q', Q, 'sigma_Q', sigma_Q,'units',{'g/m/s'});
    
%     %% optional plot (comment out to hide)
% 
%     %get predicted values
%     qz_pred_minus = qz_pred - sigma_qz_pred;
%     qz_pred_plus = qz_pred + sigma_qz_pred;
%     [z_sort, ind_sort] = sort(z_profile);
%     qz_pred_sort = qz_pred(ind_sort);
%     qz_pred_minus_sort = qz_pred_minus(ind_sort);
%     qz_pred_plus_sort = qz_pred_plus(ind_sort);
% 
%     %make plot
%     figure(1); clf;
%     errorbar(z_profile,qz_profile,sigma_qz_profile,'bx'); hold on;
%     plot(z_sort,qz_pred_sort,'k');
%     plot(z_sort,qz_pred_minus_sort,'r--',z_sort,qz_pred_plus_sort,'r--');
%     set(gca,'yscale','log');
%     xlabel('z (m)');
%     ylabel('q (g/m^2/s)');
%     set(gca,'FontSize',16);
%     pause;
end

%reshape FluxBSNE so it is a column vector like all other structured arrays
FluxBSNE = FluxBSNE';