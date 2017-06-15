%% function to take info about BSNE trap weights ('WeightBSNE') and grain sizes ('GrainSize_BSNE')
% and aggregate this into mass and grain-size profiles
% Dependencies: IntervalsStartsWithin, qzCalc, BSNE_profilefit_exponential,
% qz_profilefit_exponential.m

function FluxBSNE = ProcessBSNEs(WeightBSNE,GrainSize_BSNE,zq_estimated_Site_m)

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
 
    %get BSNE lateral position
    y = WeightBSNE.Spanwise_m(indices_Intervals);
    FluxBSNE(i).y = y;
    
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
    d = struct('d_10',d_10_profile,'d_50',d_50_profile,'d_90',d_90_profile,'IsCombined',d_IsCombined_profile,'StartTime',d_StartTime_profile','EndTime',d_EndTime_profile','units',{'mm'});
    FluxBSNE(i).d = d;
end

%% Iteratively compute BSNE heights, flux profiles, and uncertainties for each interval
for i = 1:N_Intervals
    %get site name
    Site = FluxBSNE(i).Site;
    
    if strcmp(Site,'Oceano')
        %get indices of profile values on each side of tower
        ind_1 = find(FluxBSNE(i).y>=0); %indices of BSNEs to right of Wenglors (looking from upwind perspective)
        ind_2 = find(FluxBSNE(i).y<0); %indices of BSNEs to left of Wenglors (looking from upwind perspective)
        
        %separate flux and height profiles for fitting    
        qz_profile_1 = FluxBSNE(i).qz.qz(ind_1);
        qz_profile_2 = FluxBSNE(i).qz.qz(ind_2);
        z_bottom_profile_1 = FluxBSNE(i).z.bottom(ind_1); %bottom heights of BSNEs
        z_bottom_profile_2 = FluxBSNE(i).z.bottom(ind_2); %bottom heights of BSNEs
        z_trapheight_profile_1 = FluxBSNE(i).z.trapheight(ind_1); %trap heights of BSNEs
        z_trapheight_profile_2 = FluxBSNE(i).z.trapheight(ind_2); %trap heights of BSNEs
        sigma_qz_profile_1 = FluxBSNE(i).qz.sigma(ind_1);
        sigma_qz_profile_2 = FluxBSNE(i).qz.sigma(ind_2);
        sigma_z_profile_1 = FluxBSNE(i).z.sigma_bottom(ind_1);
        sigma_z_profile_2 = FluxBSNE(i).z.sigma_bottom(ind_2);

        %fit profile 1
        [z_profile_1,q0_1,zq_1,Q_1,sigma_q0_1,sigma_zq_1,sigma_Q_1,qz_pred_1,sigma_qz_pred_1,sigma_logqz_pred_1,sigma2_q0zq_1,...
            z_profile_geomean_1,q0_geomean_1,zq_geomean_1,Q_geomean_1] = ...
            BSNE_profilefit_exponential(qz_profile_1, z_bottom_profile_1, z_trapheight_profile_1, sigma_qz_profile_1, sigma_z_profile_1,zq_estimated_Site_m);

        %fit profile 2
        [z_profile_2,q0_2,zq_2,Q_2,sigma_q0_2,sigma_zq_2,sigma_Q_2,qz_pred_2,sigma_qz_pred_2,sigma_logqz_pred_2,sigma2_q0zq_2,...
            z_profile_geomean_2,q0_geomean_2,zq_geomean_2,Q_geomean_2] = ...
            BSNE_profilefit_exponential(qz_profile_2, z_bottom_profile_2, z_trapheight_profile_2, sigma_qz_profile_2, sigma_z_profile_2,zq_estimated_Site_m);
        
        %compute averaged values and their uncertainties
        if (~isnan(q0_1))&&(~isnan(q0_2))
            %q0 = mean([q0_1 q0_2]);
            %sigma_q0 = sqrt(sigma_q0_1.^2+sigma_q0_2.^2);
            q0 = sqrt(q0_1*q0_2); %geometric mean
            sigma_q0 = (1/2)*sqrt((q0_2/q0_1)*sigma_q0_1.^2+(q0_1/q0_2)*sigma_q0_2.^2);
            zq = mean([zq_1 zq_2]);
            sigma_zq = (1/2)*sqrt(sigma_zq_1.^2+sigma_zq_2.^2);
            %Q = mean([Q_1 Q_2]);
            %sigma_Q = sqrt(sigma_Q_1.^2+sigma_Q_2.^2);
            Q = sqrt(Q_1*Q_2); %geometric mean
            sigma_Q = (1/2)*sqrt((Q_2/Q_1)*sigma_Q_1.^2+(Q_1/Q_2)*sigma_Q_2.^2);
            %sigma2_q0zq = sigma2_q0zq_1 + sigma2_q0zq_2;
            sigma2_q0zq = sqrt(sigma2_q0zq_1*sigma2_q0zq_2);
        elseif ~isnan(q0_1)
            q0 = q0_1;
            sigma_q0 = sigma_q0_1;
            zq = zq_1;
            sigma_zq = sigma_zq_1;
            Q = Q_1;
            sigma_Q = sigma_Q_1;
            sigma2_q0zq = sigma2_q0zq_1;
        elseif ~isnan(q0_2)
            q0 = q0_2;
            sigma_q0 = sigma_q0_2;
            zq = zq_2;
            sigma_zq = sigma_zq_2;
            Q = Q_2;
            sigma_Q = sigma_Q_2;
            sigma2_q0zq = sigma2_q0zq_2;
        end
            
        %combine profile values based on original ordering
        N_profile = length(ind_1)+length(ind_2);
        qz_pred = zeros(N_profile,1);
        sigma_qz_pred = zeros(N_profile,1);
        sigma_logqz_pred = zeros(N_profile,1);
        z_profile = zeros(N_profile,1);
        sigma_z_profile = zeros(N_profile,1);
        qz_pred(ind_1) = qz_pred_1; qz_pred(ind_2) = qz_pred_2;
        sigma_qz_pred(ind_1) = sigma_qz_pred_1; sigma_qz_pred(ind_2) = sigma_qz_pred_2;
        sigma_logqz_pred(ind_1) = sigma_logqz_pred_1; sigma_logqz_pred(ind_2) = sigma_logqz_pred_2;
        z_profile(ind_1) = z_profile_1; z_profile(ind_2) = z_profile_2;
        sigma_z_profile(ind_1) = sigma_z_profile_1; sigma_z_profile(ind_2) = sigma_z_profile_2;
        
        %add combined values to structured array
        FluxBSNE(i).qz.q0 = q0;
        FluxBSNE(i).qz.sigma_q0 = sigma_q0;
        FluxBSNE(i).z.zq = zq;
        FluxBSNE(i).z.sigma_zq = sigma_zq;
        FluxBSNE(i).Q = struct('Q', Q, 'sigma_Q', sigma_Q, 'sigma2_q0zq', sigma2_q0zq, 'units',{'g/m/s'});
        FluxBSNE(i).qz.qz_pred = qz_pred;
        FluxBSNE(i).qz.sigma_qz_pred = sigma_qz_pred;
        FluxBSNE(i).qz.sigma_logqz_pred = sigma_logqz_pred;
        FluxBSNE(i).z.z = z_profile;
        FluxBSNE(i).z.sigma_z = sigma_z_profile;
        
        %add subprofile values to structured array
        FluxBSNE(i).qz.qz_1 = qz_profile_1;
        FluxBSNE(i).qz.qz_2 = qz_profile_2;
        FluxBSNE(i).qz.sigma_1 = sigma_qz_profile_1;
        FluxBSNE(i).qz.sigma_2 = sigma_qz_profile_2;
        FluxBSNE(i).qz.q0_1 = q0_1;
        FluxBSNE(i).qz.q0_2 = q0_2;
        FluxBSNE(i).qz.sigma_q0_1 = sigma_q0_1;
        FluxBSNE(i).qz.sigma_q0_2 = sigma_q0_2;
        FluxBSNE(i).z.zq_1 = zq_1;
        FluxBSNE(i).z.zq_2 = zq_2;
        FluxBSNE(i).z.sigma_zq_1 = sigma_zq_1;
        FluxBSNE(i).z.sigma_zq_2 = sigma_zq_2;
        FluxBSNE(i).Q.Q_1 = Q_1;
        FluxBSNE(i).Q.Q_2 = Q_2;
        FluxBSNE(i).Q.sigma_Q_1 = sigma_Q_1;
        FluxBSNE(i).Q.sigma_Q_2 = sigma_Q_2;
        FluxBSNE(i).Q.sigma2_q0zq_1 = sigma2_q0zq_1;
        FluxBSNE(i).Q.sigma2_q0zq_2 = sigma2_q0zq_2;
        FluxBSNE(i).qz.qz_pred_1 = qz_pred_1;
        FluxBSNE(i).qz.qz_pred_2 = qz_pred_2;
        FluxBSNE(i).qz.sigma_qz_pred_1 = sigma_qz_pred_1;
        FluxBSNE(i).qz.sigma_qz_pred_2 = sigma_qz_pred_2;
        FluxBSNE(i).qz.sigma_logqz_pred_1 = sigma_logqz_pred_1;
        FluxBSNE(i).qz.sigma_logqz_pred_2 = sigma_logqz_pred_2;
        FluxBSNE(i).z.z_1 = z_profile_1;
        FluxBSNE(i).z.z_2 = z_profile_2;
        FluxBSNE(i).z.sigma_z_1 = sigma_z_profile_1;
        FluxBSNE(i).z.sigma_z_2 = sigma_z_profile_2;
        
        %values for geometric mean computation
        FluxBSNE(i).qz.q0_geomean_1 = q0_geomean_1;
        FluxBSNE(i).qz.q0_geomean_2 = q0_geomean_2;
        FluxBSNE(i).z.zq_geomean_1 = zq_geomean_1;
        FluxBSNE(i).z.zq_geomean_2 = zq_geomean_2;
        FluxBSNE(i).Q.Q_geomean_1 = Q_geomean_1;
        FluxBSNE(i).Q.Q_geomean_2 = Q_geomean_2;
        FluxBSNE(i).z.z_geomean_1 = z_profile_geomean_1;
        FluxBSNE(i).z.z_geomean_2 = z_profile_geomean_2;
    else
        %get flux and height profiles for fitting    
        qz_profile = FluxBSNE(i).qz.qz;
        z_bottom_profile = FluxBSNE(i).z.bottom; %bottom heights of BSNEs
        z_trapheight_profile = FluxBSNE(i).z.trapheight; %trap heights of BSNEs
        sigma_qz_profile = FluxBSNE(i).qz.sigma;
        sigma_z_profile = FluxBSNE(i).z.sigma_bottom;

        %fit profile
        [z_profile,q0,zq,Q,sigma_q0,sigma_zq,sigma_Q,qz_pred,sigma_qz_pred,sigma_logqz_pred,sigma2_q0zq,...
            z_profile_geomean,q0_geomean,zq_geomean,Q_geomean] = ...
            BSNE_profilefit_exponential(qz_profile, z_bottom_profile, z_trapheight_profile, sigma_qz_profile, sigma_z_profile, zq_estimated_Site_m);

        %add to structured array
        FluxBSNE(i).qz.q0 = q0;
        FluxBSNE(i).qz.sigma_q0 = sigma_q0;
        FluxBSNE(i).z.zq = zq;
        FluxBSNE(i).z.sigma_zq = sigma_zq;
        FluxBSNE(i).Q = struct('Q', Q, 'sigma_Q', sigma_Q, 'sigma2_q0zq', sigma2_q0zq, 'units',{'g/m/s'});
        FluxBSNE(i).qz.qz_pred = qz_pred;
        FluxBSNE(i).qz.sigma_qz_pred = sigma_qz_pred;
        FluxBSNE(i).qz.sigma_logqz_pred = sigma_logqz_pred;
        FluxBSNE(i).z.z = z_profile;
        FluxBSNE(i).z.sigma_z = sigma_z_profile;
        
        %values for geometric mean computation
        FluxBSNE(i).qz.q0_geomean = q0_geomean;
        FluxBSNE(i).z.zq_geomean = zq_geomean;
        FluxBSNE(i).Q.Q_geomean = Q_geomean;
        FluxBSNE(i).z.z_geomean = z_profile_geomean;
    end
end

%reshape FluxBSNE so it is a column vector like all other structured arrays
FluxBSNE = FluxBSNE';