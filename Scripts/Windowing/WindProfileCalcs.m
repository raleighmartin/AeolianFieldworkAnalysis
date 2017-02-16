%% SCRIPT TO GENERATE WINDOWS OF FLUX ACTIVITY, SALTATION FLUX, AND STRESS VALUES FOR ANALYSIS

%% initialize
clearvars;
close all;

%% parameter values
u_sigma_max = 5; %maximum standard deviation in total wind for error detection
rho_a = [1.16, 1.22, 1.22]; %air density kg/m^3 (assumes T~30 C at Jeri and ~15 C at Rancho and Oceano)
kappa = 0.4; %von Karman parameter
g = 9.8; %gravity m/s^2
zU_max = 2.5; %maximum anemometer height for profile fit (m)

%% information about where to load/save data, plots, and functions
folder_TimeData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_ProcessedData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for windowing data - restricted windows
LoadTimeData_Path = strcat(folder_TimeData,'TimeWindows_30min_Restricted'); %path for loading time windows
SaveDataWindows_Path = strcat(folder_ProcessedData,'WindProfiles_30min_Restricted'); %path for saving output data
min_datetime = []; %no restrictions on dates

% %% Specific information for windowing data - unrestricted windows
% LoadTimeData_Path = strcat(folder_TimeData,'TimeWindows_30min_Unrestricted'); %path for loading time windows
% SaveDataWindows_Path = strcat(folder_ProcessedData,'WindProfiles_30min_Unrestricted'); %path for saving output data
% min_datetime = []; %no restrictions on dates
  
%% load time windows
load(LoadTimeData_Path);

%% load processed data for each site, add to cell arrays of all sites 
WindData_profile = cell(N_Sites,1);
N_Anemometers = zeros(N_Sites,1); %number of anemometers

for i = 1:N_Sites
    ProcessedData_Path = strcat(folder_ProcessedData,'ProcessedData_',Sites{i});
    load(ProcessedData_Path); %load processed data
    
    %get wind data
    WindData = ProcessedData.(AnemometerType{i}); 
    clear ProcessedData; %remove 'ProcessedData' to clear up memory
    
    %get wind profile data
    N_Anemometers(i) = length(AnemometerName_profile{i}); %number of anemometers
    WindData_profile{i} = cell(N_Anemometers(i),1);
    for k = 1:length(AnemometerName_profile{i})
        WindData_profile{i}{k} = WindData.(AnemometerName_profile{i}{k}); %data for anemometer profile
    end
    clear WindData; %remove 'WindData' to clear up memory
end

%% initialize variable lists
theta_profile_window = cell(N_Sites,1); %wind angle for time window
u_profile_window = cell(N_Sites,1); %rotated streamwise wind for time window
v_profile_window = cell(N_Sites,1); %rotated transverse wind for time window
w_profile_window = cell(N_Sites,1); %rotated vertical wind for time window
Temp_profile_window = cell(N_Sites,1); %temperature for time window
ustRe_profile_window = cell(N_Sites,1); %Reynolds shear velocity at all heights of wind profile
sigma_ustRe_profile_window = cell(N_Sites,1); %Reynolds shear velocity at all heights of wind profile - uncertainty
tauRe_profile_window = cell(N_Sites,1); %Reynolds shear stress at all heights of wind profile
sigma_tauRe_profile_window = cell(N_Sites,1); %Reynolds shear stress at all heights of wind profile - uncertainty
zL_profile_window = cell(N_Sites,1); %stability parameter at all heights of wind profile

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
        
    %% get start times and end times
    WindowStartTimes = StartTime_window{i};
    WindowEndTimes = EndTime_window{i};
    N_Windows = length(WindowStartTimes);

    %% initialize profile values
    theta_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %wind angles for anemometer
    u_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %rotated streamwise wind for time window
    v_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %rotated transverse wind for time window
    w_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %rotated vertical wind for time window
    Temp_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %temperature for time window
    ustRe_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %Reynolds shear velocity at all heights of wind profile
    sigma_ustRe_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %Reynolds shear velocity at all heights of wind profile - uncertainty
    tauRe_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %Reynolds shear stress at all heights of wind profile
    sigma_tauRe_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %Reynolds shear stress at all heights of wind profile - uncertainty
    zL_profile_window{i} = zeros(N_Windows,N_Anemometers(i))*NaN; %stability parameter at all heights of wind profile    
        
    %% go through time windows
    for j = 1:N_Windows

        %display processing status
        processing_status = [Sites{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]

        %get specific start and end time
        StartTime = WindowStartTimes(j);
        EndTime = WindowEndTimes(j);

        %get duration of interval in seconds
        T_interval = seconds(EndTime-StartTime);
                               
        %% got through anemometers in profile
        for k = 1:N_Anemometers(i)

            %only perform extraction if there are data in window
            if ~isnan(zU_profile_window{i}(j,k))

                %extract time interval
                [~, ~, IntervalN, IntervalInd] = ExtractVariableTimeInterval(WindData_profile{i}{k},StartTime,EndTime,'u','int','int');

                %use only longest interval
                ind_longest = find(cellfun(@length,IntervalInd)==max(cellfun(@length,IntervalInd)));
                IntervalN = IntervalN(ind_longest); %number of interval
                IntervalInd = IntervalInd{ind_longest}; %indices for interval

                %further reduce IntervalInd based on eliminating error times
                [~, ErrInd, ~] = intersect(WindData_profile{i}{k}(IntervalN).t.int,WindData_profile{i}{k}(IntervalN).t.err);
                if strcmp(AnemometerType{i},'Sonic') %add in additional points from diagnostic flag if instrument is sonic
                    diag_ind = find(WindData_profile{i}{k}(IntervalN).diag.raw~=0); %get diagnostic points for interval
                    ErrInd = union(ErrInd,diag_ind); %add these to error points
                end
                IntervalInd_noerr = setdiff(IntervalInd,ErrInd); %indices for interval excluding error times

                %get velocity values using no error points
                u = WindData_profile{i}{k}(IntervalN).u.int(IntervalInd_noerr);
                v = WindData_profile{i}{k}(IntervalN).v.int(IntervalInd_noerr);
                w = WindData_profile{i}{k}(IntervalN).w.int(IntervalInd_noerr);
                Temp = WindData_profile{i}{k}(IntervalN).T.int(IntervalInd_noerr);

                %remove additional error points based on large deviations in wind velocity
                u_total = sqrt(u.^2+v.^2+w.^2); %total wind
                u_total_max = mean(u_total)+u_sigma_max*std(u_total); %maximum total wind based on multiple of std dev
                ind_good_wind = find(u_total<=u_total_max); %get indices of points with total wind below upper limit
                u = u(ind_good_wind); %keep only good u
                v = v(ind_good_wind); %keep only good v
                w = w(ind_good_wind); %keep only good w
                Temp = Temp(ind_good_wind); %keep only good temperatures

                %compute wind angle
                theta = atan(mean(v./u))*180/pi; %wind angle

                %rotate instrument, call these 'rot' values
                [u_rot, v_rot, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

                %get mean values
                ubar = mean(u_rot);
                vbar = mean(v_rot);
                wbar = mean(w_rot);
                Tempbar = mean(Temp);

                %add to window lists - mean values
                theta_profile_window{i}(j,k) = theta; %wind angle
                u_profile_window{i}(j,k) = ubar; %rotated streamwise wind
                v_profile_window{i}(j,k) = vbar; %rotated transverse wind
                w_profile_window{i}(j,k) = wbar; %rotated vertical wind
                Temp_profile_window{i}(j,k) = Tempbar; %temperature

                %make shear stress computation computation
                uw = (u_rot-mean(u_rot)).*(w_rot-mean(w_rot)); %u'w' product
                tauRe_kernal = mean(uw);
                if tauRe_kernal<=0
                    ustRe = sqrt(-tauRe_kernal);
                    tauRe = -rho_a(i)*tauRe_kernal;
                else
                    tauRe = NaN;
                    ustRe = NaN;
                end

                %compute stress uncertainty using Marcelo's code - use only error-free data
                [sigma_uw_bar, ~] = random_error(uw, length(uw), 1/dt_wind_window(i), zU_profile_window{i}(j,k), mean(u_rot), 0, 0, 0, 1); %use script from Salesky et al (2012)

                %add shear velocity values to list
                ustRe_profile_window{i}(j,k) = ustRe; %u*
                sigma_ustRe_profile_window{i}(j,k) = sigma_uw_bar/(2*ustRe); %uncertainty in u*

                %add shear stress values to list
                tauRe_profile_window{i}(j,k) = tauRe; %Reynolds shear stress at all heights of wind profile
                sigma_tauRe_profile_window{i}(j,k) = rho_a(i)*sigma_uw_bar; %Reynolds shear stress at all heights of wind profile - uncertainty

                %make stability computation
                heatflux = mean((Temp-Tempbar).*(w-wbar));
                zL = (-(g./(Tempbar+273.15)).*heatflux)./(ustRe.^3./(kappa*zU_profile_window{i}(j,k)));

                %add stability parameter at all heights to list
                zL_profile_window{i}(j,k) = zL;
            end
        end
    end
end

% SAVE DATA
save(SaveDataWindows_Path,'Site*','N_Sites','AnemometerName_profile','*window'); %save primary data reduced file in ProcessedData folder