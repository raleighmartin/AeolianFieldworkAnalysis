%% SCRIPT TO ANALYZE THRESHOLD USING TFEM METHOD

%% initialize
clearvars;

%% parameter values
u_ft = [9, 7.9, 7.3]; %estimated fluid threshold wind for each site
u_it = [7.2, 6.3, 5.8]; %estimated impact threshold wind for each site
rho_a = 1.18; %air density kg/m^3
g = 9.8; %gravity m/s^2
kappa = 0.4; %von Karman parameter
z0 = 1e-4; %aerodynamic roughness length (m)
fD_min = 0.005; %minimum detection rate for zero flux

%% information about where to load/save data, plots, and functions
folder_LoadData = '../../../../Google Drive/Data/AeolianFieldwork/Processed/'; %folder for retrieving processed data
folder_SaveData = '../../AnalysisData/Windowing/'; %folder for outputs of this analysis
folder_Functions = '../Functions/'; %folder with functions
addpath(folder_Functions); %point MATLAB to location of functions

%% Specific information for window data
LoadData_Path = strcat(folder_LoadData,'WindowAverageWindows_30min'); %path for loading window data
SaveData_Path = strcat(folder_SaveData,'ActivityWindows_30min.mat'); %path for saving output data

%% load window data
load(LoadData_Path);

%% initialize variable list
fD_window = cell(N_Sites,1); %Wenglor detection frequency list for each window
fQ_window = cell(N_Sites,1); %Wenglor transport frequency matrix for each window
uth_fQ_window = cell(N_Sites,1); %threshold wind estimated from histogram of winds and fQ
ustth_fQ_window = cell(N_Sites,1); %threshold shear velocity estimated uth and law-of-the-wall 
tauth_fQ_window = cell(N_Sites,1); %threshold wind stress estimated uth and law-of-the-wall 
fplus_window = cell(N_Sites,1); %fraction of time with u_avg above fluid threshold
fminus_window = cell(N_Sites,1); %fraction of time with u_avg below impact threshold
fint_window = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone
fint_down_window = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone from above
fint_up_window = cell(N_Sites,1); %fraction of time with u_avg in intermediate zone from below

%% PERFORM ANALYSIS FOR EACH SITE
for i = 1:N_Sites
    
    N_Windows = length(n_avg_window{i}); %get number of windows
    
    %% initialize lists of values
    fD_window{i} = zeros(N_Windows,N_avg_window); %Wenglor detection frequency list for each window
    fQ_window{i} = zeros(N_Windows,N_avg_window); %Wenglor transport frequency matrix for each window
    uth_fQ_window{i} = zeros(N_Windows,N_avg_window); %threshold wind estimated from histogram of winds and fQ
    ustth_fQ_window{i} = zeros(N_Windows,N_avg_window); %threshold shear velocity estimated uth and law-of-the-wall 
    tauth_fQ_window{i} = zeros(N_Windows,N_avg_window); %threshold wind stress estimated uth and law-of-the-wall 
    fplus_window{i} = zeros(N_Windows,N_avg_window); %fraction of time with u_avg above fluid threshold
    fminus_window{i} = zeros(N_Windows,N_avg_window); %fraction of time with u_avg below impact threshold
    fint_window{i} = zeros(N_Windows,N_avg_window); %fraction of time with u_avg in intermediate zone
    fint_down_window{i} = zeros(N_Windows,N_avg_window); %fraction of time with u_avg in intermediate zone from above
    fint_up_window{i} = zeros(N_Windows,N_avg_window); %fraction of time with u_avg in intermediate zone from below
    
    %% go through time windows
    for j = 1:N_Windows
        
        %display processing status
        processing_status = [SiteNames{i},', ',int2str(j),' of ',int2str(N_Windows),', ',datestr(now)]
        
        %% go through averaging times
        for k = 1:N_avg_window
            
            %extract relevant info
            t = t_avg_window{i}{j,k}; %times for window average
            n = n_avg_window{i}{j,k}; %window-averaged count rates values
            u = u_avg_window{i}{j,k}; %window-averaged u values
            ind_flux_err = ind_flux_err_avg_window{i}{j,k}; %indices of error points in flux window-average
            ind_wind_err = ind_wind_err_avg_window{i}{j,k}; %indices of error points in wind window-average

            %determine which points to include in calculation
            ind_err = union(ind_flux_err,ind_wind_err); %indices of all error points
            ind_noerr = setdiff((1:length(t))',ind_err); %indices of nonerror points
            T_noerr = length(ind_noerr); %number of timesteps
            t_noerr = t(ind_noerr); %times for nonerror points
            n_noerr = n(ind_noerr,:); %counts for nonerror points
            u_noerr = u(ind_noerr); %winds for nonerror points

            %calculate flux detection frequency and mean particle counts for no error points
            D = sum(n_noerr')'>0; %particle detections (1 if detected, 0 if not)
            fD = sum(D)/T_noerr; %calculate detection rate (only non-error points)
            n_avg = mean(n_noerr); %mean particle counts rate profile

            %Determine flux frequency iteratively
            if fD<=fD_min %if no (or negligible) flux, set frequencies to zero
                fQ = 0;
            else %otherwise, do calculations
                %estimate fQ iteratively 
                fQ_prev = 0; %previous estimate of fQ (arbitrarily set to 0 for while loop functioning) 
                fQ = fD; %start by assuming fQ=fD
                while abs(fQ-fQ_prev)>0.001 %compare current and previous estimations of fQ for while loop
                    fQ_prev = fQ; %update fQ_prev with last fQ
                    lambda = sum(n_avg)/(fQ*seconds(delta_t_avg_window(k))); %estimated particle passage rate during transport
                    fQ = fD/(1-exp(-lambda)); %calculate new fQ
                end

                %remove points with fQ>1
                if fQ>1
                    fQ = NaN;
                end
            end

            %determine wind corresponding to fQ
            u_sort = sort(u_noerr); %sort u's
            ind_uth_fQ = round((1-fQ)*T_noerr); %get index in list of u's corresponding to threshold
            if ind_uth_fQ==0||isnan(ind_uth_fQ)
                uth_fQ = NaN;
                ustth_fQ = NaN;
                tauth_fQ = NaN;
            else
                uth_fQ = u_sort(ind_uth_fQ); %threshold wind speed
                ustth_fQ = (kappa*uth_fQ)/log(zU_window{i}(j)/z0); %threshold shear velocity
                tauth_fQ = rho_a*ustth_fQ^2; %threshold shear stress
            end

            %determine fractions of u in different regions
            ind_fplus = find(u_noerr>=u_ft(i)); %indices for u above u_ft
            ind_fminus = find(u_noerr<u_it(i)); %indices for u below u_it
            ind_fint = intersect(find(u_noerr>=u_it(i)),find(u_noerr<u_ft(i))); %indices for u in intermediate zone
            N_fplus = length(ind_fplus); %number of u above u_ft
            N_fminus = length(ind_fminus); %number of u below u_it
            N_fint = length(ind_fint); %number of u in intermediate zone
            fplus = N_fplus/T_noerr; %fraction of time with u above fluid threshold
            fminus = N_fminus/T_noerr; %fraction of time with u below impact threshold
            fint = N_fint/T_noerr; %fraction of time with u in intermediate zone

            %look at hysteresis for u in intermediate zone
            N_fint_up = 0; %initialize N_fint_up (number of upcrossings to hysteresis zone)
            N_fint_down = 0; %initialize N_fint_down (number of downcrossings to hysteresis zone)
            N_fint_unknown = 0; %initialize N_fint_isempty (number of events in hysteresis zone with no known history)
            for l = 1:N_fint;
                ind_last_fminus = find(ind_fminus<ind_fint(l), 1, 'last' );
                ind_last_fplus = find(ind_fplus<ind_fint(l), 1, 'last' );
                if(isempty(ind_last_fplus))
                    if(isempty(ind_last_fminus))
                        N_fint_unknown = N_fint_unknown+1;
                    else
                        N_fint_up = N_fint_up+1;
                    end
                elseif(isempty(ind_last_fminus))
                    if(isempty(ind_last_fplus))
                        N_fint_unknown = N_fint_unknown+1;
                    else
                        N_fint_down = N_fint_down+1;
                    end
                else
                    if(ind_last_fminus>ind_last_fplus)
                        N_fint_up = N_fint_up+1;
                    elseif(ind_last_fplus>ind_last_fminus)
                        N_fint_down = N_fint_down+1;
                    end
                end
            end
            N_fint_known = N_fint_up+N_fint_down;
            fint_up = (N_fint_up/N_fint_known)*fint;
            fint_down = (N_fint_down/N_fint_known)*fint;

            %add to lists of all values
            fD_window{i}(j,k) = fD; %detection frequency
            fQ_window{i}(j,k) = fQ; %flux frequency
            uth_fQ_window{i}(j,k) = uth_fQ; %threshold wind from flux frequency
            tauth_fQ_window{i}(j,k) = tauth_fQ; %threshold wind from flux frequency
            ustth_fQ_window{i}(j,k) = ustth_fQ; %threshold wind from flux frequency
            fplus_window{i}(j,k) = fplus; %fraction of time with u_avg above fluid threshold
            fminus_window{i}(j,k) = fminus; %fraction of time with u_avg below impact threshold
            fint_window{i}(j,k) = fint; %fraction of time with u_avg in intermediate zone
            fint_down_window{i}(j,k) = fint_down; %fraction of time with u_avg in intermediate zone with downward crossing
            fint_up_window{i}(j,k) = fint_up; %fraction of time with u_avg in intermediate zone with upward crossing
        end
    end
end

% SAVE DATA
save(SaveData_Path,'SiteNames','N_Sites','*window'); %save reduced file in GitHub folder