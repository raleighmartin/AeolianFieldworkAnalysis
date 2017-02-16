%% SCRIPT TO:
% 1. import wind data
% 2. find and remove error points
% 3. reorient wind vector along streamwise coordinate
% 4. compute values (mean wind, mean temperature, Reynolds stress, shear velocity, stability parameter)

%% PARAMETERS
%specify physical parameters
g = 9.8; %gravity acceleration, m/s^2
rho_a = 1.22; %air density, kg/m^3
zU = 0.64; %anemometer height, m
kappa = 0.4; %von Karman parameter

%specify filename for import
filename = 'Oceano_20150515_1610_1640_Wind_S1_z640mm.txt';

%specify year, month, and date
Year = 2015;
Month = 5;
Date = 15;

%specify maximum # of std devs in total wind away from mean wind for error detection
u_sigma_max = 5;

%% SCRIPT
%import raw wind data from file
[Hour,Minute,Second,u_raw,v_raw,w_raw,T_raw,diag] = import_wind(filename);

%generate times in datetime format
t_raw = datetime(Year,Month,Date,Hour,Minute,Second);

%determine which points have errors (and which ones don't)
[ind_wind_error, ind_wind_noerror] = find_ind_wind_error(u_raw, v_raw, w_raw, diag, u_sigma_max);

%keep only points with no errors
u = u_raw(ind_wind_noerror); %streamwise wind
v = v_raw(ind_wind_noerror); %lateral wind
w = w_raw(ind_wind_noerror); %vertical wind
t = t_raw(ind_wind_noerror); %times
T = T_raw(ind_wind_noerror); %temperatures

%reorient wind
[u_rot, v_rot, w_rot] = reorient_anemometers_vanboxel2004(u, v, w); %rotate instrument

%make computations
ubar = mean(u_rot) %mean streamwise wind
wbar = mean(w_rot); %mean vertical wind (should be ~0 after anemometer reorientation)
Tbar = mean(T); %mean temperature
uw = (u_rot-ubar).*(w_rot-wbar); %u'w' product
tauRe_kernal = mean(uw); %mean of u'w' product
if tauRe_kernal<=0
    tauRe = -rho_a*tauRe_kernal %Reynolds shear stress, Pa
    ustRe = sqrt(-tauRe_kernal) %shear velocity from Reynolds shear stress, m/s
else %if kernal is positive, then shear stress and velocity are undefined
    ustRe = NaN
    tauRe = NaN
end
Tw = mean((T-Tbar).*(w_rot-wbar)); %T'w' product
zL = (-(g./(Tbar+273.15)).*Tw)./(ustRe.^3./(kappa*zU)) %determine stability parameter, z/L