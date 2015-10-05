%information about beam and background flux
d_b = 0.6; %beam diameter, mm
f_s = 0.9; %fractional sensitivity
n_arrival = 1000; %total arrival rate for particle size, #/mm^2/s

%particle diameters (mm)
N_d = 50; %number of diameters to use
mu_lognormal = -1.25; %approximate lognormal mu from d50
sigma_lognormal = 0.35; %approximate lognormal sigma from difference in these diameters
d_p = logspace(log10(0.1),log10(1),N_d)'; %generate log-spaced list of diameters with +/-3 sigma
d_pdf = (1./(d_p*sigma_lognormal*sqrt(2*pi))).*exp(-(log(d_p)-mu_lognormal).^2/(2*sigma_lognormal^2)); %pdf values for diameters
d_f = d_pdf/sum(d_pdf); %fractional values for diameters
d10 = min(d_p(cumsum(d_f)>=0.1))
d50 = min(d_p(cumsum(d_f)>=0.5))
d90 = min(d_p(cumsum(d_f)>=0.9))

%initialize list of detection rates
n_detect = zeros(N_d,1)
for i = 1:N_d
    delta_test = 0:
    %delta_max = ??; %maximum distance of particle center from beam center for detection
    A_detect = (pi/4)*delta_max.^2; %detection area, mm^2
    n_detect = n_arrival*A_detect; %detection rate for particle size

% plot detection rate versus particle diameter
figure(1); clf;
plot(d_p,n_detect);

GOAL: number of particles detected
% n_b = ??
% f_p = 