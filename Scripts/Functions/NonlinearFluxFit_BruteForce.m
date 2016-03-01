%% INPUTS
% tau = shear stress (Pa)
% Q = saltation flux (g/m/s)
% sigma_Q = uncertainty in saltation flux (g/m/s)

%% OUTPUTS
% Q = C u*^n (tau - tauit)
% n = u* scaling exponent (m/s)
% tauit = impact threshold stress (Pa)
% C = scaling coefficient

function [Chi2, n, tauit, C, sigma_n, sigma_tauit, sigma_C] = NonlinearFluxFit_BruteForce(tau, Q, sigma_tau, sigma_Q, Site, N_fit)

%BASIC PARAMETERS
if nargin == 5
    N_fit = 400;
end
nu = length(Q)-3; %nu, #degrees of freedom
rho_a = 1.18; %air density (kg/m^3)
ust = sqrt(tau/rho_a); %u* (m/s)
delta_Chi2_uncertainty = 1; %difference in Chi2 for uncertainty estimation

%FITTING PARAMATERS
%ust exponent
n_fit_min = -2;
n_fit_max = 2;
n_fit = linspace(n_fit_min,n_fit_max,N_fit);
%threshold stress
tauit_fit_min = 0;
tauit_fit_max = 0.2;
tauit_fit = linspace(tauit_fit_min,tauit_fit_max,N_fit);
%scaling value
C_fit_min = 0;
C_fit_max = 300;
C_fit = linspace(C_fit_min,C_fit_max,N_fit);

%COMPUTE CHI2
%initialize Chi2 array
Chi2_fit = zeros(length(n_fit),length(tauit_fit),length(C_fit));

%get Chi2 for all parameter combinations
for i=1:length(n_fit) %cycling over n_fit values
    for j=1:length(tauit_fit) %cycling over tauit_fit values
        for k=1:length(C_fit) %cycling over C_fit values
                        
            %compute chi-square for this set of parameters
            [~, Chi2] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_fit(i), tauit_fit(j), C_fit(k));
            Chi2_fit(i,j,k) = Chi2;
        end 
    end
end

%determine indices that minimize chi_square
[Chi2_min,k_min]=min(min(min(Chi2_fit)));
[~,j_min]=min(min(Chi2_fit(:,:,k_min)));
[~,i_min]=min(Chi2_fit(:,j_min,k_min));

if (Chi2_fit(i_min,j_min,k_min)-min(min(min(Chi2_fit)))~=0)
    error('something wrong')
end

%COMPUTE FIT VALUES
%determine values corresponding to mimima for chi_square
n = n_fit(i_min);
tauit = tauit_fit(j_min);
C = C_fit(k_min);
Chi2 = Chi2_min;

%determine maximum Chi2 for uncertainty
Chi2_max = Chi2+delta_Chi2_uncertainty;

%determine range of n
Chi2_n = zeros(N_fit,1);
for i = 1:N_fit
    Chi2_n(i) = min(min(Chi2_fit(i,:,:)));
end
n_range = [min(n_fit(Chi2_n<=Chi2_max)), max(n_fit(Chi2_n<=Chi2_max))];
sigma_n = range(n_range)/2;

%determine range of tauit
Chi2_tauit = zeros(N_fit,1);
for j = 1:N_fit
    Chi2_tauit(j) = min(min(Chi2_fit(:,j,:)));
end
tauit_range = [min(tauit_fit(Chi2_tauit<=Chi2_max)), max(tauit_fit(Chi2_tauit<=Chi2_max))];
sigma_tauit = range(tauit_range)/2;

%determine range of C
Chi2_C = zeros(N_fit,1);
for k = 1:N_fit
    Chi2_C(k) = min(min(Chi2_fit(:,:,k)));
end
C_range = [min(C_fit(Chi2_C<=Chi2_max)), max(C_fit(Chi2_C<=Chi2_max))];
sigma_C = range(C_range)/2;

%save
savename = strcat('NonlinearFluxFit_',Site);
save(savename, 'C', 'C_fit', 'sigma_C', 'C_range',...
    'n', 'n_fit', 'sigma_n', 'n_range',...
    'tauit', 'tauit_fit', 'sigma_tauit', 'tauit_range',...
    'Chi2','Chi2_fit');