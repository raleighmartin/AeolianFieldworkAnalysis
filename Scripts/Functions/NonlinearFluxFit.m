%% INPUTS
% tau = shear stress (Pa)
% Q = saltation flux (g/m/s)
% sigma_Q = uncertainty in saltation flux (g/m/s)

%% OUTPUTS
% Q = C u*^n (tau - tauit)
% n = u* scaling exponent (m/s)
% tauit = impact threshold stress (Pa)
% C = scaling coefficient

function [n, tauit, C] = NonlinearFluxFit(tau, Q, sigma_tau, sigma_Q, N_fit)

%BASIC PARAMETERS
if nargin == 4
    N_fit = 50;
end
nu = length(Q)-3; %nu, #degrees of freedom
rho_a = 1.18; %air density (kg/m^3)
ust = sqrt(tau/rho_a); %u* (m/s)
sigma_ust = sigma_tau./(2*rho_a*ust); %uncertainty in u* (m/s)

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
C_fit_max = 500;
C_fit = linspace(C_fit_min,C_fit_max,N_fit);

%COMPUTE CHI2
%initialize Chi2 array
Chi2_fit = zeros(length(n_fit),length(tauit_fit),length(C_fit));

%get Chi2 for all parameter combinations
for i=1:length(n_fit) %cycling over n_fit values
    for j=1:length(tauit_fit) %cycling over tauit_fit values
        for k=1:length(C_fit) %cycling over C_fit values
            
            %compute Q for parameter combination
            Q_fit = C_fit(k)*ust.^(n_fit(i)).*(tau-tauit_fit(j));
            
            %compute total uncertainty in Q for parameter combination
            sigma_Q_ust = n_fit(i)*(sigma_ust./ust).*Q_fit;
            sigma_Q_tau = C_fit(k)*ust.^(n_fit(i)).*sigma_tau;
            sigma_Q_fit = sqrt(sigma_Q.^2+sigma_Q_ust.^2+sigma_Q_tau.^2);
            
            %compute chi-square from differences
            Chi2_fit(i,j,k) = Chi2Calculation(Q, sigma_Q_fit, Q_fit, nu);
        end 
    end
end

%determine indices that minimize chi_square
[~,k_min]=min(min(min(Chi2_fit)));
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