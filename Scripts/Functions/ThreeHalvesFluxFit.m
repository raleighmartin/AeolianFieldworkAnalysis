%% NONLINEAR FLUX LAW FITTING
% Q_fit = C u*^n (tau - tauit)
% n = u* scaling exponent (m/s)
% tauit = impact threshold stress (Pa)
% C = scaling coefficient
% tau = shear stress (Pa)
% Q = saltation flux (g/m/s)
% Chi2_min = min difference between Q and Q_fit

%% initialize
function [tauit, tauit_sigma, C, C_sigma, Chi2_min, Chi2_contributions, Q_pred] = ...
    ThreeHalvesFluxFit(ust, sigma_ust, tau, sigma_tau, Q, sigma_Q)

%% Information for fitting
%threshold stress
tauit_fit_min = 0;
tauit_fit_max = 0.2;
tauit_fit_N = 201;
tauit_fit = linspace(tauit_fit_min,tauit_fit_max,tauit_fit_N);
%scaling value
C_fit_min = 0;
C_fit_max = 600;
C_fit_N = 601;
C_fit = linspace(C_fit_min,C_fit_max,C_fit_N);

%% Basic info
rho_a = 1.18; %air density (kg/m^3)

%% Information for fitting
delta_Chi2_sigma = 1; %difference in Chi2 for sigma interval

%% initialize Chi2 array
Chi2_matrix = zeros(tauit_fit_N, C_fit_N);

%get Chi2 for all parameter combinations
for k=1:tauit_fit_N %cycling over tauit_fit values
    for l=1:C_fit_N %cycling over C_fit values

        %predict value of Q for parameter combination
        Q_pred = C_fit(l)*ust.*(tau-tauit_fit(k)); 

        %compute Q uncertainty due to tau
        sigma_Q_tau = sigma_ust.*C_fit(l).*abs(3*tau-tauit_fit(k)); %uncertainty in Q due to tau

        %compute total Q uncertainty
        sigma_Q_total = sqrt(sigma_Q.^2+sigma_Q_tau.^2); %compute total uncertainty in Q as combination of Q and tau uncertainty

        %compute chi-square for this set of parameters
        Chi2 = sum(((Q-Q_pred)./sigma_Q_total).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
        Chi2_matrix(k,l) = Chi2; %add this chi-square to array
    end 
end

%determine indices that minimize chi_square
[Chi2_min,l_min]=min(min(Chi2_matrix));
[~,k_min]=min(Chi2_matrix(:,l_min));

if (Chi2_matrix(k_min,l_min)-min(min(Chi2_matrix))~=0)
    error('something wrong')
end

%determine values corresponding to mimima for Chi2
tauit = tauit_fit(k_min);
C = C_fit(l_min);

%determine contributions of each point to Chi2 value
Q_pred = C*ust.*(tau-tauit); %predict value of Q for parameter combination

%determine maximum Chi2 for sigma range calc
Chi2_sigma = Chi2_min + delta_Chi2_sigma;

%determine range of tauit
Chi2_tauit = zeros(tauit_fit_N,1);
for k = 1:tauit_fit_N
    Chi2_tauit(k) = min(Chi2_matrix(k,:));
end
tauit_sigma = range([min(tauit_fit(Chi2_tauit<=Chi2_sigma)), max(tauit_fit(Chi2_tauit<=Chi2_sigma))])/2;

%determine range of C
Chi2_C = zeros(C_fit_N,1);
for l = 1:C_fit_N
    Chi2_C(l) = min(Chi2_matrix(:,l));
end
C_sigma = range([min(C_fit(Chi2_C<=Chi2_sigma)), max(C_fit(Chi2_C<=Chi2_sigma))])/2;

%% compute Chi2 contributions based on sigma_Q_total
sigma_Q_tau = sigma_ust.*C.*abs(3*tau-tauit); %uncertainty in Q due to tau
sigma_Q_total = sqrt(sigma_Q.^2+sigma_Q_tau.^2); %compute total uncertainty in Q as combination of Q and tau uncertainty
Chi2_contributions = ((Q-Q_pred)./sigma_Q_total).^2; %compute Chi2 contributions