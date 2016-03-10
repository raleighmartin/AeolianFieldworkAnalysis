%% NONLINEAR FLUX LAW FITTING
% Q_fit = C u*^n (tau - tauit)
% n = u* scaling exponent (m/s)
% tauit = impact threshold stress (Pa)
% C = scaling coefficient
% tau = shear stress (Pa)
% Q = saltation flux (g/m/s)
% Chi2_min = min difference between Q and Q_fit

%% initialize
function [n, n_range_confidence, tauit, tauit_range_confidence, C, C_range_confidence, Chi2_min, Chi2_contributions, Q_pred] = ...
    NonlinearFluxFit(ust, sigma_ust, tau, sigma_tau, Q, sigma_Q)

%% Information for fitting
%ust exponent
n_fit_min = -4;
n_fit_max = 2;
n_fit_N = 601;
n_fit = linspace(n_fit_min,n_fit_max,n_fit_N);
%threshold stress
tauit_fit_min = 0;
tauit_fit_max = 0.2;
tauit_fit_N = 201;
tauit_fit = linspace(tauit_fit_min,tauit_fit_max,tauit_fit_N);
%scaling value
C_fit_min = 0;
C_fit_max = 300;
C_fit_N = 301;
C_fit = linspace(C_fit_min,C_fit_max,C_fit_N);

%% Basic info
rho_a = 1.18; %air density (kg/m^3)

%% Information for fitting
delta_Chi2_confidence = 2; %difference in Chi2 for 95% confidence interval

%% initialize Chi2 array
Chi2_matrix = zeros(n_fit_N, tauit_fit_N, C_fit_N);

%get Chi2 for all parameter combinations
for j=1:n_fit_N %cycling over n_fit values
    for k=1:tauit_fit_N %cycling over tauit_fit values
        for l=1:C_fit_N %cycling over C_fit values

            %predict value of Q for parameter combination
            Q_pred = C_fit(l)*ust.^n_fit(j).*(tau-tauit_fit(k)); 

            %compute Q uncertainty due to tau
            sigma_Q_tau = sigma_ust.*...
                abs(C_fit(l)*(rho_a*(n_fit(j)+2)*ust.^(n_fit(j)+1)-...
                n_fit(j)*ust.^(n_fit(j)-1)*tauit_fit(k))); %uncertainty in Q due to tau

            %compute total Q uncertainty
            sigma_Q_total = sqrt(sigma_Q.^2+sigma_Q_tau.^2); %compute total uncertainty in Q as combination of Q and tau uncertainty

            %compute chi-square for this set of parameters
            Chi2 = sum(((Q-Q_pred)./sigma_Q_total).^2); %compute Chi2 (Bevington and Robinson, Eq. 8.4)
            Chi2_matrix(j,k,l) = Chi2; %add this chi-square to array
        end 
    end
end

%determine indices that minimize chi_square
[Chi2_min,k_min]=min(min(min(Chi2_matrix)));
[~,j_min]=min(min(Chi2_matrix(:,:,k_min)));
[~,i_min]=min(Chi2_matrix(:,j_min,k_min));

if (Chi2_matrix(i_min,j_min,k_min)-min(min(min(Chi2_matrix)))~=0)
    error('something wrong')
end

%determine values corresponding to mimima for Chi2
n = n_fit(i_min);
tauit = tauit_fit(j_min);
C = C_fit(k_min);

%determine contributions of each point to Chi2 value
Q_pred = C*ust.^n.*(tau-tauit); %predict value of Q for parameter combination

%determine maximum Chi2 for confidence range calc
Chi2_confidence = Chi2_min + delta_Chi2_confidence;

%determine range of n
Chi2_n = zeros(n_fit_N,1);
for j = 1:n_fit_N
    Chi2_n(j) = min(min(Chi2_matrix(j,:,:)));
end
n_range_confidence = [min(n_fit(Chi2_n<=Chi2_confidence)), max(n_fit(Chi2_n<=Chi2_confidence))];

%determine range of tauit
Chi2_tauit = zeros(tauit_fit_N,1);
for k = 1:tauit_fit_N
    Chi2_tauit(k) = min(min(Chi2_matrix(:,k,:)));
end
tauit_range_confidence = [min(tauit_fit(Chi2_tauit<=Chi2_confidence)), max(tauit_fit(Chi2_tauit<=Chi2_confidence))];

%determine range of C
Chi2_C = zeros(C_fit_N,1);
for l = 1:C_fit_N
    Chi2_C(l) = min(min(Chi2_matrix(:,:,l)));
end
C_range_confidence = [min(C_fit(Chi2_C<=Chi2_confidence)), max(C_fit(Chi2_C<=Chi2_confidence))];

%% compute Chi2 contributions based on sigma_Q_total
sigma_Q_tau =  sigma_ust.*abs(C*(rho_a*(n+2)*ust.^(n+1)-n*ust.^(n-1)*tauit)); %get uncertainty in Q due to tau
sigma_Q_total = sqrt(sigma_Q.^2+sigma_Q_tau.^2); %compute total uncertainty in Q as combination of Q and tau uncertainty
Chi2_contributions = ((Q-Q_pred)./sigma_Q_total).^2; %compute Chi2 contributions