%% INPUTS
% tau = shear stress (Pa)
% Q = saltation flux (g/m/s)
% sigma_Q = uncertainty in saltation flux (g/m/s)

%% OUTPUTS
% Q = C u*^n (tau - tauit)
% n = u* scaling exponent (m/s)
% tauit = impact threshold stress (Pa)
% C = scaling coefficient

function [n_fit, tauit_fit, C_fit, Chi2_fit] = NonlinearFluxFit_Iterative(tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess)

%PHYSICAL CALCULATIONS
rho_a = 1.18; %air density (kg/m^3)
ust = sqrt(tau/rho_a); %u* (m/s)

%FITTING PARAMETERS
if nargin == 4
    n_guess = 0; %initial n for fitting
    tauit_guess = 0.1; %initial tauit for fitting
    C_guess = 200; %initial C for fitting
end
delta_n = 0.0005; %n increment for fitting
delta_tauit = 0.0001; %tauit increment for fitting
delta_C = 0.05; %C increment for fitting
delta_Chi2_min = 0.001; %change in Chi2 for convergence
delta_Chi2_uncertainty = 1; %difference in Chi2 for uncertainty estimation
nu = length(Q)-3; %nu, #degrees of freedom

%COMPUTE INITIAL CHI2
[~, Chi2_guess] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess); %get value of chi2 for initial guess
delta_Chi2 = Inf; %initialize delta_Chi2, the change in Chi2 from previous iteration
delta2_Chi2 = Inf; %initialize delta2_Chi2, the change in Chi2 from two iterations ago
Chi2_previous = Inf; %initialize Chi2_previous, record of Chi2 in last iteration

%RECORD PARAMETER VALUES
n_all = n_guess;
tauit_all = tauit_guess;
C_all = C_guess;
Chi2_all = Chi2_guess;

%ITERATE FITTING FOR EACH PARAMETER UNTIL MINIMUM CHI2 IS REACHED
while (delta2_Chi2~=0)&&... %iteration will finish if stuck in loop (i.e., chi2 is same as two iterations ago)
    ((abs(delta_Chi2)>delta_Chi2_min)||... %or if chi2 reduction for iteration is sufficiently small
    (delta_Chi2>0));... %and this reduction is negative
    Chi2_previous2 = Chi2_previous; %record Chi2 of 2 iterations ago
    Chi2_previous = Chi2_guess; %record previous Chi2 for later comparison
    
    %FIND N THAT MINIMIZES CHI2
    %determine direction of n search to minimize chi2
    [~, Chi2_n_plus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess+delta_n, tauit_guess, C_guess); %guess value of Q for parameter combination with n+delta_n
    [~, Chi2_n_minus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess-delta_n, tauit_guess, C_guess); %guess value of Q for parameter combination with n-delta_n
    sign_delta_n = -sign(Chi2_n_plus-Chi2_n_minus); %direction of search for value of n
    
    %increment n until Chi2 reaches minimum
    Chi2_previous_n = Chi2_guess; %record previous Chi2 for later comparison
    Chi2_guess = min([Chi2_n_plus,Chi2_n_minus]); %get new Chi2 guess from min of adjacent Chi2
    delta_Chi2_n = Chi2_guess - Chi2_previous_n; %calculate Chi2 reduction
    while delta_Chi2_n < 0
        Chi2_previous_n = Chi2_guess; %record previous Chi2 for later comparison
        n_guess = n_guess+delta_n*sign_delta_n; %increment n_guess
        [~, Chi2_guess] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess); %determine Chi2 for new n_guess
        delta_Chi2_n = Chi2_guess - Chi2_previous_n; %calculate Chi2 reduction
    end
    
    %FIND TAUIT THAT MINIMIZES CHI2
    %determine direction of tauit search to minimize chi2
    [~, Chi2_tauit_plus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess+delta_tauit, C_guess); %guess value of Q for parameter combination with tauit+delta_tauit
    [~, Chi2_tauit_minus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess-delta_tauit, C_guess); %guess value of Q for parameter combination with tauit-delta_tauit
    sign_delta_tauit = -sign(Chi2_tauit_plus-Chi2_tauit_minus); %direction of search for value of tauit
    
    %increment tauit until Chi2 reaches minimum
    Chi2_previous_tauit = Chi2_guess; %record previous Chi2 for later comparison
    Chi2_guess = min([Chi2_tauit_plus,Chi2_tauit_minus]); %get new Chi2 guess from min of adjacent Chi2
    delta_Chi2_tauit = Chi2_guess - Chi2_previous_tauit; %calculate Chi2 reduction
    while delta_Chi2_tauit < 0
        Chi2_previous_tauit = Chi2_guess; %record previous Chi2 for later comparison
        tauit_guess = tauit_guess+delta_tauit*sign_delta_tauit; %increment tauit_guess
        [~, Chi2_guess] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess); %determine Chi2 for new tauit_guess
        delta_Chi2_tauit = Chi2_guess - Chi2_previous_tauit; %calculate Chi2 reduction
    end
    
    %FIND C THAT MINIMIZES CHI2
    %determine direction of C search to minimize chi2
    [~, Chi2_C_plus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess+delta_C); %guess value of Q for parameter combination with C+delta_C
    [~, Chi2_C_minus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess-delta_C); %guess value of Q for parameter combination with C-delta_C
    sign_delta_C = -sign(Chi2_C_plus-Chi2_C_minus); %direction of search for value of C
    
    %increment C until Chi2 reaches minimum
    Chi2_previous_C = Chi2_guess; %record previous Chi2 for later comparison
    Chi2_guess = min([Chi2_C_plus,Chi2_C_minus]); %get new Chi2 guess from min of adjacent Chi2
    delta_Chi2_C = Chi2_guess - Chi2_previous_C; %calculate Chi2 reduction
    while delta_Chi2_C < 0
        Chi2_previous_C = Chi2_guess; %record previous Chi2 for later comparison
        C_guess = C_guess+delta_C*sign_delta_C; %increment C_guess
        [~, Chi2_guess] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_guess, tauit_guess, C_guess); %determine Chi2 for new C_guess
        delta_Chi2_C = Chi2_guess - Chi2_previous_C; %calculate Chi2 reduction
    end
    
    %COMPUTE TOTAL REDUCTION IN CHI2 FOR THIS ITERATION
    delta_Chi2 = Chi2_guess - Chi2_previous;
    
    %COMPUTE REDUCTION IN CHI2 FOR TWO ITERATIONS AGO
    delta2_Chi2 = Chi2_guess - Chi2_previous2;
    
    %RECORD PARAMETER VALUES
    n_all = [n_all, n_guess];
    tauit_all = [tauit_all, tauit_guess];
    C_all = [C_all, C_guess];
    Chi2_all = [Chi2_all, Chi2_guess];
end

%SET FIT VALUES
%set final guess values as fit values
n_fit = n_guess;
tauit_fit = tauit_guess;
C_fit = C_guess;
Chi2_fit = Chi2_guess;

% %INCREMENT n AWAY FROM BEST FIT VALUE TO DETERMINE MAXIMUM n RANGE FOR delta_CHI2 <= 1
% delta_Chi2 = 0; %set initial Chi2 reduction as 0 to force at least one iteration of while loop 
% n_plus = n_fit; %set initial value for positive n search
% n_minus = n_fit; %set initial value for negative n search
% while delta_Chi2 <= delta_Chi2_uncertainty
%     n_plus = n_plus + delta_n; %increment value for positive n search
%     n_minus = n_minus - delta_n; %increment value for negative n search
%     
%     %COMPUTE VALUES OF TAUIT AND C TO MINIMIZE CHI2 WITH THESE N VALUES
%     [~, Chi2_n_plus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_plus, tauit_guess, C_guess); %get Chi2 for positive n increment
%     [~, Chi2_n_minus] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n_minus, tauit_guess, C_guess); %get Chi2 for negative n increment
% end

%PLOT FITTING PROCEDURE
iteration = 1:length(n_all);
figure;
subplot(5,1,1);
plot(iteration,Chi2_all-min(Chi2_all));
set(gca,'YScale','log');
ylabel('$${\chi}^2-min({\chi^2})$$','Interpreter','Latex');

subplot(5,1,2);
plot(iteration(2:end),-diff(Chi2_all));
set(gca,'yscale','log');
ylim([1e-3 1e3]);
ylabel('$$-{\Delta}{\chi}^2$$','Interpreter','Latex');

subplot(5,1,3);
plot(iteration,n_all);
ylabel('$$n$$','Interpreter','Latex');

subplot(5,1,4);
plot(iteration,tauit_all);
ylabel('$$\tau_{it}$$','Interpreter','Latex');

subplot(5,1,5);
plot(iteration,C_all);
ylabel('$$C$$','Interpreter','Latex');
xlabel('iteration');
