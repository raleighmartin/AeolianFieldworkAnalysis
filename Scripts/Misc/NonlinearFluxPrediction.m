function [Q_pred, Chi2_pred] = NonlinearFluxPrediction(nu, ust, tau, Q, sigma_tau, sigma_Q, n, tauit, C)

Q_pred = C*ust.^n.*(tau-tauit); %predict value of Q for parameter combination
sigma_Q_tau = C*ust.^n.*sigma_tau; %compute uncertainty in Q due to tau (ignore uncertainty due to u*, since it is correlated with tau)
sigma_Q_total = sqrt(sigma_Q.^2+sigma_Q_tau.^2); %compute total uncertainty in Q as combination of Q and tau uncertainty
Chi2_pred = Chi2Calculation(Q, sigma_Q_total, Q_pred, nu); %compute Chi2 for guess values