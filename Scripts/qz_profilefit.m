%% function to compute saltation flux profile assuming exponential

%%INPUTS
%z - observation heights (m)
%qz - vertical fluxes (g/m^2/s)
%sigma_z = height uncertainty
%sigma_qz = flux uncertainty

%%OUTPUTS
%q0 - profile fit (g/m^2/s)
%ze - profile fit e-folding height (m)
%qz_fit - predicted q from profile fit
%sigma_q0 - uncertainty in q0
%sigma_ze - uncertainty in ze
%sigma_qz_fit - uncertainty in predictions of qz (varies with qz)
%sigma_logqz_fit - uncertainty in predictions of log(qz)

%% FITTING TO:
% q = q0*exp(-z/zq)
% logq = a+b*z
% q0 = exp(a)
% zq = -1/b

%use method of bevington and robinson (p. 105)
function [q0,zq,sigma_q0,sigma_zq,qz_fit,sigma_qz_fit,sigma_logqz_fit] = qz_profilefit(qz, z, sigma_qz, sigma_z)

%set input error values to 0 if not given
if nargin == 2
    sigma_qz = 0;
    sigma_z = 0;
end

%compute log of flux and associated error
logqz = log(qz); %[log(q)]
sigma_logqz = sqrt((sigma_qz.^2./qz.^2)+(sigma_z.^2./z.^2)); %[log(q)] (combines height and flux error)

%perform fit if sufficient number of observations
if length(z)>=3
    if (max(sigma_qz)>0)&&(max(sigma_z)>0); %calculations if error is included
        delta = sum(1./sigma_logqz.^2).*sum(z.^2./sigma_logqz.^2)-(sum(z./sigma_logqz.^2)).^2; %[m^2/log(q)^4] (Bevington and Robinson, Eq. 6.12c)
        a = (1/delta)*(sum(z.^2./sigma_logqz.^2)*sum(logqz./sigma_logqz.^2)-sum(z./sigma_logqz.^2)*sum(z.*logqz./sigma_logqz.^2)); %[log(q)] (Bevington and Robinson, Eq. 6.12a)
        b = (1/delta)*(sum(1./sigma_logqz.^2)*sum(z.*logqz./sigma_logqz.^2)-sum(z./sigma_logqz.^2)*sum(logqz./sigma_logqz.^2)); %[1/m] (Bevington and Robinson, Eq. 6.12b)
        sigma_a = sqrt((1/delta)*(sum(z.^2./sigma_logqz.^2))); %[log(q)]
        sigma_b = sqrt((1/delta)*(sum(1./sigma_logqz.^2))); %[1/m]
        da_dlogqz = (1./(delta*sigma_logqz.^2)).*(sum(z.^2./sigma_logqz.^2)-z*sum(z./sigma_logqz.^2)); %[] (Bevington and Robinson, Eq. 6.20)
        db_dlogqz = (1./(delta*sigma_logqz.^2)).*(z*sum(1./sigma_logqz.^2)-sum(z./sigma_logqz.^2)); %[1/m] (Bevington and Robinson, Eq. 6.20)
        sigma2_ab = sum(sigma_logqz.^2.*da_dlogqz.*db_dlogqz); %[log(q)^2/m] (Bevington and Robinson, Eq. 7.23)
        q0 = exp(a); %[g/m^2/s]
        zq = -1/b; %[m]
        sigma_q0 = sigma_a*q0; %[g/m^2/s]
        sigma_zq = sigma_b*zq.^2; %[m]
        qz_fit = q0*exp(-z/zq); %prediction of qz [g/m^2/s] from least squares fit
        sigma_logqz_fit = sqrt(sigma_a^2+sigma_b^2*z.^2+2*sigma2_ab*z); %[log(q)] uncertainty in prediction of log(qz) (Kok et al. 2014, Eq. A19)
        sigma_qz_fit = sigma_logqz_fit.*qz_fit;
        %sigma_logqz_fit = sqrt((1/(length(z)-2))*sum((logqz-a-b*z).^2)); %[log(q)], Bevington and Robinson (2003) Eq. 6.15
        %sigma_qz_fit = sigma_logqz_fit*qz_fit;
    else %calculations if no error included
        P = polyfit(z,logqz,1); %perform linear fit
        a = P(2); b = P(1); %get fit parameters
        q0 = exp(a); %[g/m^2/s]
        zq = -1/b; %[m]
        sigma_q0 = NaN;
        sigma_zq = NaN;
        qz_fit = q0*exp(-z/zq);
        sigma_qz_fit = NaN;
        sigma_logqz_fit = NaN;
    end
else
    q0 = NaN;
    zq = NaN;
    sigma_q0 = NaN;
    sigma_zq = NaN;
    qz_fit = NaN;
    sigma_qz_fit = NaN;
    sigma_logqz_fit = NaN;
end

% %optional plot
% [z_sort ind_sort] = sort(z);
% qz_fit_sort = qz_fit(ind_sort);
% sigma_qz_fit_sort = sigma_qz_fit(ind_sort);
% figure(1); clf;
% errorbar(z,qz,sigma_logqz.*qz,'bx'); hold on; %plot raw values with error bars
% plot(z_sort,qz_fit_sort,'k'); %plot predicted values
% plot(z_sort,qz_fit_sort-sigma_qz_fit_sort,'k--',z_sort,qz_fit_sort+sigma_qz_fit_sort,'k--'); %plot confidence ranges
% xlabel('z (m)');
% ylabel('q (g/m^s/s)');
% set(gca,'yscale','log');
% set(gca,'FontSize',16);
% pause;