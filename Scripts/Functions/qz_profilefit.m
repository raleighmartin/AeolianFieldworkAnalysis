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

%set input error values arbitrarily as one if input is not given
if nargin == 2
    sigma_qz = ones(size(qz)); 
    sigma_z = ones(size(z));
end

%combine values that are at the same height
z_unique = unique(z);
[n, ~] = histc(z, z_unique);
ind_repeat = find(n>1);
ind_norep = find(n==1);
if ~isempty(ind_repeat)
    z_new = z_unique;
    qz_new = zeros(size(z_new));
    sigma_z_new = zeros(size(z_new));
    sigma_qz_new = zeros(size(z_new));
    for i = 1:length(ind_norep)
        ind = ind_norep(i);
        qz_new(ind) = qz(z==z_new(ind));
        sigma_z_new(ind) = sigma_z(z==z_new(ind));
        sigma_qz_new(ind) = sigma_qz(z==z_new(ind));
    end
    for i = 1:length(ind_repeat)
        ind = ind_repeat(i);
        qz_new(ind) = mean(qz(z==z_new(ind)));
        sigma_z_new(ind) = mean(sigma_z(z==z_new(ind)));
        sigma_qz_new(ind) = (max(qz(z==z_new(ind))+sigma_qz(z==z_new(ind)))-...
            min(qz(z==z_new(ind))-sigma_qz(z==z_new(ind))))...
            /sqrt(length(find(z==z_new(ind))));
    end
    %rename values
    z_old = z;
    z = z_new;
    qz = qz_new;
    sigma_z = sigma_z_new;
    sigma_qz = sigma_qz_new;
end

%compute log of flux and associated error
logqz = log(qz); %[log(q)]
sigma_logqz = sqrt((sigma_qz.^2./qz.^2)+(sigma_z.^2./z.^2)); %[log(q)] (combines height and flux error)

%perform fit if sufficient number of observations
if length(z)>=3
    if (max(sigma_qz)>0)&&(max(sigma_z)>0); %calculations if error is nonzero
        %call linearfit.m to get fitting parameters
        if nargin == 2; %if errors were not provided, exclude them in call to linearfit
            [a, b, sigma_a, sigma_b, ~, sigma_logqz_fit] = linearfit(z,logqz);
        else %otherwise, use errors
            [a, b, sigma_a, sigma_b, ~, sigma_logqz_fit] = linearfit(z,logqz,sigma_z,sigma_logqz);
        end
        q0 = exp(a); %[g/m^2/s]
        zq = -1/b; %[m]
        sigma_q0 = sigma_a*q0; %[g/m^2/s]
        sigma_zq = sigma_b*zq.^2; %[m]
        qz_fit = q0*exp(-z/zq); %prediction of qz [g/m^2/s] from least squares fit
        sigma_qz_fit = sigma_logqz_fit.*qz_fit;
    else %calculations if error is zero
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

%update "fit" values if there are repeats
if ~isempty(ind_repeat)
    %initialize fit values
    qz_fit_new = zeros(size(z_old));
    sigma_qz_fit_new = zeros(size(z_old));
    sigma_logqz_fit_new = zeros(size(z_old));
    
    %go through and add fit values for repeat entries
    for i = 1:length(z)
        qz_fit_new(z_old==z(i)) = qz_fit(i);
        sigma_qz_fit_new(z_old==z(i)) = sigma_qz_fit(i);
        sigma_logqz_fit_new(z_old==z(i)) = sigma_logqz_fit(i);
    end
    
    %rename values
    qz_fit = qz_fit_new;
    sigma_qz_fit = sigma_qz_fit_new;
    sigma_logqz_fit = sigma_logqz_fit_new;
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