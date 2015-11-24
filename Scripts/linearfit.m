%% general script for linear fit

%% INPUTS
% x - independent variable values
% y - dependent variable values
% sigma_x - independent variable uncertainty (optional, if not provided
% assume 0)
% sigma_y - dependent variable uncertainty (optional, if not provided
% assume 0)

%% OUTPUTS
% yfit = a + b*x
% a - fitting intercept
% b - fitting slope
% sigma_a - uncertainty in a
% sigma_b - uncertainty in b
% yfit = predicted values of y based on linear fit with x-values as input
% sigma_yfit - corresponding uncertainty in predictions of y

%use method of bevington and robinson (p. 105)
function [a, b, sigma_a, sigma_b, yfit, sigma_yfit] = linearfit(x, y, sigma_x, sigma_y)

%set input error values to 0 if not given
if nargin == 2
    sigma_x = zeros(size(x));
    sigma_y = zeros(size(y));
end

%perform fit if sufficient number of observations
if length(x)>=3
    if (max(sigma_y)>0)&&(max(sigma_x)>0); %calculations if error is included
    %if (min(sigma_y)>0)&&(min(sigma_x)>0); %calculations if error is included
        delta = sum(1./sigma_y.^2).*sum(x.^2./sigma_y.^2)-(sum(x./sigma_y.^2)).^2; %(Bevington and Robinson, Eq. 6.12c)
        a = (1/delta)*(sum(x.^2./sigma_y.^2)*sum(y./sigma_y.^2)-sum(x./sigma_y.^2)*sum(x.*y./sigma_y.^2)); %(Bevington and Robinson, Eq. 6.12a)
        b = (1/delta)*(sum(1./sigma_y.^2)*sum(x.*y./sigma_y.^2)-sum(x./sigma_y.^2)*sum(y./sigma_y.^2)); %(Bevington and Robinson, Eq. 6.12b)
        sigma_a = sqrt((1/delta)*(sum(x.^2./sigma_y.^2))); %slope parameter
        sigma_b = sqrt((1/delta)*(sum(1./sigma_y.^2))); %intercept parameter
        da_dy = (1./(delta*sigma_y.^2)).*(sum(x.^2./sigma_y.^2)-x*sum(x./sigma_y.^2)); %(Bevington and Robinson, Eq. 6.20)
        db_dy = (1./(delta*sigma_y.^2)).*(x*sum(1./sigma_y.^2)-sum(x./sigma_y.^2)); %(Bevington and Robinson, Eq. 6.20)
        sigma2_ab = sum(sigma_y.^2.*da_dy.*db_dy); %(Bevington and Robinson, Eq. 7.23)
        yfit = a+b*x;
        sigma_yfit = sqrt(sigma_a^2+sigma_b^2*x.^2+2*sigma2_ab*x); %uncertainty in prediction of y (Kok et al. 2014, Eq. A19)
    else %calculations if no error included
        P = polyfit(x,y,1); %perform linear fit
        a = P(2); b = P(1); %get fit parameters
        sigma_b = NaN;
        sigma_a = NaN;
        yfit = a+b*x;
        sigma_yfit = NaN*zeros(size(yfit));
    end
else
    a = NaN;
    b = NaN;
    sigma_b = NaN;
    sigma_a = NaN;
    yfit = NaN;
    sigma_yfit = NaN;
end

% %% optional plot (comment out to hide this plot)
% figure(1); clf;
% errorbar(x,y,sigma_y,'bx'); hold on; %plot raw values with error bars
% [x_sort, ind_sort] = sort(x); %sort out x-values for confidence plot
% yfit_sort = yfit(ind_sort); %get corresponding sorting of yfit values
% sigma_yfit_sort = sigma_yfit(ind_sort); %get corresponding sorting of yfit uncertainties
% plot(x_sort,yfit_sort,'k'); %plot predicted values
% plot(x_sort,yfit_sort-sigma_yfit_sort,'k--',x_sort,yfit_sort+sigma_yfit_sort,'k--'); %plot confidence ranges
% xlabel('x');
% ylabel('y');
% set(gca,'FontSize',16);