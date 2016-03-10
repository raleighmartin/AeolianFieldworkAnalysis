%% general script for linear fit

%% INPUTS
% x - independent variable values
% y - dependent variable values
% sigma_y - dependent variable uncertainty (optional, if not provided assume 0)

%% OUTPUTS
% yfit = a + b*x
% a - fitting intercept
% b - fitting slope
% sigma_a - uncertainty in a
% sigma_b - uncertainty in b
% sigma2_ab - covarying uncertainty in a and b
% yfit = predicted values of y based on linear fit with x-values as input
% sigma_yfit - corresponding uncertainty in predictions of y

%use method of bevington and robinson (p. 105)
function [a, b, sigma_a, sigma_b, yfit, sigma_yfit, sigma2_ab] = linearfit(x, y, sigma_y)

N = length(x); %number of observations
N_min = 2; %minimum number of observations for fit

%perform fit with uniform errors if none are provided
if (nargin==2 && N>=N_min)
	delta = N*sum(x.^2)-(sum(x))^2;
    a = (1/delta)*(sum(x.^2)*sum(y)-sum(x)*sum(x.*y));
    b = (1/delta)*(N*sum(x.*y)-sum(x)*sum(y));
    sigma_y = sqrt((1/(N-2))*sum((y-a-b*x).^2));
    sigma_a = sqrt(sigma_y^2*sum(x.^2)/delta);
    sigma_b = sqrt(N*sigma_y^2/delta);
	da_dy = (1/delta)*(sum(x.^2)-x*sum(x));
    db_dy = (1/delta)*(N*x-sum(x));
    sigma2_ab = sum(sigma_y.^2.*da_dy.*db_dy); %(Bevington and Robinson, Eq. 7.23)
	yfit = a+b*x;
	sigma_yfit = sqrt(sigma_a^2+sigma_b^2*x.^2+2*sigma2_ab*x); %uncertainty in prediction of y (Kok et al. 2014, Eq. A19)

%perform fit with variable errors if sufficient number of observations and nonzero erros
elseif (nargin==3 && N>=N_min) && (max(sigma_y)>0)
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

%otherwise perform standard polyfit and set errors to NaN
elseif N>=N_min
	P = polyfit(x,y,1); %perform linear fit
	a = P(2); b = P(1); %get fit parameters
	sigma_b = NaN;
	sigma_a = NaN;
    sigma2_ab = NaN;
	yfit = a+b*x;
	sigma_yfit = NaN*zeros(size(yfit));

%if not enough observations, set everything to NaN
else
    a = NaN;
    b = NaN;
    sigma_b = NaN;
    sigma_a = NaN;
    sigma2_ab = NaN;
    yfit = NaN;
    sigma_yfit = NaN;
end