%% INPUTS
% X is timeseries vector
% dt is time increment in data vector

%% OUTPUTS
% R is autocorrelation values
% lags are times associated with autocorrelation values

function [R lags] = autocorrelation(X,dt,maxlag)

%set maxlag (# of time steps) if not specified by user
if nargin == 2
  maxlag = 500;
end

%initialize output variables
R = zeros(maxlag+1,1);
lags = (0:maxlag)*dt;

%compute autocorrelation for each lag
Xbar = mean(X);
Xvar = var(X);
for i = 0:maxlag
    X0 = X(1:(end-i))-Xbar;
    X1 = X((1+i):end)-Xbar;
    R(i+1) = mean(X0.*X1)/Xvar;
end