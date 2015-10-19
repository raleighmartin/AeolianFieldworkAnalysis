%function to compute pdf value(s) associated with vector X based on
%Gaussian distribution with mean "mu" and standard deviation, "sigma"

function Y = normcdf(X,mu,sigma)

Y = (1/2)*(1+erf((X-mu)/(sigma*sqrt(2*pi))));