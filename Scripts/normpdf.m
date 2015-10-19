%function to compute pdf value(s) associated with vector X based on
%Gaussian distribution with mean "mu" and standard deviation, "sigma"

function Y = normpdf(X,mu,sigma)

Y = (1/sigma*sqrt(2*pi))*exp(-(X-mu).^2/(2*sigma^2));

