function [mu, sigma] = LogNormalGrainSize(d10,d50,d90)

%calculate mu and sigma
mu = log(d50);
sigma_d10 = (log(d10)-mu)/(sqrt(2)*erfinv(2*0.1-1)); %sigma based on d10
sigma_d90 = (log(d90)-mu)/(sqrt(2)*erfinv(2*0.9-1)); %sigma based on d90
sigma = mean([sigma_d10,sigma_d90]); %get average of sigmas
