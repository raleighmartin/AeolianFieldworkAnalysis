%ntotal = total counts rate per second
%deltat_s = sampling interval time in seconds
%fD_min = mininum detection frequency (set to 0 if below this value)

function [fD,fQ,lambda] = CalculateFluxActivity(ntotal,deltat_s,fD_min)

ntotal_bar = mean(ntotal); %mean total counts per second
T = length(ntotal); %number of data points
fD = sum(ntotal>0)/T; %detection rate

%estimate particle arrival rate per sampling interval
if fD<=fD_min %set to zero if below detection limit
    lambda = 0;
else %otherwise estimate arrival rate per sampling interval based on fD
    lambda = ntotal_bar*deltat_s/fD;
end

%estimate flux activity from flux detection rate and particle arrival rate per sampling interval
if lambda==0 %if no (or negligible) flux, set frequencies to zero
    fQ = 0;
elseif fD==1 %if fD = 1, set fQ to 1
    fQ = 1;
else %otherwise, estimate fQ based on other parameters
    fQ = fD/(1-exp(-lambda)); %calculate fQ
    if fQ>1
        fQ = 1; %if correction gives fQ>1, just set as 1
    end
end