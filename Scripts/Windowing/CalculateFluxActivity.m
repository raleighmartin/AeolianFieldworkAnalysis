%fD = detection frequency
%N_s = total counts rate second (across all Wenglors)
%dt_s = time increment (seconds)

function [fQ] = CalculateFluxActivity(fD,N_s,dt_s,fD_min,lambda)

%if not given, set assumed particle arrival rate according to fD
if nargin == 4
    lambda = N_s*dt_s/fD; %estimated particle passage rate during transport
end

%Determine flux frequency iteratively
if fD<=fD_min %if no (or negligible) flux, set frequencies to zero
    fQ = 0;
elseif fD==1 %if fD = 1, set fQ to 1
    fQ = 1;
else %otherwise, estimate fQ based on other parameters
    fQ = fD/(1-exp(-lambda)); %calculate fQ
end
