%fD = detection frequency
%n_s = average counts rate profile (number per second)
%dt_s = time increment (seconds)

function [fQ] = CalculateFluxActivity(fD,n_s,dt_s,fD_min)

%set minimum detection rate if not given
if nargin == 3
    fD_min = 0.005; %minimum detection rate for zero flux
end

%Determine flux frequency iteratively
if fD<=fD_min %if no (or negligible) flux, set frequencies to zero
    fQ = 0;
elseif fD==1 %if fD = 1, set fQ to 1
    fQ = 1;
else %otherwise, do calculations, estimating fQ iteratively 
    fQ_prev = 0; %previous estimate of fQ (arbitrarily set to 0 for while loop functioning) 
    fQ = fD; %start by assuming fQ=fD
    while abs(fQ-fQ_prev)>0.001 %compare current and previous estimations of fQ for while loop
        fQ_prev = fQ; %update fQ_prev with last fQ
        lambda = sum(n_s)*dt_s/fQ; %estimated particle passage rate during transport - correcting error and putting dt in numerator
        %lambda = sum(n_s)/(fQ*dt_s); %estimated particle passage rate during transport
        fQ = fD/(1-exp(-lambda)); %calculate new fQ
        
        %if fQ>1, set it as NaN and cut out of cycle
        if fQ>1
            fQ = NaN;
            break;
        end
    end

    %remove points with fQ>1
    if fQ>1
        fQ = NaN;
    end
end
% else %otherwise estimate fQ from fD and mean particle passage rate during transport
%     lambda = sum(n_s)*dt_s/fD; %estimated particle passage rate during transport - now correcting error and putting dt in numerator
%     %lambda = sum(n_s)/(fD*dt_s); %estimated particle passage rate during transport - now using fD
%     fQ = fD/(1-exp(-lambda)); %calculate fQ
% end
