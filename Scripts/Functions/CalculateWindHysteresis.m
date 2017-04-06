%% determine fractions of u in different hysteresis regions

%INPUTS
% u = wind velocity time series
% u_ft = fluid threshold wind speed
% u_it = impact threshold wind speed
% init_state = describes most recent state outside hysteresis zone
% -1 if most recent time was below it
% +1 if most recent time was above ft
% 0 if most recent time unknown

%OUTPUTS
% fplus = fraction of time above fluid threshold
% fminus = fraction of time below impact threshold
% fint = fraction of time between thresholds
% fint_up = fraction of time between thresholds from above
% fint_down = fraction of time between thresholds from below
% final_state = final state of wind for increment (-1 if fminus or fint_up, +1 if fplus or fint_down)

%% function script
function [fplus,fminus,fint,fint_up,fint_down,...
    ind_fplus,ind_fminus,ind_fint,ind_fint_up,ind_fint_down,final_state] = CalculateWindHysteresis(u,u_ft,u_it,init_state)

if nargin==3
    init_state = 0;
end

%get number of time steps
T = length(u);

%get fplus, fminus, and fint
ind_fplus = find(u>=u_ft); %indices for u above u_ft
ind_fminus = find(u<u_it); %indices for u below u_it
ind_fint = intersect(find(u>=u_it),find(u<u_ft)); %indices for u in intermediate zone
N_fplus = length(ind_fplus); %number of u above u_ft
N_fminus = length(ind_fminus); %number of u below u_it
N_fint = length(ind_fint); %number of u in intermediate zone
fplus = N_fplus/T; %fraction of time with u above fluid threshold
fminus = N_fminus/T; %fraction of time with u below impact threshold
fint = N_fint/T; %fraction of time with u in intermediate zone

%look at hysteresis for u in intermediate zone
N_fint_up = 0; %initialize N_fint_up (number of upcrossings to hysteresis zone)
N_fint_down = 0; %initialize N_fint_down (number of downcrossings to hysteresis zone)
N_fint_unknown = 0; %initialize N_fint_isempty (number of events in hysteresis zone with no known history)
ind_fint_up = []; %initialize indices of int values from above
ind_fint_down = []; %initialize indices of int values from below
ind_fint_unknown = []; %initialize indices of int values from unknown state
for l = 1:N_fint
    ind_last_fplus = ind_fplus(find(ind_fplus<ind_fint(l), 1, 'last' ));
    ind_last_fminus = ind_fminus(find(ind_fminus<ind_fint(l), 1, 'last' ));
    if(isempty(ind_last_fplus))
        if(isempty(ind_last_fminus))
            if init_state == -1
                N_fint_up = N_fint_up+1;
                ind_fint_up = [ind_fint_up; ind_fint(l)];
            elseif init_state == 1
                N_fint_down = N_fint_down+1;
                ind_fint_down = [ind_fint_down; ind_fint(l)];
            else
                N_fint_unknown = N_fint_unknown+1;
                ind_fint_unknown = [ind_fint_unknown; ind_fint(l)];
            end
        else
            N_fint_up = N_fint_up+1;
            ind_fint_up = [ind_fint_up; ind_fint(l)];
        end
    elseif(isempty(ind_last_fminus))
        if(isempty(ind_last_fplus))
            if init_state == -1
                N_fint_up = N_fint_up+1;
                ind_fint_up = [ind_fint_up; ind_fint(l)];
            elseif init_state == 1
                N_fint_down = N_fint_down+1;
                ind_fint_down = [ind_fint_down; ind_fint(l)];
            else
                N_fint_unknown = N_fint_unknown+1;
                ind_fint_unknown = [ind_fint_unknown; ind_fint(l)];
            end
        else
            N_fint_down = N_fint_down+1;
            ind_fint_down = [ind_fint_down; ind_fint(l)];
        end
    else
        if(ind_last_fminus>ind_last_fplus)
            N_fint_up = N_fint_up+1;
            ind_fint_up = [ind_fint_up; ind_fint(l)];
        elseif(ind_last_fplus>ind_last_fminus)
            N_fint_down = N_fint_down+1;
            ind_fint_down = [ind_fint_down; ind_fint(l)];
        end
    end
end

%tabulate fractions of up and down hysteresis events
N_fint_known = N_fint_up+N_fint_down;
if N_fint_known == 0
    fint_up = 0;
    fint_down = 0;
else
    fint_up = (N_fint_up/N_fint_known)*fint;
    fint_down = (N_fint_down/N_fint_known)*fint;
end

%set final state
ind_transport_max = max([max(ind_fplus),max(ind_fint_down),0]); %maximum index of transport state
ind_nontransport_max = max([max(ind_fminus),max(ind_fint_up),0]); %maximum index of non-transport state

if ind_transport_max > ind_nontransport_max
    final_state = 1;
elseif ind_nontransport_max > ind_transport_max
    final_state = -1;
end