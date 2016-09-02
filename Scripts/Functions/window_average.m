%t_out gives timeseries of first time for each window

function [data_out, t_out] = window_average(data_in, t_in, dt_avg)

L_data = length(data_in);
dt = mode(diff(t_in));

%save some time by checking if dt = dt_avg
if dt == dt_avg
    data_out = data_in;
    t_out = t_in;
else
    L_chunk = dt_avg/dt;
    N_chunk = floor(L_data/L_chunk);

    data_out = zeros(N_chunk,1);
    for n=1:N_chunk
        ind_start = floor(L_chunk*(n-1)+1); %first index for window average
        ind_end = floor(L_chunk*n); %last index for window average
        data_out(n) = mean(data_in(ind_start:ind_end)); %get window average
    end

    t_out = (t_in(1)+(0:(N_chunk-1))*dt_avg)';
    
%     L_chunk = round(dt_avg/dt);
%     N_chunk = floor(L_data/L_chunk);
% 
%     data_out = zeros(N_chunk,1);
%     for n=1:N_chunk
%         data_out(n) = mean(data_in(L_chunk*(n-1)+(1:L_chunk)));
%     end
% 
%     t_out = (t_in(1)+(0:(N_chunk-1))*dt_avg)';
end