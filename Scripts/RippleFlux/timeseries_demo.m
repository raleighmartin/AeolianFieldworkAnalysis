for T_detrend = 1000:1000:5000;

    %T_detrend = 2000; %detrending time-scale in seconds
    T_smooth = 50; %smoothing time-scale in seconds

    %generate a synthetic timeseries
    t = 1:1:10000; %seconds
    z_raw = sin(t/500)+rand(1,10000)+t/2000;

    %perform detrending
    z_detrend = z_raw - moving_average(z_raw,t,T_detrend); %using moving average
    %z_detrend = detrend(z_raw); %using built-in "detrend"

    %perform smoothing
    z_smooth = moving_average(z_detrend,t,T_smooth);

    %plot things
    figure; clf;
    subplot(3,1,1);
    plot(t,z_raw);
    subplot(3,1,2);
    plot(t,z_detrend);
    subplot(3,1,3);
    plot(t,z_smooth);

end