%estimate for impact threshold from Lapotre et al. (2016), Eq. S1

function [ustit_ms] = ustit_Lapotreetal2016(d_m)

ustit_ms = exp((0.0408*log(d_m).^4)+...
    (1.2371*log(d_m).^3)+...
    (13.9837*log(d_m).^2)+...
    (70.3581*log(d_m))+...
    132.6241);