function [mu sigma d10 d50 d90] = LogNormalGrainSize(d_lower, d_upper, dV)

%cumulative distribution by volume
x = cumsum(dV); 

%get d10
ind_d10 = find(x>=0.1, 1);
d_below = d_lower(ind_d10);
d_above = d_upper(ind_d10);
wt_below = (x(ind_d10)-0.1)/dV(ind_d10);
wt_above = 1-wt_below;
d10 = exp(wt_below*log(d_below)+wt_above*log(d_above));

%get d50
ind_d50 = find(x>=0.5, 1);
d_below = d_lower(ind_d50);
d_above = d_upper(ind_d50);
wt_below = (x(ind_d50)-0.5)/dV(ind_d50);
wt_above = 1-wt_below;
d50 = exp(wt_below*log(d_below)+wt_above*log(d_above));

%get d90
ind_d90 = find(x>=0.9, 1);
d_below = d_lower(ind_d90);
d_above = d_upper(ind_d90);
wt_below = (x(ind_d90)-0.9)/dV(ind_d90);
wt_above = 1-wt_below;
d90 = exp(wt_below*log(d_below)+wt_above*log(d_above));

%get centers of bins
d = sqrt(d_lower.*d_upper);

%calculate mu and sigma
mu = log(d50);
sigma_d10 = (log(d10)-mu)/(sqrt(2)*erfinv(2*0.1-1));
sigma_d90 = (log(d90)-mu)/(sqrt(2)*erfinv(2*0.9-1));
sigma = mean([sigma_d10,sigma_d90]);

figure;
dVdlogd = dV./log(d_upper./d_lower);
plot(d,dV); hold on;
plot(d,lognpdf(d,mu,sigma)./sum(lognpdf(d,mu,sigma)));
set(gca,'xscale','log');
xlabel('d');
ylabel('dVdlogd');
