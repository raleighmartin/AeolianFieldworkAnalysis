%given volume fractions by bin with range d_lower - d_upper, compute d10,
%d50, and d90

function [d10, d50, d90] = ReferenceGrainSizes(dV, d_lower, d_upper)

%get cumulative distribution
CV = cumsum(dV);
CV = CV+(1-max(CV)); %if numbers don't add to 100%, adjust so that top value is 1

ind_d10 = find(CV>=0.1, 1);
d_below = d_lower(ind_d10);
d_above = d_upper(ind_d10);
wt_below = (CV(ind_d10)-0.1)/dV(ind_d10);
wt_above = 1-wt_below;
d10 = exp(wt_below*log(d_below)+wt_above*log(d_above));

%get d50
ind_d50 = find(CV>=0.5, 1);
d_below = d_lower(ind_d50);
d_above = d_upper(ind_d50);
wt_below = (CV(ind_d50)-0.5)/dV(ind_d50);
wt_above = 1-wt_below;
d50 = exp(wt_below*log(d_below)+wt_above*log(d_above));

%get d90
ind_d90 = find(CV>=0.9, 1);
d_below = d_lower(ind_d90);
d_above = d_upper(ind_d90);
wt_below = (CV(ind_d90)-0.9)/dV(ind_d90);
wt_above = 1-wt_below;
d90 = exp(wt_below*log(d_below)+wt_above*log(d_above));