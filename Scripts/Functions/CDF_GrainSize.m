%given volume fractions by bin with range d_lower - d_upper, compute CDF value for grain size

%INPUTS
%dV = fraction of grain size distribution (by volume) in each bin
%d_lower = lower limit of each bin
%d_upper = upper limit of each bin
%d_CDF = reference grain size associated with f_CDF

%OUTPUTS
%f_CDF = fractional value in CDF for reference grain size

function [f_CDF] = CDF_GrainSize(dV, d_lower, d_upper, d_CDF)

%get cumulative distribution
CV = cumsum(dV);
CV = CV+(1-max(CV)); %if numbers don't add to 100%, adjust so that top value is 1

%determine f_CDF
ind_d_CDF = find(d_CDF>=d_lower&d_CDF<=d_upper);
d_below = d_lower(ind_d_CDF);
d_above = d_upper(ind_d_CDF);
wt_partialbin = (log10(d_CDF)-log10(d_below))/(log10(d_above)-log10(d_below));
f_CDF = CV(ind_d_CDF-1)+dV(ind_d_CDF)*wt_partialbin;