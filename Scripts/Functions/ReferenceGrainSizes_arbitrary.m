%given volume fractions by bin with range d_lower - d_upper, compute
%arbitrary CDF values

%INPUTS
%dV = fraction of grain size distribution (by volume) in each bin
%d_lower = lower limit of each bin
%d_upper = upper limit of each bin
%f_CDF = fractional value in CDF for reference grain size

%OUTPUTS
%d_CDF = reference grain size associated with f_CDF

function [d_CDF] = ReferenceGrainSizes_arbitrary(dV, d_lower, d_upper, f_CDF)

%get cumulative distribution
CV = cumsum(dV);
CV = CV+(1-max(CV)); %if numbers don't add to 100%, adjust so that top value is 1

%determine d_CDF
ind_d_CDF = find(CV>=f_CDF, 1);
d_below = d_lower(ind_d_CDF);
d_above = d_upper(ind_d_CDF);
wt_below = (CV(ind_d_CDF)-f_CDF)/dV(ind_d_CDF);
wt_above = 1-wt_below;
d_CDF = exp(wt_below*log(d_below)+wt_above*log(d_above));