% dV = volume fraction by bin
% d_distribution_lower = diameters for lower end of each grain size volume fraction
% d_distribution_upper = diameters for upper end of each grain size volume fraction
% d_bin_lower = lower limit of imposed grain size bins
% d_bin_upper = upper limit of imposed grain size bins

function f_d = GrainSizeBinFraction(dV, d_distribution_lower, d_distribution_upper, d_bin_lower, d_bin_upper)

N_d_bin = length(d_bin_lower); %get number of bins

f_d = zeros(1,N_d_bin); %initialize list of fractions by bin

for k = 1:N_d_bin
    ind_full_bin = find(d_distribution_lower>=d_bin_lower(k) & d_distribution_upper<=d_bin_upper(k)); %grain size bins completely in bin
    if ~isempty(ind_full_bin)
        ind_bottom_bin = min(ind_full_bin)-1; %lower (partial) bin
        ind_top_bin = max(ind_full_bin)+1; %upper (partial) bin
    else
        ind_bottom_bin = find(d_distribution_upper>d_bin_lower(k), 1 ); %lower (partial) bin
        ind_top_bin = find(d_distribution_lower<d_bin_upper(k), 1, 'last' ); %upper (partial) bin
    end
    wt_bottom_bin = (log(d_distribution_upper(ind_bottom_bin))-log(d_bin_lower(k)))/...
        (log(d_distribution_upper(ind_bottom_bin))-log(d_distribution_lower(ind_bottom_bin))); %fraction of partial bin for this bin - log scale
    wt_top_bin = (log(d_bin_upper(k))-log(d_distribution_lower(ind_top_bin)))/...
        (log(d_distribution_upper(ind_top_bin))-log(d_distribution_lower(ind_top_bin))); %fraction of partial bin for this bin - log scale
    f_d(k) = sum(dV(ind_full_bin))+...
            wt_bottom_bin*dV(ind_bottom_bin)+...
            wt_top_bin*dV(ind_top_bin); %compute fraction of grains by volume in this bin
end