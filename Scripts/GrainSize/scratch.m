Site_N = 3; %number of site

dates_all = [FluxBSNE_all{Site_N}.Date]';
dates_unique = unique(dates_all);
N_dates = length(dates_unique);
z_min = zeros(N_dates,1);
z_max = zeros(N_dates,1);

for i = 1:N_dates
    ind_date = find(dates_all == dates_unique(i));
    z_min_date = zeros(size(ind_date));
    z_max_date = zeros(size(ind_date));
    for j = 1:length(ind_date)
        z_date = FluxBSNE_all{Site_N}(ind_date(j)).z.z;
        z_min_date(j) = min(z_date(z_date>0));
        z_max_date(j) = max(z_date(z_date>0));
    end
    z_min(i) = min(z_min_date);
    z_max(i) = max(z_max_date);
end