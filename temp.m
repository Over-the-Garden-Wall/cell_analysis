C = get_constants;

bins = C.j_bins;
bin_size = bins(2)-bins(1);
pseudo_bins = [bins(1)-bin_size bins bins(end)+bin_size];

total_hist = zeros(1,length(bins));

for c = C.type.j
    cell_dat = cell_data(c);
    p = cell_dat.get_surface;
    cell_mid = cell_dat.get_midpoint(true);
    p(:,2) = p(:,2) - cell_mid(2);
    p(:,3) = p(:,3) - cell_mid(3);
    d = p(:,2)*cell_dat.dist_axis(1) + p(:,3)*cell_dat.dist_axis(2);
    d = d/1000;
    my_hist = hist(d,pseudo_bins);
    total_hist = total_hist + my_hist(2:end-1);
end

figure; plot(bins, total_hist, 'LineWidth', 2);

prep_figure(gcf,gca,'xlabel', 'Distance from Soma (microns)', 'ylabel', 'Surface Area (surface voxels)');
