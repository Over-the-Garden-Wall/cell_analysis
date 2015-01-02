function histogram_by_depth_and_angle(cell_num, no_angle_bins, depth_bin_size)
    C = get_constants;

    cell_dat = cell_data(cell_num);
    
    p = cell_dat.get_surface;
    d = C.f(p(:,1));
    
    depth_bins = min(d):depth_bin_size:max(d);
    no_depth_bins = length(depth_bins);
    
    angle_bins = (0:no_angle_bins-1)/no_angle_bins * 2 * pi;
    angle_bin_size = angle_bins(2)-angle_bins(1);
    
    cell_center = cell_dat.get_midpoint(true);
    for n = 2:3
        p(:,n) = p(:,n) - cell_center(n);
    end
    
    theta = atan2(p(:,3),p(:,2)) + pi;
    
    histogram_values = zeros(no_depth_bins, no_angle_bins);
    
    ab = 1+floor(theta / 2 / pi * no_angle_bins);
    db = 1+floor((d-min(d)) / depth_bin_size);
    
    for n = 1:size(p,1)
        
        histogram_values(db(n), ab(n)) = histogram_values(db(n), ab(n))+1;
    
    end
    
    for n = 1:size(histogram_values,1);
        histogram_values(n, :) = histogram_values(n, :) / sqrt(sum(histogram_values(n, :).^2));
    end
    
    
    figure;
    imagesc(histogram_values);
    
end
    
    