function data = analyze_skeleton(cell_name)

    res = [16.5 16.5 25];
    allowed_outlier_percentage = .0001;

    coords = get_coords_for_sac(cell_name);
    for k = 1:3
        coords(:,k) = coords(:,k)*res(k);
    end
    
    [Q P valid_coords] = find_planar_rotation_iterative(coords, allowed_outlier_percentage);    

    rotated_valid_coords = valid_coords*Q';
    rotated_coords = coords*Q';
        
    data.P = P;
    data.valid_mean_depth = mean(rotated_valid_coords(:,1));
    data.mean_depth = mean(rotated_coords(:,1));
    data.x_spread = std(rotated_coords(:,1));    
    data.y_spread = std(rotated_coords(:,2));
    data.z_spread = std(rotated_coords(:,3));
    data.coords = coords;
    data.rotated_coords = rotated_coords;
    
    
    
end