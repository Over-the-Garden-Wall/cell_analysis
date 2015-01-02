function contact_type_by_angle(cell_num, target_types, no_angle_bins, distance_from_soma_threshold)
    C = get_constants;

    cell_dat = cell_data(cell_num);
    
    
    conts = double(cell_dat.contacts);
    
    num_types = length(target_types);
    type_num = zeros(1, size(conts,2));
    
    cell_mid = cell_dat.get_midpoint(true);
    for n = 1:num_types
        for k = 1:size(conts,2)
            d = sqrt(sum((conts(4:5,k)'-cell_mid(2:3)).^2));
            if d > distance_from_soma_threshold && any(conts(1,k)==C.type.(target_types{n})) 
               type_num(k) = n;
            end
        end
    end
    
    angle_histogram = zeros(no_angle_bins, num_types);
    
    for n = find(type_num);
        tn = type_num(n);
        cont_loc = conts(3:5,n);
        
        
        relative_loc = cont_loc - cell_mid';

        theta = atan2(relative_loc(3),relative_loc(2)) + pi;
        theta_bin = ceil(no_angle_bins * theta / 2 / pi);
        
        angle_histogram(theta_bin, tn) = angle_histogram(theta_bin, tn) + 1;
        
    end
    
    figure; plot(angle_histogram);
end
        
        
        