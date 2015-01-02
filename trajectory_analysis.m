function trajectory_analysis(cell_num, target_types, no_trajectory_bins, cont_neighborhood_size, distance_from_soma_threshold)
    C = get_constants;

    cell_dat = cell_data(cell_num);
    p = cell_dat.get_surface;
    
    
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
    
    trajectory_histogram = zeros(no_trajectory_bins, num_types);
    
    for n = find(type_num);
        tn = type_num(n);
        cont_loc = conts(3:5,n);
        d = sqrt((p(:,1)-cont_loc(1)).^2 + (p(:,2)-cont_loc(2)).^2 + (p(:,3)-cont_loc(3)).^2);
        rel_p = p(d <= cont_neighborhood_size,:);
        
        cov_mat = cov(rel_p);
        try
            
            [V,D]=eigs(cov_mat,1);

            relative_loc = cont_loc - cell_mid';

            phi = asin(V(1));

            if relative_loc(2)*V(2) + relative_loc(3)*V(3) < 0
                phi = -phi;
            end

            phi = phi+pi/2;
            phi_bin = ceil(phi/pi * no_trajectory_bins);

            trajectory_histogram(phi_bin, tn) = trajectory_histogram(phi_bin, tn) + conts(2,n);
        catch ME
            disp(ME.message)
        end
    end
    
    figure; plot(trajectory_histogram);
end
        
        
        