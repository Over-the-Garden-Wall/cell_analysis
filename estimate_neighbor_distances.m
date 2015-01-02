function neighbor_distances = estimate_neighbor_distances(cell_nums)

    num_cells = length(cell_nums);
    
    dist = ones(num_cells)*Inf;
    is_neighbor = false(num_cells);
    
    midpoints = zeros(num_cells,3);
    hulls = cell(num_cells,1);
    
    for k = 1:num_cells
        cell_dat = cell_data(cell_nums(k));
        midpoints(k,:) = cell_dat.get_midpoint(false);
        hulls{k} = cell_dat.hull_2d;
    end
    
    
    for k = 1:num_cells
        
        for l = k+1:num_cells

            
            dist(k,l) = sqrt(sum((midpoints(k,2:3)-midpoints(l,2:3)).^2));
            
            is_neighbor(k,l) = ~isempty(polybool('intersection',hulls{k}(:,1),hulls{k}(:,2),hulls{l}(:,1),hulls{l}(:,2)));
            
        end
    end
    
    max_good_dist = max(dist(is_neighbor(:)));
    
%     is_neighbor = dist <= max_good_dist;
    
    neighbor_distances = dist(is_neighbor(:));
end