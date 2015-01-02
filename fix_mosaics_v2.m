function final_clusters = fix_mosaics(initial_clusters, improvement_threshold)

    C = get_constants;

    num_clusters = length(initial_clusters);
    all_cell_nums = [];
    for k = 1:num_clusters
        all_cell_nums = [all_cell_nums initial_clusters{k}];
    end
    total_cells = length(all_cell_nums);

    perc10 = zeros(total_cells,1);
    cell_volume = zeros(total_cells,1);
    hulls = cell(total_cells,1);
    neighbors = cell(total_cells,1);
    centers = zeros(total_cells,2);
    radii = zeros(total_cells,1);
    current_membership = zeros(total_cells,1);
    cell_area = zeros(total_cells,1);
    dist_to_edge = Inf*ones(total_cells,1);
    nearest_edge = zeros(total_cells,1);
    
    perc10_edge = 23;
    volume_edge = 9.6*10^5;
    
    for c = 1:total_cells;
        cell_dat = cell_data(all_cell_nums(c));
        p = cell_dat.get_surface;
        d = C.f(p(:,1));
        d(d<0) = [];
        d = sort(d);
        
        perc10(c) = d(ceil(length(d)*.1));
        
        cell_volume(c) = length(d);
        [hulls{c}(:,1), hulls{c}(:,2)] = poly2cw(cell_dat.hull_2d(:,1), cell_dat.hull_2d(:,2));
        cm = cell_dat.get_midpoint(false);
        centers(c,:) = cm(2:3);
        radii(c) = max(sqrt((cm(2)-hulls{c}(:,1)).^2 + (cm(3)-hulls{c}(:,2)).^2));
        cell_area(c) = poly_area(hulls{c});
        
        for k = 1:num_clusters
            if any(all_cell_nums(c)==initial_clusters{k})
                current_membership(c) = k;
            end
        end
        
    end
    
    initial_membership = current_membership;
    
    vol_min = min(cell_volume);
    vol_max = max(cell_volume);
    
    perc_min = min(perc10);
    perc_max = max(perc10);
    
    perc10 = (perc10-perc_min)/(perc_max-perc_min);
    perc10_edge = (perc10_edge-perc_min)/(perc_max-perc_min);
    
    cell_volume = (cell_volume-vol_min)/(vol_max-vol_min);
    volume_edge = (volume_edge-vol_min)/(vol_max-vol_min);
    
    for c = 1:total_cells
        if current_membership(c)==1
            nearest_edge(c) = (cell_volume(c) > volume_edge) + 2;
            dist_to_edge(c) = abs(perc10(c)-perc10_edge);
        elseif current_membership(c)==2
            dist_to_1edge = abs(perc10(c)-perc10_edge);
            dist_to_3edge = abs(cell_volume(c)-volume_edge);
            if dist_to_1edge<dist_to_3edge
                nearest_edge(c) = 1;
                dist_to_edge(c) = dist_to_1edge;
            else
                nearest_edge(c) = 3;
                dist_to_edge(c) = dist_to_3edge;
            end
        else
            dist_to_1edge = abs(perc10(c)-perc10_edge);
            dist_to_2edge = abs(cell_volume(c)-volume_edge);
            if dist_to_1edge<dist_to_3edge
                nearest_edge(c) = 1;
                dist_to_edge(c) = dist_to_1edge;
            else
                nearest_edge(c) = 2;
                dist_to_edge(c) = dist_to_2edge;
            end
        end
        
    end
    
    
    for c = 1:total_cells-1
        for d = c+1:total_cells
            if sqrt(sum((centers(c,:)-centers(d,:)).^2)) < radii(c) + radii(d)
                neighbors{c} = [neighbors{c} d];                
                neighbors{d} = [neighbors{d} c];
            end
        end
    end
        
    improvement = zeros(total_cells,1);

    ol = zeros(length(initial_clusters),1);
    cov = zeros(length(initial_clusters),1);
    
    while 1
        for c = 1:total_cells
            
            for k = 1:length(initial_clusters)
                k_neighbors = neighbors{c}(current_membership(neighbors{c})==k);
                [ol(k) cov(k)] = get_ol_and_cov(hulls{c}, hulls(k_neighbors));
            
            end
            
            ck = current_membership(c);
            if ck == initial_membership(c)
                pk = nearest_edge(c);
            else
                pk = initial_membership(c);                
            end
            
                improvement(c) = real(((cov(pk) - ol(pk)) - (cov(ck) - ol(ck)))/cell_area(c));
            
        end
        
        if any(improvement > improvement_threshold);

            is_good = improvement > improvement_threshold;
            closest_edge = min(dist_to_edge(is_good));
            
            is_best = find(improvement > improvement_threshold & dist_to_edge == closest_edge,1,'first');
            ck = current_membership(is_best);
            if ck == initial_membership(is_best)
                pk = nearest_edge(is_best);
            else
                pk = initial_membership(is_best);                
            end
        
                disp(['moving ' num2str(all_cell_nums(is_best)) ' to cluster ' num2str(pk) ' from ' num2str(ck) ' with improvement factor ' num2str(improvement(is_best))]);
                current_membership(is_best) = pk;
        
        else
            break
        end
        
    end
    
    
    final_clusters = cell(length(initial_clusters),1);
    for k = 1:length(initial_clusters)
        final_clusters{k} = all_cell_nums(current_membership==k);
    end
    
    
end
    
function [current_overlap current_coverage] = get_ol_and_cov(my_hull, neighbor_hulls)    
    
    current_overlap = 0;
    non_overlap_hull = my_hull;
    for k = 1:length(neighbor_hulls)
        overlap_hull = [];
        
        a = polybool('intersection', my_hull(:,1), my_hull(:,2),neighbor_hulls{k}(:,1), neighbor_hulls{k}(:,2));
        if ~isempty(a)
            [overlap_hull(:,1), overlap_hull(:,2)] = polybool('intersection', my_hull(:,1), my_hull(:,2),neighbor_hulls{k}(:,1), neighbor_hulls{k}(:,2));
            current_overlap = current_overlap + poly_area(overlap_hull);

            temp_hull = [];
            temp = polybool('-', my_hull(:,1), my_hull(:,2),neighbor_hulls{k}(:,1), neighbor_hulls{k}(:,2));
            if isempty(temp)
                temp_hull = [];
            else
                [temp_hull(:,1), temp_hull(:,2)] = polybool('-', my_hull(:,1), my_hull(:,2),neighbor_hulls{k}(:,1), neighbor_hulls{k}(:,2));
            end
            non_overlap_hull = temp_hull;
        end

    end
    current_coverage = poly_area(non_overlap_hull);

end