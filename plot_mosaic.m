function plot_mosaic(cell_type)

    C = get_constants;
    cell_nums = C.type.(cell_type);
    num_cells = length(cell_nums);
    
    single_hulls = cell(num_cells,1);
    double_hulls = [];
    triple_hulls = [];
    
    for k = 1:num_cells
        cell_dat = cell_data(cell_nums(k), true);
        single_hulls{k} = cell_dat.hull_2d;
    end
    
    p = 1;
    for k = 1:num_cells-1
        for j = k+1:num_cells
            hull_intersects = get_hull_intersection(single_hulls{k}, single_hulls{j});
            for n = 1:length(hull_intersects)
                double_hulls{p} = hull_intersects{n};
                p = p+1;
            end                        
        end
    end
    num_doubles = p-1;
    
    p = 1;
    for k = 1:num_doubles-1
        for j = k+1:num_doubles                
            hull_intersects = get_hull_intersection(double_hulls{k}, double_hulls{j});
            for n = 1:length(hull_intersects)
                triple_hulls{p} = hull_intersects{n};
                p = p+1;
            end   
        end
    end
    num_triples = p-1;
    
    
            
    figure; hold on
    
    for k = 1:num_cells        
        fill(single_hulls{k}(:,1), single_hulls{k}(:,2), [0 0 1]);
    end
    for k = 1:num_doubles        
        fill(double_hulls{k}(:,1), double_hulls{k}(:,2), [1 1 0], 'EdgeColor', 'none');
    end
    for k = 1:num_triples        
        fill(triple_hulls{k}(:,1), triple_hulls{k}(:,2), [1 0 0], 'EdgeColor', 'none');
    end
    
    
    
    
    prep_figure(gcf,gca);
    set(gca, 'XLim', [.4*10^5 2.4*10^5])
    set(gca, 'YLim', [.4*10^5 2.4*10^5])
    
end


function hull_int = get_hull_intersection(hullA, hullB)
    t = polybool('intersection', hullA(:,1), hullA(:,2), hullB(:,1), hullB(:,2));
    if isempty(t)
        hull_int = [];
        return
    end
    
    [full_hull(:,1) full_hull(:,2)] = polybool('intersection', hullA(:,1), hullA(:,2), hullB(:,1), hullB(:,2));
                    
    nan_locs = [0; find(isnan(t)); length(t)+1];
    num_pieces = length(nan_locs)-1;
    
    hull_int = cell(num_pieces,1);
    for n = 1:num_pieces
        hull_int{n} = full_hull(nan_locs(n)+1:nan_locs(n+1)-1,:);
    end
end