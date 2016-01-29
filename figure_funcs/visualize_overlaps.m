function visualize_overlaps(cell_nos, labels)

    num_cells = length(cell_nos);
    
    max_overlap_level = 3;
    
        
    if ~exist('labels', 'var')
        labels = cell(length(cell_nos),1);
        for n = 1:num_cells
            labels{n} = num2str(cell_nos(n));
        end
    end
    
    if isempty(labels);
        labels = cell(length(cell_nos),1);
        for n = 1:num_cells
            labels{n} = ' ';
        end
    end

    hull_points = cell(num_cells,1);

    num_h_points = 0;
    for k = 1:length(cell_nos)
        c_d = cell_data(cell_nos(k));
%         p = c_d.get_surface;
%         p = p(1:100:end,2:3);
%         hull_points{k} = make_locally_convex_hull(p, rad, true);
        [hull_points{k}(:,1), hull_points{k}(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
%         num_h_points = num_h_points + 
    end
    
%     for k = 1:length(hull_points)
        
    
%     polygons = poly_from_lines(
    
    
    
%     intersection_levels = cell(num_cells,1);
    intersection_levels = cell(max_overlap_level,1);
    intersection_ids = cell(max_overlap_level,1);
    
    intersection_levels{1} = hull_points;
    intersection_ids{1} = (1:length(hull_points))';
    
    for n = 2:max_overlap_level %num_cells
        [intersection_levels{n}, intersection_ids{n}] = hull_intersect(intersection_levels{n-1}, intersection_ids{n-1}, intersection_levels{1});
        if isempty(intersection_levels{n})
%             intersection_levels = intersection_levels(1:n-1);
            break
        end
    end
    
    ilev = max_overlap_level;
    
    c = colormap('jet');
    figure; hold all
    for k = 1:ilev
%     figure; hold all
        for m = 1:length(intersection_levels{k})
            fill(intersection_levels{k}{m}(:,1), intersection_levels{k}{m}(:,2), c(round(k/ilev*size(c,1)),:), 'lineStyle', 'none');
        end
    end
    
    for k = 1:length(intersection_levels{1});
        plot(intersection_levels{1}{k}(:,1), intersection_levels{1}{k}(:,2), 'Color', [0 0 0]);
    end
    
    for k = 1:length(cell_nos)
        text( mean(hull_points{k}(:,1)), mean(hull_points{k}(:,2)), labels{k});
    end
        
end

function [intersection_hulls, intersection_ids] = hull_intersect(hulls, hull_ids, base_hulls)

    intersection_hulls = cell(length(hulls), length(base_hulls));
    is_nonempty_hull = false(size(intersection_hulls));
    
    for k = 1:length(hulls)
        for l = 1:length(base_hulls)
            
            if ~any(l == hull_ids(k,:))
                is_nonempty_hull(k,l) = ~isempty( ...
                    polybool('intersection', hulls{k}(:,1), hulls{k}(:,2), ...
                    base_hulls{l}(:,1), base_hulls{l}(:,2)) );

                if is_nonempty_hull(k,l)
                    [intersection_hulls{k,l}(:,1), intersection_hulls{k,l}(:,2)] ...
                        = polybool('intersection', ...
                        hulls{k}(:,1), hulls{k}(:,2), ...
                        base_hulls{l}(:,1), base_hulls{l}(:,2));
                end
            end
        end
    end
    intersection_hulls = intersection_hulls(is_nonempty_hull(:));
    
    nonempty_inds = find(is_nonempty_hull(:));
    [nonempty_row, nonempty_col] = ind2sub(size(is_nonempty_hull), nonempty_inds);
    intersection_ids = [hull_ids(nonempty_row,:) nonempty_col];
end