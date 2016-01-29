function conns = contacts_near_location(cell_num, distance_threshold, guide_region)

    fh = figure; plot_cells(cell_num, 1, .01, [0 0 0]);
    ph = gca;
    
    if exist('guide_region', 'var') && ~isempty(guide_region)
        plot(guide_region(:,1), guide_region(:,2), 'r', 'lineWidth', 2);
    end
    
    x = 1;
    all_coords = zeros(0,3);
    curr_xlim = get(ph, 'XLim');
    curr_ylim = get(ph, 'YLim');
    
    
    c_d = cell_data(cell_num);
    p = c_d.get_surface;
    
    while 1;
        title('Define Zoom Box');
        [x, y] = getpts(fh);
        if length(x) < 2
            break
        end
        x = x(end-1:end);
        y = y(end-1:end);
        
        set(ph, 'XLim', [min(x) max(x)]);
        set(ph, 'YLim', [min(y) max(y)]);
        
        title('Select Points!');
        [x, y] = getpts(fh);
        if ~isempty(x)
            title('finding points...');
            coords = zeros(length(x),3);
            
            for n = 1:length(x)
                d = sqrt((p(:,2) - x(n)).^2 + (p(:,3) - y(n)).^2);
                [dummy, min_ind] = min(d);
                coords(n,:) = p(min_ind,:);
            end

            scatter(coords(:,2),coords(:,3), '*', 'markerEdgeColor', [.8 .8 0]);
            all_coords = [all_coords; coords];   
        
        end
        
        set(ph, 'XLim', curr_xlim);
        set(ph, 'YLim', curr_ylim);
        
    end
    
    all_conns = double(c_d.contacts);
    
    is_near = false(1, size(all_conns,2));
    for n = 1:size(coords,1);
        d = sqrt((all_conns(3,:) - coords(n,1)).^2 ...
            + (all_conns(4,:) - coords(n,2)).^2 ...
            + (all_conns(5,:) - coords(n,3)).^2);
        is_near = is_near | d < distance_threshold;
    end
    
    conns = all_conns(:,is_near);
    
    scatter(conns(4,:), conns(5,:));
    
end