function plot_node_diameters(cell_num)

    C = get_constants;
    
    load([C.skele_dir, 's', num2str(cell_num), '.mat']);
    
    fh = figure; scatter(nodes(:,2), nodes(:,3), 1);
    
    
    [x y] = getpts(fh);
    
    x = [x; x(1)];
    y = [y; y(1)];
    
    [x y] = poly2cw(x,y);
    
    close all;
    
    fh = figure; hold all; scatter(nodes(:,2), nodes(:,3), 1);
    plot(x,y, 'lineWidth', 2); title('center point?');
    [soma_x soma_y] = getpts(fh);
    
    is_in_hull = inpolygon(nodes(:,2), nodes(:,3), x, y);
    
    
    close all;
    fh = figure; hold all; scatter(nodes(:,2), nodes(:,3), 1, 'markerEdgeColor', [1 1 1]*.7);   
    
    nodes = nodes(is_in_hull,2:3);
    scatter(nodes(:,1), nodes(:,2), 1, 'markerEdgeColor', [1 0 0]);
    scatter(soma_x, soma_y, 10, 'markerEdgeColor', [0 0 1]);
    plot(x,y, 'lineWidth', 2, 'Color', [0 0 0]);
    
    node_diameter = node_diameter(is_in_hull);
    
    node_dist = sqrt((nodes(:,1)-soma_x).^2 + (nodes(:,2)-soma_y).^2);
    
    [node_dist, sort_ind] = sort(node_dist);
    node_diameter = node_diameter(sort_ind);
    figure; plot(node_dist, node_diameter*50, 'lineWidth', 2);
    xlabel('distance from soma');
    ylabel('node radius');
    
end
    
    
    
    
    