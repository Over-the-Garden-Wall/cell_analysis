function display_node_cloud(cell_num)

    
    
    C = get_constants;
    
    load([C.skele_dir 's' num2str(cell_num) '.mat']);
    nodes(:,2) = nodes(:,2) - mean(nodes(:,2));
    nodes(:,3) = nodes(:,3) - mean(nodes(:,3));
    
    num_nodes = size(nodes,1);
    is_branch = false(num_nodes,1);
    is_leaf = false(num_nodes,1);
    
    
    for k = 1:size(nodes,1)
        is_branch(k) = sum(edges(:)==k) >= 3;
        is_leaf(k) = sum(edges(:)==k) == 1;
    end
    
    minvals = min(nodes(:,2:3));
    maxvals = max(nodes(:,2:3));
        nodes(:,1) = C.f(nodes(:,1));
        
        zvals = [min(nodes(:,1)), 100];
        
        diff = max(maxvals - minvals);
        maxvals = minvals + diff;
        
    
    figure; set(gcf, 'Position', [0 0 1100 300])
    
    subplot(1,3,1); hold on
    title([num2str(cell_num) ' x-y']);
    scatter(nodes(is_branch,2), nodes(is_branch,3), '*', 'markerEdgeColor', [1 0 0]);
    scatter(nodes(is_leaf,2), nodes(is_leaf,3), '*', 'markerEdgeColor', [0 0 1]);
    set(gca, 'XLim', [minvals(1) maxvals(1)], 'YLim', [minvals(2) maxvals(2)]);
    
    subplot(1,3,2); hold on
    title('z-y');
    scatter(nodes(is_branch,1), nodes(is_branch,3), '*', 'markerEdgeColor', [1 0 0]);
    scatter(nodes(is_leaf,1), nodes(is_leaf,3), '*', 'markerEdgeColor', [0 0 1]);
    set(gca, 'XLim', zvals, 'YLim', [minvals(2) maxvals(2)]);
    
    subplot(1,3,3); hold on
    title('z-x');    
    scatter(nodes(is_branch,1), nodes(is_branch,2), '*', 'markerEdgeColor', [1 0 0]);
    scatter(nodes(is_leaf,1), nodes(is_leaf,2), '*', 'markerEdgeColor', [0 0 1]);
    set(gca, 'XLim', zvals, 'YLim', [minvals(1) maxvals(1)]);

end
    