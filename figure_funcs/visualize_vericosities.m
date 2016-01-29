function visualize_vericosities(cell_num, vericosity_threshold, pixel_size, ground_truth)

    C = get_constants;

    skele = load([C.skele_dir, 's', num2str(cell_num), '.mat']);
    
    
    
    resamp_nodes = skele.nodes(:,2:3);
    
    n_min = min(resamp_nodes)-1;
    
    for d = 1:2
        resamp_nodes(:,d) = ceil((resamp_nodes(:,d) - n_min(d))/pixel_size(d));
    end;
    
    [resamp_nodes, sort_ind] = sortrows(resamp_nodes);
    skele.node_diameter = skele.node_diameter(sort_ind);
    
    n_max = max(resamp_nodes);
    
    im = zeros(n_max);
    im_count = zeros(n_max);
    
    for n = 1:size(resamp_nodes,1)
        im(resamp_nodes(n,1), resamp_nodes(n,2)) = ...
            im(resamp_nodes(n,1), resamp_nodes(n,2)) + skele.node_diameter(n);        
        
        im_count(resamp_nodes(n,1), resamp_nodes(n,2)) = ...
            im_count(resamp_nodes(n,1), resamp_nodes(n,2)) + 1;
    end
    
    im = im./im_count;
    im(isnan(im)) = 0;
    
    figure; h = imagesc(im, [0 vericosity_threshold]);
    
    c = zeros(100,3);
    c(1:100,1) = (1:100)'/100;
    c(end,:) = [1 1 1];
    
    colormap(c); %colorbar;
    
    hold on;
    ground_truth = ground_truth(:,2:3);
    for d = 1:2
        ground_truth(:,d) = ceil((ground_truth(:,d) - n_min(d))/pixel_size(d));
    end
    
    scatter(ground_truth(:,2), ground_truth(:,1), 'markerEdgeColor', [.5 .5 1]);
    
    
end
    