function conns = connections_within_region(cell_num, guide_region)

    fh = figure; plot_cells(cell_num, 1, .01, [0 0 0]);
    
    if exist('guide_region', 'var') && ~isempty(guide_region)
        plot(guide_region(:,1), guide_region(:,2), 'r', 'lineWidth', 2);
    end
    
    [x, y] = getpts(fh);
    
    [x, y] = poly2cw(x,y);
    
    c_d = cell_data(cell_num);
    all_conns = double(c_d.contacts);
    
    is_in = inpolygon(all_conns(4,:)', all_conns(5,:), x, y);
    
    conns = all_conns(:,is_in);
    
    scatter(conns(4,:), conns(5,:));
end