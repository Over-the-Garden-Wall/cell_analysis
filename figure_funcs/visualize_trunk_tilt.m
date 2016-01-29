function visualize_trunk_tilt(cell_nums)

    figure; a = gca; hold on;
    
    C = get_constants;
    colmp = make_colormap(length(cell_nums), 5)*.9;
    
    for n = 1:length(cell_nums)
        c_d = cell_data(cell_nums(n));
        
        p = c_d.get_surface;
        d = C.f(p(:,1));
        minval = min(d);
        
        is_trunky = d<minval + 5;
        def_not_trunky = d > median(d);
        
        tp = mean(p(is_trunky,:));
        midp = mean(p(def_not_trunky,:));
        
        
        plot(a, [midp(2); tp(2)], [midp(3); tp(3)], 'LineWidth', 2, 'Color', colmp(n,:));
        scatter(a, midp(2), midp(3), 'markerEdgeColor', colmp(n,:));
%         plot(a, c_d.hull_2d(:,1), c_d.hull_2d(:,2), 'LineWidth', 2, 'Color', colmp(n,:));        
        
        text(a, midp(2), midp(3), num2str(cell_nums(n)), 'Color', colmp(n,:)*.7);
        
        figure; plot_cells(cell_nums(n), 2, .01, [0 0 0]);
        hold on;
        scatter(midp(1), midp(3), '*', 'markerEdgeColor', [1 0 0]);
        scatter(tp(1), tp(3), '*', 'markerEdgeColor', [1 0 1]);
        figure; plot_cells(cell_nums(n), 3, .01, [0 0 0]);
        hold on;
        scatter(midp(1), midp(2), '*', 'markerEdgeColor', [1 0 0]);
        scatter(tp(1), tp(2), '*', 'markerEdgeColor', [1 0 1]);
        
    end
    
end