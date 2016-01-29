function visualize_cell_locations(cell_nums)

    figure; hold on;
    
%     C = get_constants;
    colmp = make_colormap(length(cell_nums), 5);
    
    for n = 1:length(cell_nums)
        c_d = cell_data(cell_nums(n));
        midp = c_d.get_midpoint;
        
        plot(c_d.hull_2d(:,1), c_d.hull_2d(:,2), 'LineWidth', 2, 'Color', colmp(n,:));
        
        
        text(midp(2), midp(3), num2str(cell_nums(n)), 'Color', colmp(n,:)*.7);
    end
    
end