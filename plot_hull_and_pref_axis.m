function plot_hull_and_pref_axis(cell_num)
    c = 150000;
    cell_dat = cell_data(cell_num);
    
    figure; hold all
    
    plot(cell_dat.hull_2d(:,1), cell_dat.hull_2d(:,2), 'LineWidth', 2);
    
    if ~isempty(cell_dat.dist_axis);
        
        cell_mid = cell_dat.get_midpoint(true);
        da = cell_dat.dist_axis;
        cm = cell_mid(2:3);
        
        plot([cm(1); cm(1)+c*da(1)], [cm(2); cm(2)+c*da(2)], 'LineWidth', 2);
        
    end
    
    ax = gca;
    xlims = get(ax, 'XLim');
    ylims = get(ax, 'YLim');
    
    xmid = sum(xlims)/2;
    ymid = sum(ylims)/2;
    
    max_dist = max(xlims(2)-xlims(1), ylims(2)-ylims(1));
    
    set(ax, 'XLim', [xmid-max_dist/2, xmid+max_dist/2]);
    set(ax, 'YLim', [ymid-max_dist/2, ymid+max_dist/2]);
    
    prep_figure(gcf,gca);