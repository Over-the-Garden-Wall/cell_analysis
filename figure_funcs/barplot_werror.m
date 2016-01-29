function h = barplot_werror(bar_data, error_data)

    bh = bar(bar_data); hold on;

    eh = zeros(1,length(bh));
    for n = 1:length(bh);
        
        bar_xverts = get(get(bh(n), 'Children'), 'XData');
        bar_width = bar_xverts(3,1) - bar_xverts(1,1);
        error_x = mean(bar_xverts([1 3],:));
        
        eh(n) = errorbar(error_x, bar_data(:,n)', error_data(:,n)', 'lineWidth', 2, 'color', [0 0 0], 'lineStyle', 'none');
        set(eh(n), 'LData', bar_data(:,n));
        
        error_children = get(eh(n), 'Children');
        error_lines = get(error_children(2), 'XData');
        if error_lines(5) - error_lines(4) > bar_width*.6;
            for t = 1:9:length(error_lines)
                error_lines(t + [3 6]) = error_lines(t) - bar_width*.3;
                error_lines(t + [4 7]) = error_lines(t) + bar_width*.3;
            end
            set(error_children(2), 'XData', error_lines);
        end
        
        
%         t = get(eh(n), 'LData');
        
        
%         'LData', zeros(1,size(error_data,1))
        
    end
    
    delete(bh);
    bh = bar(bar_data, 'BarWidth', .8);
    
    h = [bh; eh];

end