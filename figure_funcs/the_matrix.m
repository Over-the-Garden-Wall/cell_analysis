function the_matrix

    C = get_constants;

    types = {'sure_off_sac', 't1','t2','t3a','t3b','t4'};
    type_lbl = {'Off SAC', 'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'};
    
    
    num_types = length(types);
    num_per_type = zeros(num_types,1);
    all_nums = [];
    
    for k = 1:num_types
        num_per_type(k) = length(C.type.(types{k}));
        all_nums = [all_nums C.type.(types{k})];
    end
    
    num_cells = length(all_nums);
    M = zeros(num_cells);
    
    for k = 1:num_cells
        cell_dat = cell_data(all_nums(k));
        for j = 1:num_cells
            if cell_dat.contact_map.isKey(all_nums(j));
                M(k,j) = cell_dat.contact_count( ...
                    cell_dat.contact_map(all_nums(j)));
            end
        end
    end
    
    figure; imagesc(log(M));
    
    hold on
    
    cumnums = cumsum(num_per_type);
    
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    
    
    for k = 1:num_types-1
        plot([1; xlim(2)], [cumnums(k); cumnums(k)], 'LineWidth', 2, 'Color', [.5 1 .5]);
        plot([cumnums(k); cumnums(k)], [1; ylim(2)], 'LineWidth', 2, 'Color', [.5 1 .5]);
    end
        
    set(gca, 'TickDir', 'out');
    
    set(gca, 'FontSize', 20);
    
    set(gca, 'XTick', cumnums - num_per_type/2);
    set(gca, 'XTickLabel', type_lbl);
    
    set(gca, 'YTick', cumnums - num_per_type/2);
    set(gca, 'YTickLabel', type_lbl);
    
        
    set(gcf, 'Position', [0 0 640 640]);
    
%     set(gca, 'LineStyle', ':')

end