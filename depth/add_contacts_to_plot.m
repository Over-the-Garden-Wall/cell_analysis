function add_contacts_to_plot(conns, cell_nums, target_cell, y_data, all_plot_handle, sep_plot_handle, colors)

    C = get_constants;
    
    
    num_cells = length(cell_nums);
    
    y = cell(num_cells,1);
    depth = cell(num_cells,1);
    
    for n = 1:num_cells
        if target_cell == 1
            cell_dat = cell_data(cell_nums(n));
            depth_loc = find(cumsum(cell_dat.stratification)>=.1, 1,'first');
            depth{n} = depth_loc - 1 + C.strat_x(1);
            
        else

            sub_conn = pick_conns(conns,target_cell, cell_nums(n));
            depth{n} = ceil(C.f(double(sub_conn{1}(4,:))));

            depth{n} = depth{n}(depth{n} >= C.strat_x(1) & depth{n} <= C.strat_x(end));
    %         depth = depth - C.strat_x(1) + 1;
        end
        y{n} = y_data(depth{n} - C.strat_x(1) + 1,n);
    end
        
    figure(all_plot_handle)
    hold on
    for n = 1:num_cells    
        plot(depth{n}, y{n}, '*', 'Color', colors(n,:));
    end
    hold off
        
    figure(sep_plot_handle)

    p = ceil(sqrt(num_cells));
    q = ceil(num_cells/p);
    for n = 1:num_cells
        subplot(p,q,n); hold on
        plot(depth{n}, y{n}, '*', 'Color', colors(n,:));
        hold off;
    end
    
end