function fig_handles = plot_strat_w_contacts(cell_nums, target_cell, R, color_scheme)

    if ~exist('color_scheme','var')
        color_scheme = [];
    end
    
  
    [fig_handles, colors, ~, y_data, is_valid] = plot_stratification(cell_nums, color_scheme);
    cell_nums = cell_nums(is_valid);
    
    
    add_contacts_to_plot(R, cell_nums, target_cell, y_data, fig_handles(2), fig_handles(3), colors);
    

end