function plot_individual_strats(cell_nums, xlim, label_cells)
    
    if ~exist('label_cells','var') || isempty(label_cells)
        label_cells = false;
    end

    C = get_constants;
    num_cells = length(cell_nums);
    
    x_plots = ceil(sqrt(num_cells/3));
    y_plots = ceil(num_cells/x_plots);
    
    figure
    for k = 1:num_cells
        c_d = cell_data(cell_nums(k));
        
        s = c_d.stratification;
        if length(s) > length(C.strat_x)
            s = s(1:length(C.strat_x));
        end
        
        subplot(y_plots, x_plots, k);
        plot(C.strat_x(1:length(s)), s, 'LineWidth', 2);
        set(gca, 'XLim', xlim);
        
        if label_cells
            title([num2str(cell_nums(k)) ' - ' cell_type(cell_nums(k))]);
        else
            title(num2str(cell_nums(k)));
        end
    end
    
end
    