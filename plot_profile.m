function plot_mean_strat(cell_nos, invert_flag)

    if ~exist(invert_flag,'var');
        invert_flag = false;
    end

    C =get_constants;

    num_types = length(cell_nos);
    
    strats = zeros(length(C.strat_x),num_types);
    
    for k = 1:num_types
    
        for cell_no = cell_nos{k}

            c_d = cell_data(cell_no);

            C = get_constants;

            s = c_d.stratification;
            strats(1:length(s),k) = strats(1:length(s),k) + s/length(cell_nos{k});
            
        end
        
    end
    figure; hold all;
    for k = 1:num_types
        if ~invert_flag
            plot(C.strat_x, strats(:,k), 'lineWidth', 2);
        else
            plot(strats(:,k), C.strat_x, 'lineWidth', 2);
        end
    end
    
end
end