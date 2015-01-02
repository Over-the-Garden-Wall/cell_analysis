function rel_perc = percent_strat_within_range(cell_nums, strat_range, relevant_portion)

    C = get_constants;
    
    if ~exist('relevant_portion','var') || isempty(relevant_portion)
        relevant_portion = [-Inf Inf];
    end
    
    num_cells = length(cell_nums);
    rel_perc = zeros(num_cells,1);
    
    for c = 1:num_cells
        
        c_d = cell_data(cell_nums(c));
        
        p = c_d.get_surface;
        d = C.f(p(:,1));
        
        d = d(d>relevant_portion(1) & d<relevant_portion(2));
        
        for k = 1:size(strat_range,1)
            rel_perc(c) = rel_perc(c) + sum(d>strat_range(k,1) & d<strat_range(k,2));
        end
        rel_perc(c) = rel_perc(c);%/length(d);
        
    end
    
end
            
            