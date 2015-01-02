function [perc perc_error error_rates] = find_best_percentile_split(cell_nums, cell_nums_b)

    if exist('cell_nums_b', 'var');
        cell_nums = {cell_nums, cell_nums_b};
    end

    C = get_constants;
    
    k_dat = cell(2,1);
    for k = 1:2
        k_dat{k} = zeros(100,length(cell_nums));
        
        for n = 1:length(cell_nums{k})
            c = cell_nums{k}(n);
            cell_dat = cell_data(c);
            p = cell_dat.get_surface;
            num_points = size(p,1);
            ipl_depth = sort(C.f(p(:,1)));
            
            k_dat{k}(:,n) = ipl_depth(round((0:99)*num_points/100) + 1);
        end
    end
    
    error_rates = zeros(100,1);
    
    for p = 1:100
        
        ids = [ones(1,size(k_dat{1},2)) zeros(1,size(k_dat{2},2))]; 
        vals = [k_dat{1}(p,:) k_dat{2}(p,:)];
        [vals sort_ord] = sort(vals);
        ids = ids(sort_ord);
        
        num_ones = cumsum(ids);
        num_zeros = cumsum(1-ids);
        
        num_errors = min([num_ones+num_zeros(end)-num_zeros; num_zeros+num_ones(end)-num_ones]);

        error_rates(p) = min(num_errors)/length(ids);
    
    end
    
    [perc_error perc] = min(error_rates);
    perc = perc - 1;
    
end
        