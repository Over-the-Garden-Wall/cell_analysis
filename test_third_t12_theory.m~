function C = test_third_t12_theory

    C = get_constants;

    cell_nums = [C.type.t1 C.type.t2];
    num_cells = length(cell_nums);
    
    perc_data = zeros(num_cells,1);
    for k = 1:length(cell_nums)
        
        cell_dat = cell_data(cell_nums(k));
        cumstrat = cumsum(cell_dat.stratification);
        
        cumstrat = cumstrat/cumstrat(end);
        
        perc_data(k,1) = find(cumstrat>.2,1,'first');
        perc_data(k,2) = find(cumstrat>.8,1,'first');
        
    end
    
    diff = perc_data(:,2)-perc_data(:,1);
    figure; hist(diff, 0:25);
    
   
        
    C.type.C1 = cell_nums(diff < 13)
    C.type.C2 = cell_nums(diff >= 13 & diff < 16)
    C.type.C3 = cell_nums(diff >= 16)

    














end