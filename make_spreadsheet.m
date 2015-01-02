C = get_constants;

types = {'t1','t2','t3a','t3b','t4'};
type_id = [1 2 31 32 4];
sure_types = {'sure_t3a','sure_t4'};

percentiles = 5:5:100;
percentiles(end) = percentiles(end) - .0001;


M = [];

for k = 1:length(types);
    
    cell_nums = C.type.(types{k});
    
    for n = 1:length(cell_nums);
        c = cell_nums(n);
        
        M(end+1,1) = c;
        M(end,2) = type_id(k);
        
        is_sure = 0;
        for t = 1:length(sure_types);
            if any(C.type.(sure_types{t})==c)
                is_sure = 1;
            end
        end
        
        M(end,3) = is_sure;
        
        c_d = cell_data(c);
        
        M(end,4) = c_d.SA;
        M(end,5) = c_d.V;
        M(end,6) = poly_area(c_d.hull_2d);
        
        
        
        
        
        cumstrat = cumsum(c_d.stratification);
        
        for t = 1:length(percentiles)
            f = find(cumstrat >= percentiles(t)/100,1,'first');
            M(end,6+t) = f+C.strat_x(1)-1;
        end
        
    end
        
end

dlmwrite('./cell_data.txt',M);