function M = contact_matrix(cell_nums)

    num_cells = length(cell_nums);
    
    M = zeros(num_cells);
    
    for n = 1:num_cells
        c_d = cell_data(cell_nums(n));
        conts = c_d.contacts;
        M(n,n) = polyarea(c_d.hull_2d(:,1), c_d.hull_2d(:,2))/2/10^6;
        
        for m = n+1:num_cells
            is_me = conts(1,:) == cell_nums(m);
            if any(is_me)
                M(n,m) = sum(conts(2,is_me));
            end
        end
    end
    M = M+M';
end
            
            
        
        