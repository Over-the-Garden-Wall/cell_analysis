function areas = get_contact_areas(main_cell, other_cells)

    load('/net/omicfs/home/matthew/stratification/conns.mat')
    
    to_reorder = conns(1,:) > conns(2,:);
    conns(1:2,to_reorder) = conns([2 1], to_reorder);
    
    areas = zeros(length(other_cells),1);
    
    for c = 1:length(other_cells)
        
        if main_cell < other_cells(c)
            cells = [main_cell, other_cells(c)];
        else
            cells = [other_cells(c), main_cell];
        end
        
        is_valid = conns(1,:) == cells(1) & conns(2,:) == cells(2);
        
        areas(c) = sum(conns(3,is_valid));
        
    end
end