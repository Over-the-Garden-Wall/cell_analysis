function bc_gc_connectivity_scatter(bcs, gc1, gc2)

    num_bcs = length(bcs);

    bc_cell = cell(num_bcs, 1);
    for k = 1:num_bcs
        bc_cell{k} = bcs(k);
    end


    M1 = points_within_hull_intersect(bc_cell, {gc1});
    M2 = points_within_hull_intersect(bc_cell, {gc2});
    
    T1 = contact_total(bc_cell, {gc1});
    T2 = contact_total(bc_cell, {gc2});
    
    scatter(T1./M1, T2./M2);
    
end