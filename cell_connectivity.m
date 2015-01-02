function density = cell_connectivity(ref_cells, off_cells)

    
    C = get_constants;
    
    num_refs = length(C.type.off_sac);
    
    d = cell(length(ref_cells),1);
    
    size_max = [0 0 0];
    for n = 1:num_refs
        d{n} = get_contacts_all(ref_cells(n), off_cells);
        size_max = max([size_max; size(d{n})]);
    end
    
    density = zeros(size_max);
    for n = 1:num_refs
        dSz = size(d{n});
        density(1:dSz(1), 1:dSz(2), 1:dSz(3)) = density(1:dSz(1), 1:dSz(2), 1:dSz(3)) + ...
            d{n};
    end
    
    
end
    