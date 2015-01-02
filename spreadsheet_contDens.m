function spreadsheet_contDens

    C = get_constants;
    
    BCs = unique([C.type.t1, C.type.t2, C.type.t3a, C.type.t3b, C.type.t4]);
    SACs = C.type.sure_off_sac;
    
    C = zeros(length(BCs)+1,length(SACs)+1);
    H = zeros(length(BCs)+1,length(SACs)+1);
    D = zeros(length(BCs)+1,length(SACs)+1);
    
    C(2:end,1) = BCs;
    H(2:end,1) = BCs;
    D(2:end,1) = BCs;
    
    C(1,2:end) = SACs;
    H(1,2:end) = SACs;
    D(1,2:end) = SACs;
    
    
    for s = 1:length(SACs)
        
        sc = SACs(s);
        cell_dat = cell_data(sc);
        sm = cell_dat.get_midpoint(true) / 1.15; %shrinkage
        [total_contact, total_vox_in_hull] = get_contact_density_whulls(sc, BCs);
        
        for b = 1:length(BCs);
            bc = BCs(b);
            cell_dat = cell_data(bc);
            bm = cell_dat.get_midpoint(false) / 1.15; %shrinkage
            
            C(b+1,s+1) = total_contact(b);
            H(b+1,s+1) = total_vox_in_hull(b);
            D(b+1,s+1) = sqrt((bm(2)-sm(2))^2 + (bm(3)-sm(3))^2);
        end
    end
    
    dlmwrite('./spreadsheets/contact_totals.txt',C);
    dlmwrite('./spreadsheets/hull_intersect_totals.txt',H);
    dlmwrite('./spreadsheets/distance_to_cells.txt',D);
    
end
        
    
    
    