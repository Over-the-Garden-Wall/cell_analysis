function spreadsheet_BC_info

    C = get_constants;
    
    BCs = unique([C.type.t1, C.type.t2, C.type.t3a, C.type.t3b, C.type.t4]);
    
    
    max_size = 0;
    for b = 1:length(BCs);
        cell_dat = cell_data(BCs(b));
        max_size = max(max_size, length(cell_dat.hull_2d));

    end
    
    
    H = nan(length(BCs)*2,max_size+1);
    for b = 1:length(BCs);
        H(b*2 - (0:1),1) = BCs(b);
        
        cell_dat = cell_data(BCs(b));
        
        h = cell_dat.hull_2d;
        H(b*2-1,2:length(h)+1) = h(:,1) / 1.15;
        H(b*2,2:length(h)+1) = h(:,2) / 1.15;         
        
    end
    
    
    dlmwrite('./spreadsheets/BC_hulls.txt',H);
    
end
        
    
    
    