function spreadsheet_BC_info

    C = get_constants;
    
    BCs = [C.type.t1, C.type.t2, C.type.t3a, C.type.t3b, C.type.t4];
    BC_type = [ones(size(C.type.t1)), 2*ones(size(C.type.t2)), 31*ones(size(C.type.t3a)), 32*ones(size(C.type.t3b)), 4*ones(size(C.type.t4))];
    
    M = zeros(length(BCs), 6);
    S = zeros(length(BCs),101);
    
    
        for b = 1:length(BCs);
            bc = BCs(b);
            cell_dat = cell_data(bc);
            p = cell_dat.get_surface;
            d = C.f(p(:,1));
            d(d<0) = [];
            d = sort(d);
            
            M(b,1) = bc;
            S(b,1) = bc;
            M(b,2) = BC_type(b);
            M(b,3:5) = d(ceil([.1 .25 .75]*length(d)));
            M(b,6) = cell_dat.V;
            s = cell_dat.stratification(21:end);
            S(b,2:length(s)+1) = s*cell_dat.SA;
        end
    
    dlmwrite('./spreadsheets/BC_info.txt',M);
    dlmwrite('./spreadsheets/BC_strats.txt',S);
    
end
        
    
    
    