function update_all_celldata
    
    C = get_constants;
    
    fns = dir(C.point_dir);
    
    for n = 3:length(fns)
        fn = fns(n).name;
        underscores = find(fn == '_');
        
        cell_num = str2double(fn(underscores(1)+1:underscores(2)-1));
        
        warning('hack in place')
        if cell_num >= 60463        
            c_d = cell_data(cell_num, true);
        end
    end
    
end