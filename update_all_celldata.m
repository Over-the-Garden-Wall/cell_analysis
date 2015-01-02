function update_all_celldata
    
    C = get_constants;
    types = fieldnames(C.type);
    
    cell_nums = [];
    for k = 1:length(types)
            cell_nums = [cell_nums C.type.(types{k})];
    end
    for c = unique(cell_nums);
        cell_dat = cell_data(c,true);
    end
end