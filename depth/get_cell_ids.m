function ids = get_cell_ids(identifying_flag)

    root_dir = '/omniData/e2198/completedCells/';
    
    file_dirs = dir(root_dir);
    
    is_valid = false(length(file_dirs),1);
    
    for n = 1:length(file_dirs)
        if ~isempty(strfind(file_dirs(n).name,identifying_flag))
            is_valid(n) = true;
        end
    end
    
    file_dirs = file_dirs(is_valid);
    
    ids = cell(length(file_dirs),1);
    for n = 1:length(file_dirs)
        ids{n} = file_dirs(n).name;
    end
end
    