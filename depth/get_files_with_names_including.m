function matching_fns = get_files_with_names_including(directory, find_str)
    f_info = dir(directory);
    
    fns = cell(length(f_info),1);
    is_valid = false(length(f_info),1);
    for n = 1:length(f_info)
        fns{n} = f_info(n).name;
        if ~isempty(strfind(fns{n},find_str))
            is_valid(n) = true;
        end
    end
    
    matching_fns = fns(is_valid);
end
            