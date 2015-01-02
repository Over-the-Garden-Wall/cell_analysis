function lines = get_strings_starting_with(start_str, full_str, get_length)

    start_locs = strfind(full_str,start_str);
    
    lines = cell(length(start_locs),1);
    
    for n = 1:length(start_locs)
        
        lines{n} = full_str(start_locs(n):start_locs(n)-1+get_length);
    end
    
end

    