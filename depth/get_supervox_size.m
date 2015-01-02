function get_supervox_size(cell_name, block_num)

    root_dir = ['./forSrini/' cell_name '/'];

    fns = get_files_with_names_including(root_dir, [cell_name '_' num2str(block_num)]);
    
    num_fns = length(fns);
    
    sprvox_num = zeros(num_fns,1);
    sprvox_size = zeros(num_fns,1);
    
    
    for n = 1:num_fns
        
        uscore_loc = strfind(fns{n}, '_');
        
        if length(uscore_loc) > 3
            sprvox_num(n) = str2double(fns{n}(uscore_loc(3)+1:uscore_loc(4)-1));
        

            surface_info = whos('-file', [root_dir fns{n}]);
            sprvox_size(n) = surface_info.size(1);
        end
    end
    
    [sprvox_num, sort_ord] = sort(sprvox_num);
    sprvox_size = sprvox_size(sort_ord);
    
    for n = num_fns-1:-1:1        
        if sprvox_num(n) == sprvox_num(n+1)
            sprvox_size(n) = sprvox_size(n) + sprvox_size(n+1);
            sprvox_num(n+1) = [];
            sprvox_size(n+1) = [];
        end
    end
    
    
    save([root_dir 'sizes.mat'],'sprvox_num','sprvox_size');
end
            
            