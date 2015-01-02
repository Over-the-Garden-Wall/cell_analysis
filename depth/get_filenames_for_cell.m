function [seg_list_fns, omni_file_fns] = get_filenames_for_cell(cell_id)

    seg_root = '/omniData/e2198/completedCells/';
    om_root = '/omniData/e2198/';
    om_fn_root = 'e2198_';
    

    dir_files = dir([seg_root cell_id]);
    
    is_valid = false(length(dir_files),1);
    
    for n = 1:length(dir_files)
        k = strfind(dir_files(n).name, '_');
        if ~isempty(k) && k(1) < 5
            is_valid(n) = true;
        end
        
    end
    
    dir_files = dir_files(is_valid);
    
    num_files = length(dir_files);
    
    seg_list_fns = cell(num_files,1);
    omni_file_id = cell(num_files,1);    
    omni_file_fns = cell(num_files,1);
    
    for n = 1:num_files
        seg_list_fns{n} = [seg_root cell_id '/' dir_files(n).name];
        
        under_loc = strfind(dir_files(n).name, '_');
        omni_file_id{n} = dir_files(n).name(1:under_loc-1);
    end
        
    
    
    
    dir_files = dir(om_root);
    
    
    is_valid = false(length(dir_files),1);
    
    for n = 1:length(dir_files)
        
        if ~isempty(strfind(dir_files(n).name, 'files'))
            is_valid(n) = true;
        end
        
    end
    
    dir_files = dir_files(is_valid);
    
%     try
    
    for n = 1:num_files
        
        my_flag = true;
        k = 0;
        while my_flag
            k = k+1;
            uscore_loc = strfind(dir_files(k).name,'_');            
            if strcmp(dir_files(k).name(uscore_loc(1)+1:uscore_loc(2)-1),omni_file_id{n})
                my_flag = false;
            end
                        
        end
        
        omni_file_fns{n} = [om_root dir_files(k).name];
    end
       
%     catch
%         
%         disp(omni_file_id{n})
%         disp(dir_files(k).name)
%     end
    
            
end
    