function run_length_approx0(cell_type_char, cell_type_num)

    dend_fns = get_files_with_names_including('./forSrini/', cell_type_char);
    
    
    for n = 1:length(dend_fns)
        
        
        uscore_loc = strfind(dend_fns{n}, '_');
        
        cell_id = dend_fns{n}(1:uscore_loc(1)-1);
        
        char_id_loc = strfind(cell_id, cell_type_char);
        cell_num = [cell_id(1:char_id_loc-1) cell_type_num cell_id(char_id_loc+length(cell_type_num):end)];
        
        surface_fns = get_files_with_names_including('./surface_points/', cell_num);
        
        if ~isempty(surface_fns)
            disp(dend_fns{n});
            
            if length(surface_fns)>1
                error('more than 1 matching surface_points file');
            end
            
            surface_info = whos('-file', ['./surface_points/' surface_fns{1}]);
            
            points_fns = get_files_with_names_including('./point_data/', cell_num);
            
            total_points = 0;
            for k = 1:length(points_fns)
                point_info = whos('-file', ['./point_data/' points_fns{k}]);
                
                total_points = total_points + point_info.size(1);
            end
            
            
            
            est_length = estimate_skeleton_length(total_points, surface_info.size(1), [16.5 16.5 25]);
            
            disp(est_length)
            
            omniSz = 128*[8 16 6];
            
            est_redundancy = 1*prod(omniSz-256) + ...
                2 * 2 * 128 * (prod(omniSz(1:2)-256) + prod(omniSz([1 3])-256) + prod(omniSz([2 3])-256)) + ...
                4 * 4 * 128^2 * sum(omniSz-256) + ...
                8 * 8 * 128^3;
            
            est_redundancy = est_redundancy / ...
                (prod(omniSz-256) + ...
                2 * 128 * (prod(omniSz(1:2)-256) + prod(omniSz([1 3])-256) + prod(omniSz([2 3])-256)) + ...
                4 * 128^2 * sum(omniSz-256) + ...
                8 * 128^3);
            
            est_length = est_length * est_redundancy;
            
            [total_pieces change_type] = analyze_dends(['./forSrini/' dend_fns{n}]);
            
            false_positives = cumsum(change_type(:,2));
            run_length = cumsum(change_type(:,1))/total_pieces*est_length;
            
            plot(false_positives, run_length);
            hold on
            
        end
        
        
    end
        