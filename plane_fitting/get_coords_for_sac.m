function [coords point_id] = get_coords_for_sac(sac_name)


    dir_path = ['./' sac_name '/'];
    
    fns = dir(dir_path);
    
    for n = 1:length(fns)
        if length(fns(n).name) > 5 && strcmp(fns(n).name(1:5),'e2198');
            fn = fns(n).name;
            break
        end
    end
    
    fid = fopen([dir_path fn], 'r');
    F = fread(fid, '*char' )';
    fclose(fid);
    
    relevant_lines = get_strings_starting_with('node id="', F, 80);    

    coords = zeros(length(relevant_lines),3);   
    point_id = zeros(length(relevant_lines),1);   
    
    for n = 1:length(relevant_lines)
        
        relevant_string = relevant_lines{n};
        
        quote_loc = find(relevant_string=='"');
        
%         x_loc = find(relevant_string=='x',1,'first');
%         y_loc = find(relevant_string=='y',1,'first');
%         z_loc = find(relevant_string=='z',1,'first');        
%         V_loc = find(relevant_string=='V',1,'first');

        point_id(n) = str2double(relevant_string(quote_loc(1)+1:quote_loc(2)-1));
        coords(n,1) = str2double(relevant_string(quote_loc(5)+1:quote_loc(6)-1));
        coords(n,2) = str2double(relevant_string(quote_loc(7)+1:quote_loc(8)-1));
        coords(n,3) = str2double(relevant_string(quote_loc(9)+1:quote_loc(10)-1));                
        
    end
    
end