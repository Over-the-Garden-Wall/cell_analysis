function [nodes edges] = get_skeleton(cell_name)

    
    res = [16.5 16.5 25];

    
    
    

    dir_path = ['./' cell_name '/'];
    
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

    num_points = length(relevant_lines);
    
    
    nodes = zeros(num_points,3);   
    point_id = zeros(num_points,1);   
    
    for n = 1:num_points
        
        relevant_string = relevant_lines{n};
        
        quote_loc = find(relevant_string=='"');
        
%         x_loc = find(relevant_string=='x',1,'first');
%         y_loc = find(relevant_string=='y',1,'first');
%         z_loc = find(relevant_string=='z',1,'first');        
%         V_loc = find(relevant_string=='V',1,'first');

        point_id(n) = str2double(relevant_string(quote_loc(1)+1:quote_loc(2)-1));
        nodes(n,1) = str2double(relevant_string(quote_loc(5)+1:quote_loc(6)-1));
        nodes(n,2) = str2double(relevant_string(quote_loc(7)+1:quote_loc(8)-1));
        nodes(n,3) = str2double(relevant_string(quote_loc(9)+1:quote_loc(10)-1));                
        
    end
    
    for k = 1:3
        nodes(:,k) = nodes(:,k)*res(k);
    end

    
    id2index = sparse(point_id, ones(num_points,1), 1:num_points, max(point_id),1);
    
    
    
    
    edge_str = get_strings_starting_with('edge ', F, 80);
    
%     full_len = 0;
%     figure; xyaxis = gca; hold on
%     figure; yzaxis = gca; hold on
    
    edges = zeros(length(edge_str),2);
    for n = 1:length(edge_str)
        relevant_string = edge_str{n};
        quote_loc = find(relevant_string=='"');

        src_edge = str2double(relevant_string(quote_loc(1)+1:quote_loc(2)-1));
        dest_edge = str2double(relevant_string(quote_loc(3)+1:quote_loc(4)-1));
       
        edges(n,:) = [id2index(src_edge) id2index(dest_edge)];
        
%         full_len = full_len + sqrt(sum((nodes(id2index(src_edge),:) - nodes(id2index(dest_edge),:)).^2));
%         
%         plot(xyaxis, nodes(id2index([src_edge dest_edge]),1), nodes(id2index([src_edge dest_edge]),2));
%         plot(yzaxis, nodes(id2index([src_edge dest_edge]),2), nodes(id2index([src_edge dest_edge]),3));
        
        
    end
    
    
    