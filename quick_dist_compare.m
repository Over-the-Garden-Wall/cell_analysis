function M = quick_dist_compare(nodes, threshold)

    

    num_nodes = size(nodes,1);
    node_d = size(nodes,2);

    
    sort_order = cell(1,node_d);
    
    candidate_range = zeros(num_nodes, 2, node_d);
    
    for d = 1:node_d                
        gridloc = [floor((nodes(:,d) - min(nodes(:,d)))/threshold) + 2; Inf];
        
        [sorted_grid, sort_order{d}] = sort(gridloc);
        
        lbound = 1;
        ubound = 1;
        for curr_node = 1:num_nodes
            
            while sorted_grid(curr_node) - sorted_grid(lbound) > 1
                lbound = lbound+1;
            end
            while sorted_grid(ubound) - sorted_grid(curr_node) < 1
                ubound = ubound+1;
            end
            candidate_range(curr_node, :, d) = [lbound ubound-1];
            
        end
    end
    
    M = false(num_nodes);
    for n = 1:num_nodes
        candidate_list = sort_order{1}(candidate_range(n, 1, 1):candidate_range(n, 2, 1));
        candidate_list(candidate_list <= n) = [];
        for d = 2:node_d
            candidate_list = list_intersection(candidate_list, ...
                sort_order{d}(candidate_range(n, 1, d):candidate_range(n, 2, d)));
            candidate_list(candidate_list <= n) = [];
        end
        for c = candidate_list
            dist = sqrt(sum((nodes(c,:) - nodes(n,:)).^2));
            if dist <= threshold
                M(n,c) = true;
                M(c,n) = true;
            end
        end
    end
        
end

function out_list = list_intersection(list_a, list_b)
     list_a = sort(list_a);
     list_b = sort(list_b);
     
     a_counter = 1;
     b_counter = 1;
     is_intersect = false(size(list_a));
     while a_counter <= length(list_a) && b_counter <= length(list_b)
         if list_a(a_counter) == list_b(b_counter)
             is_intersect(a_counter) = true;
             a_counter = a_counter+1;
             b_counter = b_counter+1;
         elseif list_a(a_counter) > list_b(b_counter)
             b_counter = b_counter+1;
         elseif list_a(a_counter) < list_b(b_counter)
             a_counter = a_counter+1;
         end
     end
     out_list = list_a(is_intersect);
end
         
             
     
     

end

        
    
    