function trimmed_skele = remove_trivial_skeleton_nodes(skeleton, critical_node_list, soma_node)


    if ~any(soma_node==critical_node_list)
        warning('adding soma_node to critical_node_list')
        critical_node_list = [soma_node critical_node_list];
    end

    G = skeleton2graph(skeleton);
    
    trimmed_skele.edges = [];
    
    trimmed_skele.nodes = zeros(length(critical_node_list),3);
    
    distances = graphshortestpath(G,soma_node);
    is_valid = ~isinf(distances);
    
    while any(is_valid)
        valid_list = find(is_valid);
        [~, most_distal_ind] = max(distances(is_valid));
        most_distal_node = valid_list(most_distal_ind);
        [~, gpath] = graphshortestpath(G,soma_node, most_distal_node);
        
        last_node = soma_node;
        
        for n = 2:length(gpath);
            if any(gpath(n) == critical_node_list)
                trimmed_skele.edges = [trimmed_skele.edges; last_node gpath(n)];                
                last_node = gpath(n);
            end
        end
        is_valid(gpath) = false;
    end
        
    trimmed_skele.edges = unique(trimmed_skele.edges,'rows');
        
    for n = 1:length(critical_node_list);                
        trimmed_skele.nodes(n,:) = skeleton.nodes(critical_node_list(n),:);                
    end
    
    inv_critical_node_list = sparse(critical_node_list, ones(size(critical_node_list)), 1:length(critical_node_list));
    trimmed_skele.edges = inv_critical_node_list(trimmed_skele.edges);
end
    
    
        

    