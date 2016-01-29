function skele = reduce_skeleton(skele, reduction_factor)

    num_nodes = size(skele.nodes,1);
    num_edges = size(skele.edges,1);
    node_instances = zeros(num_nodes,1);
    
    for n = 1:num_edges
        node_instances(skele.edges(n,1)) = node_instances(skele.edges(n,1)) + 1;
        node_instances(skele.edges(n,2)) = node_instances(skele.edges(n,2)) + 1;
    end
    
    is_noncritical = find(node_instances == 2);
    num_to_remove = min(round((1-reduction_factor)*num_nodes), length(is_noncritical));
    
    to_remove = is_noncritical(round((1:num_to_remove)/num_to_remove * length(is_noncritical)));
    
    node_rename = 1:num_nodes;
    for n = to_remove'
        edges_to_merge = [find(skele.edges(:,1)==n); find(skele.edges(:,2)==n)];
        new_edge = skele.edges(edges_to_merge,:);
        new_edge = new_edge(new_edge(:)~=n)';
        skele.edges(edges_to_merge(1),:) = new_edge;
        skele.edges(edges_to_merge(2),:) = [-99 -99];
        
        node_rename(n) = 0;
        node_rename(n+1:end) = node_rename(n+1:end) - 1;
    end
    
    skele.nodes(to_remove,:) = [];
    skele.node_diameter(to_remove,:) = [];
    
    skele.edges(skele.edges(:,1) == -99, :) = [];        
    skele.edges = node_rename(skele.edges);
    
    
    
end
    