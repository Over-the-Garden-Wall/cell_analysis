function tree = skele2tree(skele, root_node)

    if ~exist('root_node','var') || isempty(root_node)
        root_node = 1;
    end

    num_nodes = size(skele.nodes,1);
    tree = cell(num_nodes,1);
    
    
    
    for n = 1:num_nodes
        tree{n}.loc = skele.nodes(n,:);
        tree{n}.children = [];
    end
    
    num_edges = size(skele.edges,1);
    is_valid = true(num_edges,1);
    
    parents = root_node;
    while any(is_valid) && ~isempty(parents)
        old_parents = parents;
        parents = [];
        for p = old_parents
            
            valid_edges = is_valid & (skele.edges(:,1) == p | skele.edges(:,2) == p);
            
            if any(valid_edges)
                is_valid(valid_edges) = false;

                edge_nums = full(skele.edges(valid_edges,:));
                edge_nums(edge_nums(:)==p) = [];
                tree{p}.children = edge_nums;
                parents = [parents tree{p}.children];
            
            end
        end
    end
    
end
            
    