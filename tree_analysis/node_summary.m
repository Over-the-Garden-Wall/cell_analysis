function [is_branching is_leaf] = node_summary(skele)

    edge_count = zeros(size(skele.nodes,1),1);
    for n = 1:length(skele.edges(:));
        edge_count(skele.edges(n)) = edge_count(skele.edges(n)) + 1;
    end

    is_branching = edge_count > 2;
    is_leaf = edge_count == 1;
    
end
    


    