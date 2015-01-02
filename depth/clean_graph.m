function [new_nodes new_edges] = clean_graph(nodes, edges)

    good_nodes = unique(edges(:));
    
    new_nodes = zeros(length(good_nodes),3);
    new_edges = zeros(size(edges));
    for n = 1:length(good_nodes)
        is_node = edges==good_nodes(n);
        new_edges(is_node) = n;
        new_nodes(n,:) = nodes(good_nodes(n),:);
    end
    
end